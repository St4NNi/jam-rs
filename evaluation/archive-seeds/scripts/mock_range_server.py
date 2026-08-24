#!/usr/bin/env python3
"""Serve local objects through bounded HTTP byte ranges for remote tests.

The server is intentionally small and dependency-free.  It implements HEAD,
GET, and a single ``Range: bytes=start-end`` request, records request and byte
counters, and never follows paths outside the supplied workspace root.  The
stats file is write-once so a run cannot silently mix observations.
"""

from __future__ import annotations

import argparse
import json
import signal
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlsplit

from common import ExperimentError, ensure_new_output, workspace_path


class _Counters:
    def __init__(self) -> None:
        self.lock = threading.Lock()
        self.head_requests = 0
        self.get_requests = 0
        self.range_requests = 0
        self.full_object_requests = 0
        self.bytes_requested = 0
        self.bytes_returned = 0

    def add(self, *, method: str, length: int, ranged: bool) -> None:
        with self.lock:
            if method == "HEAD":
                self.head_requests += 1
            else:
                self.get_requests += 1
            if ranged:
                self.range_requests += 1
            else:
                self.full_object_requests += 1
            self.bytes_requested += max(0, length)
            self.bytes_returned += max(0, length)

    def snapshot(self) -> dict[str, int]:
        with self.lock:
            return {
                "head_requests": self.head_requests,
                "get_requests": self.get_requests,
                "range_requests": self.range_requests,
                "full_object_requests": self.full_object_requests,
                "bytes_requested": self.bytes_requested,
                "bytes_returned": self.bytes_returned,
            }


class _Server(ThreadingHTTPServer):
    def __init__(self, address: tuple[str, int], root: Path, counters: _Counters, max_requests: int) -> None:
        super().__init__(address, _Handler)
        self.root = root
        self.counters = counters
        self.max_requests = max_requests
        self.request_count = 0
        self.request_lock = threading.Lock()

    def count_request(self) -> bool:
        with self.request_lock:
            self.request_count += 1
            return self.max_requests > 0 and self.request_count >= self.max_requests


class _Handler(BaseHTTPRequestHandler):
    server: _Server
    protocol_version = "HTTP/1.1"

    def log_message(self, format: str, *args: object) -> None:
        # Measurements should be machine-readable; request details are in the
        # stats object and not an unbounded stderr log.
        return

    def _file(self) -> Path | None:
        name = unquote(urlsplit(self.path).path).lstrip("/")
        if not name:
            return None
        candidate = workspace_path(self.server.root / name, field="served path")
        if not candidate.is_file():
            return None
        return candidate

    def _range(self, size: int) -> tuple[int, int] | None:
        value = self.headers.get("Range")
        if value is None:
            return None
        prefix, separator, span = value.partition("=")
        if prefix.strip().lower() != "bytes" or separator != "=" or "," in span:
            return None
        start_text, separator, end_text = span.strip().partition("-")
        try:
            if not start_text:
                length = int(end_text)
                if length <= 0:
                    return None
                start = max(0, size - length)
                end = size - 1
            else:
                start = int(start_text)
                end = int(end_text) if end_text else size - 1
        except ValueError:
            return None
        if start < 0 or start >= size or end < start:
            return None
        return start, min(end, size - 1)

    def _serve(self, *, head: bool) -> None:
        path = self._file()
        if path is None:
            self.send_error(404)
            return
        size = path.stat().st_size
        selected = self._range(size)
        if selected is None and self.headers.get("Range") is not None:
            self.send_response(416)
            self.send_header("Content-Range", f"bytes */{size}")
            self.send_header("Content-Length", "0")
            self.end_headers()
            return
        if selected is None:
            start, end = 0, max(0, size - 1)
            status = 200
        else:
            start, end = selected
            status = 206
        length = max(0, end - start + 1)
        self.send_response(status)
        self.send_header("Accept-Ranges", "bytes")
        self.send_header("Content-Length", str(length))
        if selected is not None:
            self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
        self.end_headers()
        if head:
            self.server.counters.add(method="HEAD", length=0, ranged=selected is not None)
            if self.server.count_request():
                threading.Thread(target=self.server.shutdown, daemon=True).start()
            return
        with path.open("rb") as handle:
            handle.seek(start)
            payload = handle.read(length)
        self.server.counters.add(method="GET", length=length, ranged=selected is not None)
        should_stop = self.server.count_request()
        self.wfile.write(payload)
        if should_stop:
            threading.Thread(target=self.server.shutdown, daemon=True).start()

    def do_HEAD(self) -> None:  # noqa: N802 - stdlib handler API
        self._serve(head=True)

    def do_GET(self) -> None:  # noqa: N802 - stdlib handler API
        self._serve(head=False)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path, help="workspace directory containing objects")
    parser.add_argument("--bind", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=0)
    parser.add_argument("--stats", type=Path, required=True, help="new JSON stats output")
    parser.add_argument("--max-requests", type=int, default=0, help="stop after this many requests; zero means serve until interrupted")
    args = parser.parse_args()
    try:
        root = workspace_path(args.root, field="root")
        stats = ensure_new_output(args.stats)
        counters = _Counters()
        server = _Server((args.bind, args.port), root, counters, args.max_requests)
    except (ExperimentError, OSError, ValueError) as exc:
        parser.exit(1, f"mock range server failed: {exc}\n")

    def stop(_signum: int, _frame: object) -> None:
        threading.Thread(target=server.shutdown, daemon=True).start()

    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    print(json.dumps({"status": "ready", "host": args.bind, "port": server.server_port}), flush=True)
    try:
        server.serve_forever(poll_interval=0.1)
    finally:
        server.server_close()
        value = {
            "schema_version": "archive-seeds-range-server-v1",
            "root": str(root),
            "address": {"host": args.bind, "port": server.server_port},
            "metrics": counters.snapshot(),
        }
        stats.parent.mkdir(parents=True, exist_ok=True)
        try:
            with stats.open("x", encoding="utf-8") as handle:
                json.dump(value, handle, indent=2, sort_keys=True)
                handle.write("\n")
        except FileExistsError as exc:
            parser.exit(1, f"refusing to overwrite stats: {stats}: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
