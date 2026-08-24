#!/usr/bin/env python3
"""Small range-capable HTTP fixture with request and byte counters.

It is a benchmark fixture, not a production object server.  The server binds
only to the requested host and publishes counters through ``/__metrics``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import signal
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from threading import Lock
from urllib.parse import unquote, urlsplit


class Metrics:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.lock = Lock()
        self.values = {
            "schema_version": "1.0.0",
            "requests": 0,
            "head_requests": 0,
            "get_requests": 0,
            "range_requests": 0,
            "bytes_served": 0,
            "responses_206": 0,
            "responses_200": 0,
            "status_404": 0,
            "status_416": 0,
        }

    def update(self, **values: int) -> None:
        with self.lock:
            for key, value in values.items():
                self.values[key] = self.values.get(key, 0) + value
            self.path.write_text(json.dumps(self.values, indent=2) + "\n", encoding="utf-8")

    def snapshot(self) -> dict:
        with self.lock:
            return dict(self.values)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


class Handler(BaseHTTPRequestHandler):
    server_version = "jam-trace-range-fixture/1.0"

    def _root(self) -> Path:
        return self.server.root  # type: ignore[attr-defined,no-any-return]

    def _metrics(self) -> Metrics:
        return self.server.metrics  # type: ignore[attr-defined,no-any-return]

    def _path(self) -> Path | None:
        relative = unquote(urlsplit(self.path).path).lstrip("/")
        candidate = (self._root() / relative).resolve()
        try:
            candidate.relative_to(self._root().resolve())
        except ValueError:
            return None
        return candidate if candidate.is_file() else None

    def do_HEAD(self) -> None:  # noqa: N802
        self._serve(head=True)

    def do_GET(self) -> None:  # noqa: N802
        if urlsplit(self.path).path == "/__metrics":
            body = json.dumps(self._metrics().snapshot(), sort_keys=True).encode("utf-8") + b"\n"
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)
            return
        self._serve(head=False)

    def _serve(self, head: bool) -> None:
        self._metrics().update(requests=1, head_requests=int(head), get_requests=int(not head))
        path = self._path()
        if path is None:
            self._metrics().update(status_404=1)
            self.send_error(404)
            return
        size = path.stat().st_size
        etag = self.server.etags[path.resolve()]  # type: ignore[attr-defined]
        range_header = self.headers.get("Range")
        start, end = 0, size
        partial = False
        if range_header:
            try:
                unit, value = range_header.split("=", 1)
                if unit != "bytes" or "," in value:
                    raise ValueError
                first, last = value.split("-", 1)
                if first:
                    start = int(first)
                    end = int(last) + 1 if last else size
                else:
                    suffix = int(last)
                    start = max(0, size - suffix)
                    end = size
                if start < 0 or start >= end or end > size:
                    raise ValueError
                partial = True
            except (ValueError, TypeError):
                self._metrics().update(status_416=1)
                self.send_response(416)
                self.send_header("Content-Range", f"bytes */{size}")
                self.end_headers()
                return
        body_length = end - start
        self.send_response(206 if partial else 200)
        self.send_header("Accept-Ranges", "bytes")
        self.send_header("Content-Length", str(body_length))
        self.send_header("Content-Type", "application/octet-stream")
        self.send_header("ETag", etag)
        if partial:
            self.send_header("Content-Range", f"bytes {start}-{end - 1}/{size}")
            self._metrics().update(range_requests=1, responses_206=1)
        else:
            self._metrics().update(responses_200=1)
        self.end_headers()
        if not head:
            with path.open("rb") as handle:
                handle.seek(start)
                remaining = body_length
                while remaining:
                    block = handle.read(min(1024 * 1024, remaining))
                    if not block:
                        break
                    self.wfile.write(block)
                    remaining -= len(block)
            self._metrics().update(bytes_served=body_length)

    def log_message(self, _format: str, *_args: object) -> None:
        return


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--metrics", type=Path, required=True)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=0)
    args = parser.parse_args()
    if not args.root.is_dir():
        raise SystemExit(f"HTTP fixture root is not a directory: {args.root}")
    args.metrics.parent.mkdir(parents=True, exist_ok=True)
    metrics = Metrics(args.metrics)
    server = ThreadingHTTPServer((args.host, args.port), Handler)
    server.root = args.root.resolve()  # type: ignore[attr-defined]
    server.metrics = metrics  # type: ignore[attr-defined]
    server.etags = {  # type: ignore[attr-defined]
        path.resolve(): sha256_file(path)[:16]
        for path in server.root.iterdir()  # type: ignore[attr-defined]
        if path.is_file()
    }

    def stop(_signum: int, _frame: object) -> None:
        threading.Thread(target=server.shutdown, daemon=True).start()

    signal.signal(signal.SIGTERM, stop)
    signal.signal(signal.SIGINT, stop)
    print(f"PORT={server.server_address[1]}", flush=True)
    try:
        server.serve_forever(poll_interval=0.1)
    finally:
        server.server_close()
        args.metrics.write_text(json.dumps(metrics.snapshot(), indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
