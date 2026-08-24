#!/usr/bin/env python3
"""Create standard FAI sidecars for plain FASTA benchmark inputs."""

from __future__ import annotations

import argparse
from pathlib import Path


def index_fasta(path: Path, force: bool) -> Path:
    output = Path(f"{path}.fai")
    if output.exists() and not force:
        raise ValueError(f"refusing to replace existing FAI: {output}")
    records: list[tuple[str, int, int, int, int]] = []
    name: str | None = None
    length = 0
    sequence_offset = 0
    line_bases = 0
    line_width = 0
    short_line = False

    def finish() -> None:
        if name is None:
            return
        if length == 0 or line_bases == 0 or line_width < line_bases:
            raise ValueError(f"empty or invalid FASTA record {name!r} in {path}")
        records.append((name, length, sequence_offset, line_bases, line_width))

    with path.open("rb") as handle:
        while True:
            raw = handle.readline()
            if not raw:
                break
            if raw.startswith(b">"):
                finish()
                identifier = raw[1:].strip().split(None, 1)[0]
                if not identifier:
                    raise ValueError(f"empty FASTA identifier in {path}")
                name = identifier.decode("utf-8")
                length = 0
                sequence_offset = handle.tell()
                line_bases = 0
                line_width = 0
                short_line = False
                continue
            if name is None:
                raise ValueError(f"sequence precedes FASTA identifier in {path}")
            bases = raw.rstrip(b"\r\n")
            if not bases:
                raise ValueError(f"blank sequence line in {path}")
            if line_bases == 0:
                line_bases = len(bases)
                line_width = len(raw)
            elif short_line:
                raise ValueError(f"short FASTA line is not final for {name!r} in {path}")
            elif len(bases) < line_bases:
                short_line = True
            elif len(bases) != line_bases or len(raw) != line_width:
                raise ValueError(f"variable FASTA line width for {name!r} in {path}")
            length += len(bases)
    finish()
    if not records:
        raise ValueError(f"FASTA has no records: {path}")
    text = "".join(
        f"{record[0]}\t{record[1]}\t{record[2]}\t{record[3]}\t{record[4]}\n"
        for record in records
    )
    output.write_text(text, encoding="utf-8")
    return output


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", type=Path, nargs="+")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    for path in args.fasta:
        print(index_fasta(path, args.force))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
