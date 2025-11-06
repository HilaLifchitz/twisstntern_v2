#!/usr/bin/env python3
"""Compress a multi-line Newick file to `.newick.gz`.

Usage examples
--------------

    python creating_trees/convert_newick_to_gz.py \
        --input trees.newick \
        --output trees.newick.gz

    python creating_trees/convert_newick_to_gz.py \
        --input trees.newick \
        --output trees.newick.gz \
        --force

Each non-empty line is copied verbatim; the script just wraps the file in
gzip compression."""

from __future__ import annotations

import argparse
import gzip
import shutil
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Source Newick file")
    parser.add_argument("--output", type=Path, required=True, help="Destination .newick.gz file")
    parser.add_argument("--force", action="store_true", help="Overwrite output if it exists")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if args.output.exists() and not args.force:
        raise FileExistsError(f"Output file {args.output} already exists; use --force to overwrite")

    with open(args.input, "rb") as src, gzip.open(args.output, "wb") as dst:
        shutil.copyfileobj(src, dst)

    print(f"Compressed {args.input} to {args.output}")


if __name__ == "__main__":
    main()