#!/usr/bin/env python3
"""Convert a multi-line Newick file into a NEXUS file.

Each non-empty line in the input is treated as a separate tree and written
into a basic NEXUS structure under a TAXA block.

Example usage
-------------

    python creating_trees/convert_newick_to_nexus.py \
        --input trees.newick \
        --output trees.nexus

    python creating_trees/convert_newick_to_nexus.py \
        --input trees.newick \
        --output trees.nexus.gz \
        --gzip
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import Iterable


def iter_newick_lines(path: Path) -> Iterable[str]:
    with open(path, "r") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped:
                yield stripped


def write_nexus(trees: Iterable[str], output_path: Path, compress: bool = False) -> None:
    opener = gzip.open if compress else open
    with opener(output_path, "wt") as handle:
        handle.write("#NEXUS\n")
        handle.write("Begin trees;\n")
        for idx, tree in enumerate(trees, start=1):
            handle.write(f"  Tree tree{idx} = {tree}\n")
        handle.write("End;\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Source multi-line Newick file")
    parser.add_argument("--output", type=Path, required=True, help="Destination NEXUS file")
    parser.add_argument("--gzip", action="store_true", help="Compress the output with gzip")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    trees = list(iter_newick_lines(args.input))
    write_nexus(trees, args.output, compress=args.gzip)
    print(f"Converted {len(trees)} trees to {args.output}")


if __name__ == "__main__":
    main()
