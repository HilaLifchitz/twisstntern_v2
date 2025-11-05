#!/usr/bin/env python3
"""Convert a multi-line Newick file to a `.tree` (Relate-style) file.

Usage examples
--------------

Convert to a plain `.tree` file:

    python creating_trees/convert_newick_to_tree.py input.newick output.tree

Convert and gzip the result:

    python creating_trees/convert_newick_to_tree.py input.newick output.tree --gzip

Each non-empty line from the input is written verbatim to the output. The
resulting file can be consumed by tools that expect one Newick tree per line.
"""

from __future__ import annotations

import argparse
import gzip
import pathlib
from typing import Iterable


def iter_newick_lines(path: pathlib.Path) -> Iterable[str]:
    with open(path, "r") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped:
                yield stripped


def write_tree_file(lines: Iterable[str], output_path: pathlib.Path, compress: bool) -> None:
    opener = gzip.open if compress else open
    mode = "wt"
    with opener(output_path, mode) as out_handle:
        for line in lines:
            out_handle.write(f"{line}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=pathlib.Path, help="Path to the Newick file")
    parser.add_argument("output", type=pathlib.Path, help="Destination `.tree` file")
    parser.add_argument(
        "--gzip",
        action="store_true",
        help="Compress the output with gzip (writes `.gz` file)",
    )

    args = parser.parse_args()

    lines = list(iter_newick_lines(args.input))
    write_tree_file(lines, args.output, args.gzip)


if __name__ == "__main__":
    main()