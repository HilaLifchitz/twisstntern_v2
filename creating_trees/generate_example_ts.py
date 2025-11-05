#!/usr/bin/env python3
"""Generate a small multi-population tree sequence for testing.

Example usage
-------------

    python creating_trees/generate_example_ts.py \
        --output-dir ts_data \
        --sequence-length 500000 \
        --recombination-rate 1e-8 \
        --mutation-rate 1e-8 \
        --seed 42 \
        --gzip

The script will create both a binary `.ts` file and a binary `.trees`
file (both in tskit tree sequence format). When `--gzip` is provided,
compressed `.ts.gz` and `.trees.gz` files are written as well."""

from __future__ import annotations

import argparse
import gzip
import shutil
from pathlib import Path

import msprime


def build_demography(pop_sizes: list[int], migration_rate: float) -> msprime.Demography:
    dem = msprime.Demography()
    population_names = []
    for idx, size in enumerate(pop_sizes):
        name = f"pop{idx}"
        dem.add_population(name=name, initial_size=size)
        population_names.append(name)

    for i, src in enumerate(population_names):
        for j in range(i + 1, len(population_names)):
            dst = population_names[j]
            dem.set_migration_rate(src, dst, migration_rate)
            dem.set_migration_rate(dst, src, migration_rate)

    return dem


def simulate_tree_sequence(
    sequence_length: int,
    recombination_rate: float,
    mutation_rate: float,
    seed: int,
    pop_sizes: list[int],
    samples_per_population: int,
) -> msprime.TreeSequence:
    dem = build_demography(pop_sizes, migration_rate=5e-5)

    samples = [
        msprime.SampleSet(samples_per_population, population=f"pop{i}")
        for i in range(len(pop_sizes))
    ]

    ancestry_ts = msprime.sim_ancestry(
        samples=samples,
        demography=dem,
        sequence_length=sequence_length,
        recombination_rate=recombination_rate,
        ploidy=1,
        random_seed=seed,
    )

    mutated_ts = msprime.sim_mutations(
        ancestry_ts,
        rate=mutation_rate,
        random_seed=seed + 1,
    )

    return mutated_ts


def save_outputs(
    ts: msprime.TreeSequence,
    output_dir: Path,
    stem: str,
    gzip_outputs: bool = False,
) -> tuple[Path, Path, Path | None, Path | None]:
    output_dir.mkdir(parents=True, exist_ok=True)
    ts_path = output_dir / f"{stem}.ts"
    trees_path = output_dir / f"{stem}.trees"

    ts.dump(ts_path)
    ts.dump(trees_path)

    ts_gz_path = None
    trees_gz_path = None

    if gzip_outputs:
        ts_gz_path = output_dir / f"{stem}.ts.gz"
        trees_gz_path = output_dir / f"{stem}.trees.gz"
        with open(ts_path, "rb") as src, gzip.open(ts_gz_path, "wb") as dst:
            shutil.copyfileobj(src, dst)
        with open(trees_path, "rb") as src, gzip.open(trees_gz_path, "wb") as dst:
            shutil.copyfileobj(src, dst)

    return ts_path, trees_path, ts_gz_path, trees_gz_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("ts_data"),
        help="Directory to write the outputs (default:ts_data)",
    )
    parser.add_argument(
        "--stem",
        type=str,
        default="example_chromosome",
        help="Filename stem for output files (default: example_chromosome)",
    )
    parser.add_argument(
        "--sequence-length",
        type=int,
        default=500_000,
        help="Sequence length for the simulation (default: 500000)",
    )
    parser.add_argument(
        "--recombination-rate",
        type=float,
        default=1e-8,
        help="Per-base recombination rate (default: 1e-8)",
    )
    parser.add_argument(
        "--mutation-rate",
        type=float,
        default=1e-8,
        help="Per-base mutation rate (default: 1e-8)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for the ancestry simulation (default: 42)",
    )
    parser.add_argument(
        "--samples-per-pop",
        type=int,
        default=4,
        help="Number of haploid samples per population (default: 4)",
    )
    parser.add_argument(
        "--pop-sizes",
        type=int,
        nargs="*",
        default=[10000, 8000, 6000, 5000],
        help="Initial population sizes in order 0 1 2 3 (default: 10000 8000 6000 5000)",
    )
    parser.add_argument(
        "--gzip",
        action="store_true",
        help="Also write compressed `.ts.gz` and `.trees.gz` outputs",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    ts = simulate_tree_sequence(
        sequence_length=args.sequence_length,
        recombination_rate=args.recombination_rate,
        mutation_rate=args.mutation_rate,
        seed=args.seed,
        pop_sizes=args.pop_sizes,
        samples_per_population=args.samples_per_pop,
    )

    ts_path, trees_path, ts_gz_path, trees_gz_path = save_outputs(
        ts,
        args.output_dir,
        args.stem,
        gzip_outputs=args.gzip,
    )

    print(f"Saved tree sequence (.ts):    {ts_path}")
    print(f"Saved tree sequence (.trees): {trees_path}")
    if ts_gz_path and trees_gz_path:
        print(f"Saved tree sequence (.ts.gz):   {ts_gz_path}")
        print(f"Saved tree sequence (.trees.gz): {trees_gz_path}")
    print("Populations generated: 0 (outgroup), 1, 2, 3")


if __name__ == "__main__":
    main()

# example usage
"""
python creating_trees/generate_example_ts.py \
  --output-dir ts_data \
  --stem demo_with_gzip \
  --sequence-length 500000 \
  --recombination-rate 1e-8 \
  --mutation-rate 1e-8 \
  --samples-per-pop 4 \
  --seed 42 \
  --gzip 
  """