#!/usr/bin/env python
# coding: utf-8

import argparse
import os
import sys
from pathlib import Path

from twisstntern.pipeline import run_analysis


def main():
    parser = argparse.ArgumentParser(
        description="Run TWISSTNTERN analysis pipeline on tree files or CSV files."
    )
    parser.add_argument(
        "file",
        type=str,
        help="Path to the input file (tree file: .trees/.ts/.newick/.tree/.nexus or CSV file: .csv).",
    )
    parser.add_argument(
        "--granularity",
        type=str,
        default="0.1",
        help="Granularity level: 'superfine', 'fine', 'coarse', or a float (e.g., 0.1). Default: 0.1",
    )
    parser.add_argument(
        "--taxon-names",
        type=str,
        nargs="+",
        help="List of taxon names for Newick tree files (e.g., --taxon-names A B C D).",
    )
    parser.add_argument(
        "--outgroup", type=str, help="Outgroup taxon name for tree files."
    )
    parser.add_argument(
        "-o", "--output", 
        type=str,
        default="Results",
        help="Path to the output directory. If not provided, a 'Results' directory will be created."
    )

    args = parser.parse_args()
    
    # Check if input file exists
    file_path = Path(args.file)
    if not file_path.exists():
        print(f"Error: Input file not found: {file_path}")
        sys.exit(1)

    # Set output directory (default is "Results" if not specified)
    output_dir = args.output


    # Convert granularity to float if possible
    try:
        granularity = float(args.granularity)
    except ValueError:
        granularity = args.granularity

    print(f"Running analysis on {args.file} with granularity {granularity}...")

    # Call run_analysis with the new signature that returns 3 values
    results, fundamental_results, csv_file_used = run_analysis(
        file=args.file,
        granularity=granularity,
        taxon_names=args.taxon_names,
        outgroup=args.outgroup,
        output_dir=output_dir,
    )

 
    print("----------------------------------------------------------")
    print("Summary of the analysis:")
    print(f"Data file used: {csv_file_used}")
    print("\nFundamental asymmetry results:")
    print(f"n_right: {fundamental_results[0]}")
    print(f"n_left: {fundamental_results[1]}")
    print(f"D_LR: {fundamental_results[2]:.4f}")
    print(f"G-test: {fundamental_results[3]:.4f}")
    print(f"p-value: {fundamental_results[4]:.4e}")
    print(f"\nResults and plots have been saved to the '{output_dir}' directory.")


if __name__ == "__main__":
    main()
