#!/usr/bin/env python
# coding: utf-8

import argparse
import os
import sys
import time
import re
from pathlib import Path

from .pipeline import run_analysis
from .logger import (
    setup_logging,
    get_logger,
    log_system_info,
    log_analysis_start,
    log_analysis_complete,
    log_error,
)


def parse_downsample_arg(downsample_str):
    """
    Parse downsample argument in format 'N' or 'N+i' where:
    - N = sample every Nth tree/locus
    - i = starting index (offset)
    - If only N is provided, default to N+0 (start from index 0)
    - Constraint: i < N (offset must be less than the sampling interval)

    Args:
        downsample_str (str): String in format 'N' or 'N+i'

    Returns:
        tuple: (N, i) where N is the sampling interval and i is the starting index

    Raises:
        ValueError: If format is invalid or i >= N
    """
    if downsample_str is None:
        return None, None

    # Check if it's just a number (N format)
    if downsample_str.isdigit():
        N = int(downsample_str)
        if N < 1:
            raise ValueError("Downsample interval N must be >= 1")
        return N, 0  # Default to starting from index 0

    # Check if it's in N+i format
    match = re.match(r"^(\d+)\+(\d+)$", downsample_str)
    if match:
        N = int(match.group(1))
        i = int(match.group(2))

        if N < 1:
            raise ValueError("Downsample interval N must be >= 1")
        if i >= N:
            raise ValueError(f"Starting index i ({i}) must be < N ({N})")

        return N, i

    # Invalid format
    raise ValueError(
        f"Invalid downsample format: '{downsample_str}'. Use 'N' or 'N+i' (e.g., '10' or '10+3')"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Run TWISSTNTERN analysis pipeline on tree files or CSV files."
    )
    parser.add_argument(
        "file",
        type=str,
        nargs="?",
        help="Path to the input file (tree file: .trees/.ts/.newick/.tree/.nexus or CSV file: .csv).",
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        dest="file_flag",
        help="Path to the input file (alternative to positional argument).",
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
        "--topology-mapping",
        type=str,
        help='Custom topology ordering. Format: \'T1="(0,(3,(1,2)))"; T2="(0,(1,(2,3)))"; T3="(0,(2,(1,3)))";\'',
    )
    parser.add_argument(
        "--downsample",
        type=str,
        default=None,
        help="Downsample format: 'N' or 'N+i' where N=sample every Nth row, i=starting index. "
        "Examples: '10' (every 10th starting from 0), '10+1' (every 10th starting from index 1), "
        "'5+3' (every 5th starting from index 3). Constraint: i < N.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="Results",
        help="Path to the output directory for all results and plots. If not provided, a timestamped 'Results_YYYY-MM-DD_HH-MM-SS' directory will be created.",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging output.",
    )
    parser.add_argument(
        "--axis",
        type=str,
        nargs=3,
        metavar=("AXIS1", "AXIS2", "AXIS3"),
        default=["T1", "T2", "T3"],
        help="Order of axes for CSV columns. Example: --axis T2 T1 T3 (default: T1 T2 T3)",
    )

    args = parser.parse_args()

    # Determine which file argument to use
    input_file = args.file_flag if args.file_flag else args.file

    if not input_file:
        parser.error(
            "Input file must be specified either as positional argument or with -f/--file flag"
        )

    # Set output directory with timestamp if using default
    if args.output == "Results":
        # Generate timestamped directory name
        timestamp = time.strftime("%Y-%m-%d_%H-%M-%S")
        output_dir = f"Results_{timestamp}"
    else:
        # Use user-specified directory name
        output_dir = args.output

    # Setup logging
    log_file_path = setup_logging(
        output_dir=output_dir, verbose=args.verbose, console_output=True
    )

    logger = get_logger(__name__)

    # Log system information
    log_system_info()

    # Check if input file exists
    file_path = Path(input_file)
    if not file_path.exists():
        logger.error(f"Input file not found: {file_path}")
        sys.exit(1)

    # Convert granularity to float if possible
    try:
        granularity = float(args.granularity)
    except ValueError:
        granularity = args.granularity

    # Parse downsample argument
    try:
        downsample_N, downsample_i = parse_downsample_arg(args.downsample)
    except ValueError as e:
        logger.error(f"Invalid downsample argument: {e}")
        sys.exit(1)

    # Log analysis start
    log_analysis_start(
        input_file=str(file_path),
        output_dir=output_dir,
        granularity=granularity,
        taxon_names=args.taxon_names,
        outgroup=args.outgroup,
        topology_mapping=args.topology_mapping,
        verbose=args.verbose,
    )

    start_time = time.time()

    # Track files created during this run
    created_files = []

    try:
        # Call run_analysis with the new signature that returns 3 values
        results, fundamental_results, csv_file_used = run_analysis(
            file=str(file_path),
            granularity=granularity,
            taxon_names=args.taxon_names,
            outgroup=args.outgroup,
            output_dir=output_dir,
            topology_mapping=args.topology_mapping,
            downsample_N=downsample_N,
            downsample_i=downsample_i,
            colormap="viridis_r",
            axis_order=args.axis,  # <--- ADD THIS LINE
        )

        # Calculate duration
        duration = time.time() - start_time

        # Collect files created during this run
        output_path = Path(output_dir)
        if output_path.exists():
            # Get the timestamp when we started (approximate)
            start_timestamp = start_time - 1  # Subtract 1 second to be safe

            # Only include files created during this run
            for file_path in output_path.glob("*.*"):
                if file_path.stat().st_mtime >= start_timestamp:
                    created_files.append(str(file_path))

        # Log completion with only the files created in this run
        log_analysis_complete(duration, created_files)

        # Print summary to console
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

        if log_file_path:
            print(f"Log file saved to: {log_file_path}")

    except Exception as e:
        log_error(e, "main analysis")
        logger.critical("Analysis failed. Check the log for details.")
        sys.exit(1)


if __name__ == "__main__":
    main()
