#!/usr/bin/env python
# coding: utf-8

import argparse
import os
import sys
import time
from pathlib import Path

from twisstntern.pipeline import run_analysis, ensure_twisst_available
from twisstntern.logger import setup_logging, get_logger, log_system_info, log_analysis_start, log_analysis_complete, log_error


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
        type=int,
        default=None,
        help="If set, only every Nth row of the topology weights will be used for analysis (downsampling).",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="Results",
        help="Path to the output directory. If not provided, a 'Results' directory will be created.",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging output.",
    )

    args = parser.parse_args()

    # Ensure twisst.py is available before any analysis
    ensure_twisst_available()

    # Determine which file argument to use
    input_file = args.file_flag if args.file_flag else args.file
    
    if not input_file:
        parser.error("Input file must be specified either as positional argument or with -f/--file flag")

    # Set output directory (default is "Results" if not specified)
    output_dir = args.output

    # Setup logging
    log_file_path = setup_logging(
        output_dir=output_dir,
        verbose=args.verbose,
        console_output=True
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

    # Log analysis start
    log_analysis_start(
        input_file=str(file_path),
        output_dir=output_dir,
        granularity=granularity,
        taxon_names=args.taxon_names,
        outgroup=args.outgroup,
        topology_mapping=args.topology_mapping,
        verbose=args.verbose
    )

    start_time = time.time()
    
    try:
        # Call run_analysis with the new signature that returns 3 values
        results, fundamental_results, csv_file_used = run_analysis(
            file=input_file,
            granularity=granularity,
            taxon_names=args.taxon_names,
            outgroup=args.outgroup,
            output_dir=output_dir,
            topology_mapping=args.topology_mapping,
            downsample=args.downsample,
        )

        # Calculate duration
        duration = time.time() - start_time
        
        # Collect output files
        output_files = []
        output_path = Path(output_dir)
        if output_path.exists():
            output_files.extend([str(f) for f in output_path.glob("*.*")])
        
        # Log completion
        log_analysis_complete(duration, output_files)
        
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
