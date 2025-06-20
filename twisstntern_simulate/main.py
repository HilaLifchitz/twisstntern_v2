#!/usr/bin/env python
"""
Main entry point for twisstntern_simulate package.

This module provides the command-line interface for running the complete
twisstntern_simulate pipeline: simulation -> twisst processing -> analysis.

Usage:
    python -m twisstntern_simulate -c config.yaml -o Results/
    
    # Or with additional options:
    python -m twisstntern_simulate -c config.yaml -o Results/ --force-download --verbose
"""

import argparse
import sys
import os
import time
from pathlib import Path

from twisstntern_simulate.pipeline import run_pipeline

# Import twisstntern logging
from twisstntern.logger import setup_logging, get_logger, log_system_info, log_analysis_start, log_analysis_complete, log_error


def main():
    """
    Main entry point for the twisstntern_simulate package.

    Parses command-line arguments and runs the full pipeline:
    1. Load configuration from YAML file
    2. Run msprime simulation
    3. Ensure twisst is available (download if needed)
    4. Process tree sequences to generate topology weights
    5. Run twisstntern analysis and generate plots
    """

    parser = argparse.ArgumentParser(
        description="Run the complete twisstntern_simulate pipeline: simulation -> processing -> analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with config file
  python -m twisstntern_simulate -c config.yaml -o Results/
  
  # With verbose output
  python -m twisstntern_simulate -c config.yaml -o Results/ --verbose
  
  # Using default Results directory
  python -m twisstntern_simulate -c config.yaml
        """,
    )

    # Required arguments
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        required=True,
        help="Path to configuration file (YAML format). See config_template.yaml for example.",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        required=False,
        default="Results",
        help="Output directory for results. Defaults to 'Results' if not specified. Will be created if it doesn't exist.",
    )

    parser.add_argument(
        "--skip-twisst-check",
        action="store_true",
        help="Skip checking/downloading twisst (assume it's already available).",
    )

    # Twisst-related options
    parser.add_argument(
        "--force-download",
        action="store_true",
        help="Force re-download of twisst even if it already exists.",
    )

    # Output and logging options
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Enable verbose output."
    )

    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Suppress most output (only errors will be shown).",
    )

    parser.add_argument(
        "--log-file",
        type=str,
        help="Path to log file. If not specified, logs only to console.",
    )

    # Advanced simulation overrides (optional)
    parser.add_argument(
        "--seed", type=int, help="Random seed for simulation (overrides config file)."
    )



    # Analysis options
    parser.add_argument(
        "--granularity",
        type=float,
        default=0.1,
        help="Granularity for ternary analysis (default: 0.1).",
    )

    parser.add_argument(
        "--topology-mapping",
        type=str,
        help="Custom topology mapping for T1/T2/T3. "
             "Format: 'T1=(0,(1,(2,3))); T2=(0,(2,(1,3))); T3=(0,(3,(1,2)));' "
             "This allows you to specify which topology should be assigned to each axis.",
    )

    # Configuration override options
    parser.add_argument(
        "--override",
        action="append",
        help="Override specific config values. Format: 'key=value' or 'nested.key=value'. "
             "Examples: --override 'migration.p3>p2=0.3' --override 'ploidy=2' --override 'populations.p1.Ne=15000'",
    )


    args = parser.parse_args()

    # Create output directory if it doesn't exist
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Setup logging (like twisstntern)
    log_file_path = setup_logging(
        output_dir=str(output_dir),
        verbose=args.verbose,
        console_output=not args.quiet
    )
    
    logger = get_logger(__name__)

    # Log system information
    log_system_info()

    # Check if config file exists
    config_path = Path(args.config)
    if not config_path.exists():
        logger.error(f"Configuration file not found: {config_path}")
        sys.exit(1)

    # Convert granularity to float if needed
    try:
        granularity = float(args.granularity)
    except ValueError:
        granularity = args.granularity

    # Log analysis start
    log_analysis_start(
        input_file=str(config_path),
        output_dir=str(output_dir),
        granularity=granularity,
        seed_override=args.seed,
        mode_override=None,
        topology_mapping=args.topology_mapping,
        verbose=args.verbose
    )

    start_time = time.time()

    try:
        # Run the pipeline
        results = run_pipeline(
            config_path=str(config_path),
            output_dir=str(output_dir),
            skip_twisst_check=args.skip_twisst_check,
            force_download=args.force_download,
            seed_override=args.seed,
            mode_override=None,
            granularity=granularity,
            verbose=args.verbose,
            topology_mapping=args.topology_mapping,
            config_overrides=args.override,
        )

        # Calculate duration
        duration = time.time() - start_time
        
        # Collect output files
        output_files = []
        if output_dir.exists():
            output_files.extend([str(f) for f in output_dir.glob("*.*")])
        
        # Log completion
        log_analysis_complete(duration, output_files)
        
        # Print summary to console (like main twisstntern)
        print("----------------------------------------------------------")
        print("Summary of the analysis:")
        print(f"Data file used: {results['csv_file_used']}")
        print("\nFundamental asymmetry results:")
        fundamental_results = results["fundamental_results"]
        print(f"n_right: {fundamental_results[0]}")
        print(f"n_left: {fundamental_results[1]}")
        print(f"D_LR: {fundamental_results[2]:.4f}")
        print(f"G-test: {fundamental_results[3]:.4f}")
        print(f"p-value: {fundamental_results[4]:.4e}")
        print(f"\nResults and plots have been saved to the '{output_dir}' directory.")
        
        if log_file_path:
            print(f"Log file saved to: {log_file_path}")

    except Exception as e:
        log_error(e, "twisstntern_simulate main analysis")
        logger.critical("Analysis failed. Check the log for details.")
        sys.exit(1)


if __name__ == "__main__":
    main()
