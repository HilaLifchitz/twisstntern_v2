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
import logging
from pathlib import Path

from twisstntern_simulate.pipeline import run_pipeline

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


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

    parser.add_argument(
        "--mode",
        choices=["locus", "chromosome"],
        help="Simulation mode (overrides config file).",
    )

    # Analysis options
    parser.add_argument(
        "--granularity",
        type=float,
        default=0.1,
        help="Granularity for ternary analysis (default: 0.1).",
    )


    args = parser.parse_args()

    # Configure logging based on verbosity
    if args.quiet:
        logging.getLogger().setLevel(logging.ERROR)
    elif args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Add file logging if requested
    if args.log_file:
        file_handler = logging.FileHandler(args.log_file)
        file_handler.setFormatter(
            logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
        )
        logging.getLogger().addHandler(file_handler)

    # Validate arguments - no conflicts to check now

    # Check if config file exists
    config_path = Path(args.config)
    if not config_path.exists():
        logger.error(f"Configuration file not found: {config_path}")
        sys.exit(1)

    # Create output directory if it doesn't exist
    # This directory will contain: topology weights, analysis results, plots, CSV files, and tree files
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Starting twisstntern_simulate pipeline")
    logger.info(f"Config file: {config_path}")
    logger.info(f"Output directory: {output_dir}")

    try:
        # Run the pipeline
        run_pipeline(
            config_path=str(config_path),
            output_dir=str(output_dir),
            skip_twisst_check=args.skip_twisst_check,
            force_download=args.force_download,
            seed_override=args.seed,
            mode_override=args.mode,
            granularity=args.granularity,
            verbose=args.verbose,
        )

        logger.info("Pipeline completed successfully!")

    except Exception as e:
        logger.error(f"Pipeline failed with error: {str(e)}")
        if args.verbose:
            import traceback

            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
