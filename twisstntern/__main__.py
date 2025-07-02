#!/usr/bin/env python
# coding: utf-8

import os
import sys
import time
import re
from pathlib import Path
import hydra
from omegaconf import DictConfig, OmegaConf

from twisstntern.pipeline import run_analysis, ensure_twisst_available
from twisstntern.logger import setup_logging, get_logger, log_system_info, log_analysis_start, log_analysis_complete, log_error


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
    match = re.match(r'^(\d+)\+(\d+)$', downsample_str)
    if match:
        N = int(match.group(1))
        i = int(match.group(2))
        
        if N < 1:
            raise ValueError("Downsample interval N must be >= 1")
        if i >= N:
            raise ValueError(f"Starting index i ({i}) must be < N ({N})")
            
        return N, i
    
    # Invalid format
    raise ValueError(f"Invalid downsample format: '{downsample_str}'. Use 'N' or 'N+i' (e.g., '10' or '10+3')")


@hydra.main(version_base=None, config_path="../conf", config_name="config")
def main(cfg: DictConfig) -> None:
    """
    Main entry point for TWISSTNTERN analysis pipeline.
    
    Args:
        cfg: Hydra configuration object containing all parameters
    """
    # Ensure twisst.py is available before any analysis
    ensure_twisst_available()

    # Get input file from config
    input_file = cfg.file
    
    if not input_file:
        raise ValueError("Input file must be specified in configuration. Use 'file=path/to/file' in config or as command line override.")

    # Set output directory (default is "Results" if not specified)
    output_dir = cfg.output
    
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Save the final configuration to the results folder
    from omegaconf import OmegaConf
    config_save_path = Path(output_dir) / "analysis_config.yaml"
    with open(config_save_path, 'w') as f:
        OmegaConf.save(cfg, f)

    # Setup logging
    log_file_path = setup_logging(
        output_dir=output_dir,
        verbose=cfg.verbose,
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
        granularity = float(cfg.granularity)
    except (ValueError, TypeError):
        granularity = cfg.granularity

    # Parse downsample argument
    try:
        downsample_N, downsample_i = parse_downsample_arg(cfg.downsample)
    except ValueError as e:
        logger.error(f"Invalid downsample argument: {e}")
        sys.exit(1)

    # Log analysis start
    log_analysis_start(
        input_file=str(file_path),
        output_dir=output_dir,
        granularity=granularity,
        taxon_names=cfg.taxon_names,
        outgroup=cfg.outgroup,
        topology_mapping=cfg.topology_mapping,
        verbose=cfg.verbose
    )

    start_time = time.time()
    
    # Track files created during this run
    created_files = []
    
    try:
        # Call run_analysis with the new signature that returns 3 values
        results, fundamental_results, csv_file_used = run_analysis(
            file=input_file,
            granularity=granularity,
            taxon_names=cfg.taxon_names,
            outgroup=cfg.outgroup,
            output_dir=output_dir,
            topology_mapping=cfg.topology_mapping,
            downsample_N=downsample_N,
            downsample_i=downsample_i,
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
        print(f"Input file: {input_file}")
        print(f"Data file used: {csv_file_used}")
        print("\nFundamental asymmetry results:")
        print(f"n_right: {fundamental_results[0]}")
        print(f"n_left: {fundamental_results[1]}")
        print(f"D_LR: {fundamental_results[2]:.4f}")
        print(f"G-test: {fundamental_results[3]:.4f}")
        print(f"p-value: {fundamental_results[4]:.4e}")
        print(f"\nResults and plots have been saved to the '{output_dir}' directory.")
        print(f"Final configuration saved to: {config_save_path}")
        
        if log_file_path:
            print(f"Log file saved to: {log_file_path}")
            
    except Exception as e:
        log_error(e, "main analysis")
        logger.critical("Analysis failed. Check the log for details.")
        sys.exit(1)


if __name__ == "__main__":
    main()
