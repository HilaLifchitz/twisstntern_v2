#!/usr/bin/env python
"""
Main entry point for twisstntern_simulate package.

This module provides the command-line interface for running the complete
twisstntern_simulate pipeline: simulation -> twisst processing -> analysis.

Usage:
    python -m twisstntern_simulate config_file=config.yaml output=Results/
    
    # Or with additional options:
    python -m twisstntern_simulate config_file=config.yaml output=Results/ verbose=true
"""

import sys
import os
import time
import re
from pathlib import Path
import hydra
from omegaconf import DictConfig, OmegaConf

from twisstntern_simulate.pipeline import run_pipeline

# Import twisstntern logging
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


def parse_downsampleKB_arg(downsampleKB_str):
    """
    Parse downsampleKB argument in format 'Nkb' or 'Nkb+ikb' where:
    - N = sample every N kilobases
    - i = starting position (with units: kb, mb, gb)
    - If only Nkb is provided, default to Nkb+0kb (start from position 0)
    - Constraint: i < N (offset must be less than the sampling interval)
    
    Args:
        downsampleKB_str (str): String in format 'Nkb' or 'Nkb+ikb' or 'Nkb+imb' etc.
        
    Returns:
        tuple: (N, i_bp) where N is the sampling interval in kb, i_bp is starting position in base pairs
        
    Raises:
        ValueError: If format is invalid or constraint violated
    """
    if downsampleKB_str is None:
        return None, None
        
    # Handle simple case: just Nkb
    if downsampleKB_str.endswith('kb') and '+' not in downsampleKB_str:
        try:
            N = int(downsampleKB_str[:-2])  # Remove 'kb'
            if N <= 0:
                raise ValueError(f"Sampling interval must be positive, got {N}")
            return N, 0  # Default to starting from 0
        except ValueError as e:
            raise ValueError(f"Invalid downsampleKB format '{downsampleKB_str}': {e}")
    
    # Handle Nkb+ikb format
    if '+' in downsampleKB_str:
        parts = downsampleKB_str.split('+')
        if len(parts) != 2:
            raise ValueError(f"Invalid downsampleKB format '{downsampleKB_str}': expected 'Nkb+ikb'")
        
        # Parse N (sampling interval)
        N_part = parts[0].strip()
        if not N_part.endswith('kb'):
            raise ValueError(f"Invalid downsampleKB format '{downsampleKB_str}': sampling interval must end with 'kb'")
        try:
            N = int(N_part[:-2])
            if N <= 0:
                raise ValueError(f"Sampling interval must be positive, got {N}")
        except ValueError as e:
            raise ValueError(f"Invalid sampling interval '{N_part}': {e}")
        
        # Parse i (starting position)
        i_part = parts[1].strip()
        try:
            i_bp = parse_position_with_units(i_part)
        except ValueError as e:
            raise ValueError(f"Invalid starting position '{i_part}': {e}")
        
        # Convert N to base pairs for comparison
        N_bp = N * 1000
        
        # Check constraint: i < N (in base pairs)
        if i_bp >= N_bp:
            raise ValueError(f"Starting position ({i_part} = {i_bp} bp) must be less than sampling interval ({N}kb = {N_bp} bp)")
        
        return N, i_bp
    
    raise ValueError(f"Invalid downsampleKB format '{downsampleKB_str}': expected 'Nkb' or 'Nkb+ikb'")


def parse_position_with_units(position_str):
    """
    Parse a position string with units (kb, mb, gb) and return base pairs.
    
    Args:
        position_str (str): Position string like "50kb", "30mb", "1gb"
        
    Returns:
        int: Position in base pairs
        
    Raises:
        ValueError: If format is invalid
    """
    position_str = position_str.strip().lower()
    
    # Handle different units
    if position_str.endswith('kb'):
        try:
            value = float(position_str[:-2])
            return int(value * 1000)
        except ValueError:
            raise ValueError(f"Invalid kb value: {position_str}")
    elif position_str.endswith('mb'):
        try:
            value = float(position_str[:-2])
            return int(value * 1000000)
        except ValueError:
            raise ValueError(f"Invalid mb value: {position_str}")
    elif position_str.endswith('gb'):
        try:
            value = float(position_str[:-2])
            return int(value * 1000000000)
        except ValueError:
            raise ValueError(f"Invalid gb value: {position_str}")
    elif position_str.isdigit():
        # Assume base pairs if no unit specified
        return int(position_str)
    else:
        raise ValueError(f"Invalid position format '{position_str}': expected 'Nkb', 'Nmb', 'Ngb', or 'N' (base pairs)")


@hydra.main(version_base=None, config_path="../conf", config_name="simulate_config")
def main(cfg: DictConfig) -> None:
    """
    Main entry point for the twisstntern_simulate package.

    Uses Hydra configuration and runs the full pipeline:
    1. Load configuration from YAML file
    2. Run msprime simulation
    3. Ensure twisst is available (download if needed)
    4. Process tree sequences to generate topology weights
    5. Run twisstntern analysis and generate plots
    
    Args:
        cfg: Hydra configuration object containing all parameters
    """

    # Create output directory if it doesn't exist
    output_dir = Path(cfg.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save the final configuration to the results folder
    from omegaconf import OmegaConf
    config_save_path = output_dir / "simulation_config.yaml"
    with open(config_save_path, 'w') as f:
        OmegaConf.save(cfg, f)

    # Setup logging (like twisstntern)
    log_file_path = setup_logging(
        output_dir=str(output_dir),
        verbose=cfg.verbose,
        console_output=not cfg.quiet
    )
    
    logger = get_logger(__name__)

    # Log system information
    log_system_info()

    # Handle configuration - use embedded config if no external file specified
    if cfg.config_file is not None:
        # Use external configuration file
        config_path = Path(cfg.config_file)
        if not config_path.exists():
            logger.error(f"Configuration file not found: {config_path}")
            sys.exit(1)
        logger.info(f"Using external simulation configuration: {config_path}")
        use_embedded_config = False
    else:
        # Use embedded configuration from Hydra config
        logger.info("Using embedded simulation configuration from Hydra config")
        config_path = None
        use_embedded_config = True

    # Convert granularity to float if needed
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

    # Parse downsampleKB argument
    try:
        downsampleKB_N, downsampleKB_i = parse_downsampleKB_arg(cfg.downsampleKB)
    except ValueError as e:
        logger.error(f"Invalid downsampleKB argument: {e}")
        sys.exit(1)

    # Log analysis start
    config_description = "embedded Hydra config" if use_embedded_config else str(config_path)
    log_analysis_start(
        input_file=config_description,
        output_dir=str(output_dir),
        granularity=granularity,
        seed_override=cfg.seed,
        mode_override=None,
        topology_mapping=cfg.topology_mapping,
        verbose=cfg.verbose
    )

    start_time = time.time()
    
    # Track files created during this run
    created_files = []

    try:
        # Run the pipeline
        if use_embedded_config:
            # Save embedded config to a temporary file for the pipeline
            temp_config_path = output_dir / "temp_simulation_config.yaml"
            with open(temp_config_path, 'w') as f:
                OmegaConf.save(cfg.simulation, f)
            config_path_for_pipeline = str(temp_config_path)
        else:
            config_path_for_pipeline = str(config_path)
            
        results = run_pipeline(
            config_path=config_path_for_pipeline,
            output_dir=str(output_dir),
            seed_override=cfg.seed,
            mode_override=None,
            granularity=granularity,
            verbose=cfg.verbose,
            topology_mapping=cfg.topology_mapping,
            config_overrides=cfg.override,
            downsample_N=downsample_N,
            downsample_i=downsample_i,
            downsample_kb=downsampleKB_N,
            downsample_kb_i=downsampleKB_i,
        )

        # Calculate duration
        duration = time.time() - start_time
        
        # Collect files created during this run
        if output_dir.exists():
            # Get the timestamp when we started (approximate)
            start_timestamp = start_time - 1  # Subtract 1 second to be safe
            
            # Only include files created during this run
            for file_path in output_dir.glob("*.*"):
                if file_path.stat().st_mtime >= start_timestamp:
                    created_files.append(str(file_path))
        
        # Clean up temporary config file if used
        if use_embedded_config:
            temp_config_path = output_dir / "temp_simulation_config.yaml"
            if temp_config_path.exists():
                temp_config_path.unlink()
        
        # Log completion with only the files created in this run
        log_analysis_complete(duration, created_files)
        
        # Print summary to console (like main twisstntern)
        print("----------------------------------------------------------")
        print("Summary of the analysis:")
        print(f"Configuration used: {config_description}")
        print(f"Data file used: {results['csv_file_used']}")
        print(f"\nResults and plots have been saved to the '{output_dir}' directory.")
        print(f"Final configuration saved to: {config_save_path}")
        
        if log_file_path:
            print(f"Log file saved to: {log_file_path}")

    except Exception as e:
        log_error(e, "twisstntern_simulate main analysis")
        logger.critical("Analysis failed. Check the log for details.")
        sys.exit(1)


if __name__ == "__main__":
    main()
