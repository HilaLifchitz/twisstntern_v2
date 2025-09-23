"""
Main entry point for twisstntern_simulate with Hydra configuration.

This provides the same functionality as the original argparse-based system
but uses Hydra for configuration management.
"""

import time
import sys
import re
from pathlib import Path
from typing import Optional, Tuple

import hydra
from omegaconf import DictConfig
from hydra.core.config_store import ConfigStore

from .hydra_config import TwisstnternSimulateConfig
from .pipeline import run_pipeline

# Import twisstntern logging
from ..core.logger import setup_logging, get_logger, log_system_info, log_analysis_start, log_analysis_complete, log_error


# Register the config schema with Hydra
cs = ConfigStore.instance()
cs.store(name="base_config", node=TwisstnternSimulateConfig)


def parse_downsample_arg(downsample_str: Optional[str]) -> Tuple[Optional[int], Optional[int]]:
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


def parse_position_with_units(position_str: str) -> int:
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


def parse_downsampleKB_arg(downsampleKB_str: Optional[str]) -> Tuple[Optional[int], Optional[int]]:
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


def handle_get_config() -> None:
    """Handle --get-config command for downloading template configuration."""
    from .utils import download_config_template
    
    print("🔧 TWISSTNTERN_SIMULATE Configuration Template")
    print("=" * 50)
    
    if len(sys.argv) > 2:
        # Custom destination path provided
        destination = sys.argv[2]
        downloaded_path = download_config_template(destination)
    else:
        # Default to current directory
        downloaded_path = download_config_template()
    
    if downloaded_path:
        print(f"📁 Template ready for use: {downloaded_path}")
        print("\n💡 Next steps:")
        print("   1. Edit the configuration file to match your simulation needs")
        print("   2. Run: twisstntern-simulate --config-path . --config-name config_template")


@hydra.main(version_base=None, config_path="configs", config_name="config")
def hydra_main(cfg: DictConfig) -> None:
    """
    Main entry point for twisstntern_simulate with Hydra configuration.
    
    Maintains exact functional compatibility with the original argparse version.
    """
    # Setup output directory
    output_dir = Path(cfg.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Setup logging (exactly like original)
    log_file_path = setup_logging(
        output_dir=str(output_dir),
        verbose=cfg.verbose,
        console_output=not cfg.quiet
    )
    
    logger = get_logger(__name__)
    log_system_info()

    # Use Hydra config directly - no conversion needed
    config_source = "hydra_config"
    
    # Convert granularity to float if needed (exactly like original)
    granularity = cfg.analysis.granularity
    try:
        granularity = float(granularity)
    except (ValueError, TypeError):
        pass  # Keep as is if conversion fails
    
    # Parse downsample arguments (using same functions as original)
    downsample_str = cfg.analysis.downsampling.get('tree_downsample', None)
    try:
        downsample_N, downsample_i = None, None
        if cfg.analysis.downsampling.tree_interval is not None:
            downsample_N = cfg.analysis.downsampling.tree_interval
            downsample_i = cfg.analysis.downsampling.tree_start_index
        elif downsample_str:
            downsample_N, downsample_i = parse_downsample_arg(downsample_str)
    except ValueError as e:
        logger.error(f"Invalid downsample argument: {e}")
        sys.exit(1)
    
    # Parse downsampleKB arguments
    downsampleKB_str = cfg.analysis.downsampling.get('position_downsample', None)  
    try:
        downsampleKB_N, downsampleKB_i = None, None
        if cfg.analysis.downsampling.position_interval_kb is not None:
            downsampleKB_N = cfg.analysis.downsampling.position_interval_kb
            downsampleKB_i = cfg.analysis.downsampling.position_start_kb
        elif downsampleKB_str:
            downsampleKB_N, downsampleKB_i = parse_downsampleKB_arg(downsampleKB_str)
    except ValueError as e:
        logger.error(f"Invalid downsampleKB argument: {e}")
        sys.exit(1)
    
    # Log analysis start (exactly like original)
    log_analysis_start(
        input_file=config_source,
        output_dir=str(output_dir),
        granularity=granularity,
        seed_override=cfg.seed,
        mode_override=cfg.simulation.mode if hasattr(cfg, 'simulation') else None,
        topology_mapping=cfg.analysis.topology_mapping,
        verbose=cfg.verbose
    )
    
    start_time = time.time()
    
    # Track files created during this run
    created_files = []
    
    try:
        # Run the pipeline with Hydra config directly
        results = run_pipeline(
            config=cfg,
            output_dir=str(output_dir),
            seed_override=cfg.seed,
            mode_override=cfg.simulation.mode if hasattr(cfg, 'simulation') else None,
            granularity=granularity,
            verbose=cfg.verbose,
            topology_mapping=cfg.analysis.topology_mapping,
            downsample_N=downsample_N,
            downsample_i=downsample_i,
            downsample_kb=downsampleKB_N,
            downsample_kb_i=downsampleKB_i,
            heatmap_colormap=cfg.visualization.heatmap.colormap,
        )
        
        # Calculate duration (exactly like original)
        duration = time.time() - start_time
        
        # Collect files created during this run (exactly like original)
        if output_dir.exists():
            start_timestamp = start_time - 1  # Subtract 1 second to be safe
            
            for file_path in output_dir.glob("*.*"):
                if file_path.stat().st_mtime >= start_timestamp:
                    created_files.append(str(file_path))
        
        # Log completion with only the files created in this run (exactly like original)
        log_analysis_complete(duration, created_files)
        
        # Print summary to console (exactly like original)
        print("----------------------------------------------------------")
        print("Summary of the analysis:")
        print(f"Data file used: {results['csv_file_used']}")
        print(f"\nResults and plots have been saved to the '{output_dir}' directory.")
        
        if log_file_path:
            print(f"Log file saved to: {log_file_path}")
    
    except Exception as e:
        log_error(e, "twisstntern_simulate main analysis")
        logger.critical("Analysis failed. Check the log for details.")
        sys.exit(1)


def main() -> None:
    """
    Entry point that handles special commands before Hydra.
    """
    # Check for special commands first
    if len(sys.argv) > 1 and sys.argv[1] == '--get-config':
        handle_get_config()
        sys.exit(0)
    
    # Otherwise, run the normal Hydra main function
    hydra_main()


if __name__ == "__main__":
    main()
