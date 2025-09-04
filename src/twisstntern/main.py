#!/usr/bin/env python
# coding: utf-8

import hydra
from omegaconf import DictConfig, OmegaConf
from pathlib import Path
import sys
import time

from .config import TwisstnternConfig
from .pipeline import run_analysis, ensure_twisst_available
from .logger import setup_logging, get_logger, log_system_info, log_analysis_start, log_analysis_complete, log_error


@hydra.main(version_base=None, config_path="configs", config_name="default")
def main(cfg: DictConfig) -> None:
    """
    Main entry point for TWISSTNTERN analysis pipeline using Hydra configuration.
    
    Args:
        cfg: Hydra configuration object
    """
    # Ensure twisst is available before any analysis
    ensure_twisst_available()
    
    # Validate required parameters
    if not cfg.file:
        print("Error: Input file must be specified")
        print("Either set 'file' in config or use: python -m twisstntern file=path/to/file")
        sys.exit(1)
    
    # Check if input file exists
    file_path = Path(cfg.file)
    if not file_path.exists():
        print(f"Error: Input file not found: {file_path}")
        sys.exit(1)
    
    # Set up logging if requested
    logger = None
    log_file_path = None
    if cfg.output.log_file:
        log_file_path = setup_logging(
            output_dir=cfg.output.output_dir,
            verbose=cfg.output.verbose,
            console_output=True
        )
        logger = get_logger(__name__)
        log_system_info()
    
    # Log analysis start
    if logger:
        log_analysis_start(
            input_file=str(file_path),
            output_dir=cfg.output.output_dir,
            granularity=cfg.analysis.granularity,
            taxon_names=cfg.tree_processing.taxon_names,
            outgroup=cfg.tree_processing.outgroup,
            topology_mapping=cfg.tree_processing.topology_mapping,
            verbose=cfg.output.verbose
        )
    
    start_time = time.time()
    
    # Track files created during this run
    created_files = []
    
    try:
        # Run the analysis pipeline
        results, fundamental_results, csv_file_used = run_analysis(cfg)
        
        # Calculate duration
        duration = time.time() - start_time
        
        # Collect files created during this run
        output_path = Path(cfg.output.output_dir)
        if output_path.exists():
            start_timestamp = start_time - 1
            for file_path in output_path.glob("*.*"):
                if file_path.stat().st_mtime >= start_timestamp:
                    created_files.append(str(file_path))
        
        # Log completion
        if logger:
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
        print(f"\nResults and plots have been saved to the '{cfg.output.output_dir}' directory.")
        
        if log_file_path:
            print(f"Log file saved to: {log_file_path}")
            
    except Exception as e:
        if logger:
            log_error(e, "main analysis")
            logger.critical("Analysis failed. Check the log for details.")
        else:
            print(f"Analysis failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()