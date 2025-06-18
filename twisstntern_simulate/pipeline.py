"""
TWISSTNTERN Simulate Pipeline

This module provides the main pipeline functionality for running the complete
twisstntern_simulate workflow: simulation -> processing -> analysis
"""

import os
import logging
import pandas as pd
from pathlib import Path
from typing import Dict, Any, Optional

# Import simulation components
from twisstntern_simulate.config import Config
from twisstntern_simulate.simulation import run_simulation
from twisstntern_simulate.download_twisst import ensure_twisst_available
from twisstntern_simulate.ts_processing import ts_to_twisst_weights
from twisstntern_simulate.analysis import triangles_analysis, fundamental_asymmetry
from twisstntern_simulate.visualization import (
    plot_fundamental_asymmetry,
    plot,
    plot_results,
    plotting_triangle_index,
)

# Import logging from twisstntern
from twisstntern.logger import get_logger, log_analysis_start, log_analysis_complete, log_simulation_config

# Get logger (logging configured in __init__.py)
logger = logging.getLogger(__name__)


def apply_config_overrides(config, overrides_list):
    """
    Apply configuration overrides to a config object.
    
    Args:
        config: Configuration object to modify
        overrides_list: List of override strings in format 'key=value' or 'nested.key=value'
        
    Returns:
        dict: Dictionary of applied overrides for logging
    """
    if not overrides_list:
        return {}
    
    applied_overrides = {}
    
    for override_str in overrides_list:
        if '=' not in override_str:
            logger.warning(f"Invalid override format: {override_str}. Expected 'key=value'")
            continue
            
        key_path, value = override_str.split('=', 1)
        
        try:
            # Convert value to appropriate type
            if value.lower() in ['true', 'false']:
                value = value.lower() == 'true'
            elif value.replace('.', '').replace('-', '').isdigit():
                value = float(value) if '.' in value else int(value)
            
            # Handle nested keys like 'migration.p1>p2' or 'populations.p1.Ne'
            if '.' in key_path:
                parts = key_path.split('.')
                
                if parts[0] == 'migration':
                    # Handle migration overrides: migration.p1>p2=0.5
                    migration_key = parts[1]
                    if hasattr(config, 'migration') and config.migration:
                        old_value = config.migration.get(migration_key, 0.0)
                        config.migration[migration_key] = value
                        applied_overrides[f"migration.{migration_key}"] = f"{old_value} -> {value}"
                        logger.info(f"Override applied: migration.{migration_key}: {old_value} -> {value}")
                    else:
                        logger.warning(f"Migration configuration not found for override: {override_str}")
                        
                elif parts[0] == 'populations':
                    # Handle population overrides: populations.p1.Ne=15000
                    if len(parts) >= 3:
                        pop_name = parts[1]
                        pop_attr = parts[2]
                        
                        # Find the population in the list
                        for pop in config.populations:
                            if pop.name == pop_name:
                                old_value = getattr(pop, pop_attr, None)
                                setattr(pop, pop_attr, value)
                                applied_overrides[f"populations.{pop_name}.{pop_attr}"] = f"{old_value} -> {value}"
                                logger.info(f"Override applied: populations.{pop_name}.{pop_attr}: {old_value} -> {value}")
                                break
                        else:
                            logger.warning(f"Population '{pop_name}' not found for override: {override_str}")
                    else:
                        logger.warning(f"Invalid population override format: {override_str}")
                else:
                    logger.warning(f"Unsupported nested override: {override_str}")
            else:
                # Handle top-level overrides like 'ploidy=2', 'seed=1234'
                if hasattr(config, key_path):
                    old_value = getattr(config, key_path)
                    setattr(config, key_path, value)
                    applied_overrides[key_path] = f"{old_value} -> {value}"
                    logger.info(f"Override applied: {key_path}: {old_value} -> {value}")
                else:
                    logger.warning(f"Unknown configuration key: {key_path}")
                    
        except Exception as e:
            logger.error(f"Failed to apply override '{override_str}': {e}")
    
    return applied_overrides


def run_pipeline(
    config_path: str,
    output_dir: str,
    skip_twisst_check: bool = False,
    force_download: bool = False,
    seed_override: Optional[int] = None,
    mode_override: Optional[str] = None,
    granularity: float = 0.1,
    verbose: bool = False,
    topology_mapping: Optional[str] = None,
    config_overrides: Optional[list] = None,
) -> Dict[str, Any]:
    """
    Run the complete twisstntern_simulate pipeline.

    Args:
        config_path: Path to configuration YAML file
        output_dir: Directory for output files
        skip_twisst_check: If True, skip checking/downloading twisst
        force_download: If True, force re-download of twisst
        seed_override: Override random seed from config file
        mode_override: Override simulation mode from config file
        granularity: Granularity for ternary analysis
        verbose: Enable verbose output
        topology_mapping: Custom topology mapping string for T1/T2/T3
        config_overrides: List of config override strings in format 'key=value'

    Returns:
        Dictionary containing pipeline results and metadata

    Raises:
        FileNotFoundError: If config file or required inputs don't exist
        ValueError: If configuration is invalid
        RuntimeError: If any pipeline step fails
    """

    logger.info("---------------------------------------------------------")
    logger.info("Starting twisstntern_simulate pipeline")
    logger.info("---------------------------------------------------------")

    results = {
        "config_path": config_path,
        "output_dir": output_dir,
        "errors": [],
    }

    try:
        # ====================================================================
        # STEP 1: Load and validate configuration
        # ====================================================================
        logger.info("Loading configuration...")

        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Configuration file not found: {config_path}") 
        config = Config(config_path)
        
        # Apply configuration overrides first
        applied_config_overrides = apply_config_overrides(config, config_overrides)
        
        # Apply command-line overrides
        if seed_override is not None:
            logger.info(f"Overriding seed: {config.seed} -> {seed_override}")
            config.seed = seed_override

        if mode_override is not None:
            logger.info(f"Overriding simulation mode: {config.simulation_mode} -> {mode_override}")
            config.simulation_mode = mode_override

        logger.info(f"Configuration loaded - Mode: {config.simulation_mode}, Seed: {config.seed}")
        results["config"] = config

        # Log detailed configuration with all overrides
        all_overrides = applied_config_overrides.copy()
        if seed_override is not None:
            all_overrides["seed_override"] = seed_override
        if mode_override is not None:
            all_overrides["mode_override"] = mode_override
        if topology_mapping is not None:
            all_overrides["topology_mapping"] = topology_mapping
        log_simulation_config(config, all_overrides)

        # ====================================================================
        # STEP 2: Run simulation
        # ====================================================================
        logger.info("Running msprime simulation...")
        # Ensure Results directory exists
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)
        simulation_results = run_simulation(config_path, output_dir=output_dir, mode_override=config.simulation_mode)
        results["simulation_results"] = simulation_results

        logger.info(f"Simulation completed")
        
        # Log tree saving
        if "newick_file" in simulation_results:
            logger.info(f"Trees saved to: {simulation_results['newick_file']}")

        # ====================================================================
        # STEP 3: Ensure twisst is available (unless skipped)
        # ====================================================================
        if not skip_twisst_check:
            logger.info("Checking twisst availability...")
            twisst_path = ensure_twisst_available()
        else:
            logger.info("Skipped twisst check")

        # ====================================================================
        # STEP 4: Process tree sequences to generate topology weights
        # ====================================================================
        logger.info("Processing tree sequences...")
        
        # Extract only the TreeSequence objects for processing (exclude newick_file path)
        mode = config.simulation_mode
        if mode == "locus":
            ts_data = simulation_results["locus"]
        elif mode == "chromosome":
            ts_data = simulation_results["chromosome"]
        else:
            raise ValueError(f"Unknown simulation mode: {mode}")
        
        # Generate output file path for topology weights in user-specified directory
        os.makedirs(output_dir, exist_ok=True)
        weights_file = os.path.join(output_dir, f"{mode}_topology_weights.csv")

        topology_weights = ts_to_twisst_weights(
            ts_data, output_file=weights_file, verbose=verbose, topology_mapping=topology_mapping
        )
        results["topology_weights"] = topology_weights

        logger.info(f"Tree processing completed - weights saved to: {weights_file}")

        # ====================================================================
        # STEP 5: Run twisstntern analysis and visualization (integrated like main twisstntern)
        # ====================================================================
        logger.info("Running triangle analysis...")
        triangles_results = triangles_analysis(topology_weights, granularity)
        logger.info("Running fundamental asymmetry analysis...")
        fundamental_results = fundamental_asymmetry(topology_weights)
        
        # Log fundamental asymmetry results
        logger.info("="*60)
        logger.info("FUNDAMENTAL ASYMMETRY RESULTS")
        logger.info("="*60)
        logger.info(f"Data file used: {weights_file}")
        logger.info(f"n_right: {fundamental_results[0]}")
        logger.info(f"n_left: {fundamental_results[1]}")
        logger.info(f"D_LR: {fundamental_results[2]:.4f}")
        logger.info(f"G-test: {fundamental_results[3]:.4f}")
        logger.info(f"p-value: {fundamental_results[4]:.4e}")
        logger.info("="*60)

        # Generate output prefix based on simulation mode (like main twisstntern)
        output_prefix = str(output_dir / mode)
        
        logger.info("Generating visualizations...")
        # Generate all plots (exactly like main twisstntern - always generate)
        plot_fundamental_asymmetry(topology_weights, output_prefix)
        plot(topology_weights, granularity, output_prefix)
        plot_results(triangles_results, granularity, output_prefix)
        plotting_triangle_index(granularity, output_prefix)

        # ====================================================================
        # Add main sub-triangle results to the dataframe
        new_row = pd.DataFrame(
            [["main subtriangle", *fundamental_results[:5], 'NA']],
            columns=triangles_results.columns,
            index=['full dataset']
        )
        triangles_results = pd.concat([new_row, triangles_results])

        # Save results as CSV (exactly like main twisstntern naming)
        triangles_csv = output_dir / f"{mode}_triangle_analysis.csv"
        triangles_results.to_csv(triangles_csv, index=False, float_format="%.3f")
        logger.info(f"Saved triangle analysis results to: {triangles_csv}")
        
        # Store results (same format as main twisstntern)
        results["triangles_results"] = triangles_results
        results["fundamental_results"] = fundamental_results
        results["csv_file_used"] = weights_file
        
        logger.info("Analysis pipeline completed successfully!")
        logger.info(f"Pipeline completed successfully - Results in: {output_dir}")

        return results

    except Exception as e:
        error_msg = f"Pipeline failed: {str(e)}"
        logger.error(error_msg)
        results["errors"].append(str(e))
        raise RuntimeError(error_msg) from e


def _generate_summary_report(results: Dict[str, Any], output_file: str) -> None:
    """Generate a summary report of the pipeline run."""

    with open(output_file, "w") as f:
        f.write("TWISSTNTERN_SIMULATE PIPELINE SUMMARY\n")
        f.write("=" * 50 + "\n\n")

        f.write(f"Configuration file: {results['config_path']}\n")
        f.write(f"Output directory: {results['output_dir']}\n\n")

        if "config" in results:
            config = results["config"]
            f.write("CONFIGURATION:\n")
            f.write(f"  Simulation mode: {config.simulation_mode}\n")
            f.write(f"  Random seed: {config.seed}\n")
            f.write(f"  Number of populations: {len(config.populations)}\n")
            if hasattr(config, "n_loci"):
                f.write(f"  Number of loci: {config.n_loci}\n")
            f.write("\n")

        if "topology_weights" in results:
            f.write("TOPOLOGY WEIGHTS:\n")
            f.write(f"  Generated topology weights successfully\n\n")

        if "analysis_results" in results:
            f.write("ANALYSIS RESULTS:\n")
            f.write(f"  Analysis completed successfully\n\n")

        if results["errors"]:
            f.write("ERRORS:\n")
            for error in results["errors"]:
                f.write(f"  - {error}\n")
        else:
            f.write("No errors encountered.\n")

    logger.info(f"Summary report written to: {output_file}") 