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
from omegaconf import DictConfig
from .simulation import run_simulation
from .ts_processing import ts_to_twisst_weights
from ..core.analysis import triangles_analysis, fundamental_asymmetry
from . import visualization as viz
from .visualization import (
    plot_fundamental_asymmetry,
    plot,
    plot_results,
    plotting_triangle_index,
    plot_ternary_heatmap_data,
    plot_density_colored_radcount,
)


# Import logging from twisstntern
from ..core.logger import get_logger, log_analysis_start, log_analysis_complete, log_simulation_config

# Get logger (logging configured in __init__.py)
logger = logging.getLogger(__name__)


def run_pipeline(
    config: DictConfig,
    output_dir: str,
    seed_override: Optional[int] = None,
    mode_override: Optional[str] = None,
    granularity: float = 0.1,
    verbose: bool = False,
    topology_mapping: Optional[str] = None,
    downsample_N: Optional[int] = None,
    downsample_i: Optional[int] = None,
    downsample_kb: Optional[int] = None,
    downsample_kb_i: Optional[int] = None,
    heatmap_colormap: str = "viridis_r",
) -> Dict[str, Any]:
    """
    Run the complete twisstntern_simulate pipeline.

    Args:
        config: Hydra configuration object
        output_dir: Directory for output files
        seed_override: Override random seed from config file
        mode_override: Override simulation mode from config file
        granularity: Granularity for ternary analysis
        verbose: Enable verbose output
        topology_mapping: Custom topology mapping string for T1/T2/T3
        downsample_N: Downsample interval (sample every Nth tree/locus)
        downsample_i: Starting index for downsampling (offset)
        downsample_kb: Downsample interval in kilobases (sample every N kb)
        downsample_kb_i: Starting position in kilobases for KB-based downsampling (offset)
        heatmap_colormap: Colormap for the ternary heatmap. 
                         Options: 'viridis', 'viridis_r', 'plasma', 
                         'inferno', 'Blues', 'Greys'. Default: 'viridis_r'.

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

    # Sync visualization module styling with configuration so plots match legacy output
    if hasattr(config, "visualization"):
        viz.T1_color = config.visualization.colors.T1_color
        viz.T2_color = config.visualization.colors.T2_color
        viz.T3_color = config.visualization.colors.T3_color
        viz.T1_color_data = config.visualization.colors.T1_color_data
        viz.T2_color_data = config.visualization.colors.T2_color_data
        viz.T3_color_data = config.visualization.colors.T3_color_data
        viz.style_heatmap = config.visualization.heatmap.colormap

    results = {
        "config_source": "hydra_config",
        "output_dir": output_dir,
        "errors": [],
    }

    try:
        # ====================================================================
        # STEP 1: Apply configuration overrides
        # ====================================================================
        logger.info("Processing configuration...")
        
        # Apply command-line overrides (Hydra handles most overrides automatically)
        if seed_override is not None:
            logger.info(f"Overriding seed: {config.seed} -> {seed_override}")
            config.seed = seed_override

        if mode_override is not None:
            logger.info(f"Overriding simulation mode: {config.simulation.mode} -> {mode_override}")
            config.simulation.mode = mode_override

        logger.info(f"Configuration processed - Mode: {config.simulation.mode}, Seed: {config.seed}")
        results["config"] = config

        # Log detailed configuration with all overrides
        all_overrides = {}
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
        simulation_results = run_simulation(config, output_dir=output_dir, mode_override=config.simulation.mode)
        results["simulation_results"] = simulation_results

        logger.info(f"Simulation completed")
        
        # Log tree saving
        if "newick_file" in simulation_results:
            logger.info(f"Trees saved to: {simulation_results['newick_file']}")

        # ====================================================================
        # STEP 4: Process tree sequences to generate topology weights
        # ====================================================================
        logger.info("Processing tree sequences...")
        
        # Extract only the TreeSequence objects for processing (exclude newick_file path)
        mode = config.simulation.mode
        if mode == "locus":
            ts_data = simulation_results["locus"]
            n_before = len(ts_data)
            logger.info(f"Locus mode: {n_before} loci simulated.")
            ts_data = list(ts_data)
        elif mode == "chromosome":
            ts_data = simulation_results["chromosome"]
            n_before = ts_data.num_trees
            logger.info(f"Chromosome mode: {n_before} trees in simulated chromosome.")
        else:
            raise ValueError(f"Unknown simulation mode: {mode}")
        
        # Generate output file path for topology weights in user-specified directory
        os.makedirs(output_dir, exist_ok=True)
        weights_file = os.path.join(output_dir, f"{mode}_topology_weights.csv")

        # Create a mapping from numeric population IDs to descriptive labels
        population_labels = None
        if hasattr(config, 'population_labels') and config.population_labels:
            # Create mapping: numeric_id -> descriptive_label 
            # The numeric ID corresponds to the order populations were added to the demography
            population_labels = {}
            populations_with_samples = [
                pop for pop in config.simulation.populations
                if getattr(pop, 'sample_size', 0) is not None
                and getattr(pop, 'sample_size', 0) > 0
            ]
            
            for idx, pop in enumerate(populations_with_samples):
                pop_id = str(idx)  # TreeSequence uses string IDs like "0", "1", "2"
                label = config.population_labels.get(pop.name, pop.name)
                population_labels[pop_id] = label
        
        topology_weights = ts_to_twisst_weights(
            ts_data, output_file=weights_file, verbose=verbose, topology_mapping=topology_mapping, population_labels=population_labels
        )
        results["topology_weights"] = topology_weights

        # === Enhanced downsampling logic (after DataFrame is generated) ===
        n_before_downsample = len(topology_weights)
        logger.info(f"Number of data points before downsampling: {n_before_downsample}")
        topology_weights_downsampled = topology_weights
        downsample_mode = None
        
        if downsample_N is not None and downsample_N > 1:
            if downsample_i is None:
                downsample_i = 0  # Default to starting from index 0
                
            logger.info(f"Downsampling: keeping every {downsample_N}th data point starting from index {downsample_i}.")
            
            # Create the downsampled indices: start from downsample_i, then every downsample_N
            indices = list(range(downsample_i, len(topology_weights), downsample_N))
            topology_weights_downsampled = topology_weights.iloc[indices, :].reset_index(drop=True)
            
            logger.info(f"Downsampled to every {downsample_N}th data point starting from index {downsample_i}. {len(topology_weights_downsampled)} remain.")
            downsample_mode = f"every {downsample_N}th starting from index {downsample_i}"
            
        elif downsample_kb is not None and downsample_kb > 0 and mode == "chromosome":
            if "position" in topology_weights.columns:
                if downsample_kb_i is None:
                    downsample_kb_i = 0  # Default to starting from position 0
                    
                n_kb = downsample_kb * 1000
                start_pos = downsample_kb_i  # Already in base pairs from the new parser
                max_pos = topology_weights["position"].max()
                
                # Generate target positions starting from the offset
                target_positions = list(range(start_pos, int(max_pos) + 1, n_kb))
                selected_indices = []
                last_idx = -1
                
                for pos in target_positions:
                    candidates = topology_weights.index[topology_weights["position"] >= pos]
                    if len(candidates) > 0:
                        idx = candidates[0]
                        if idx != last_idx:
                            selected_indices.append(idx)
                            last_idx = idx
                            
                topology_weights_downsampled = topology_weights.loc[selected_indices].reset_index(drop=True)
                
                # Format the starting position for display
                if downsample_kb_i >= 1000000:
                    start_display = f"{downsample_kb_i/1000000:.1f}mb"
                elif downsample_kb_i >= 1000:
                    start_display = f"{downsample_kb_i/1000:.0f}kb"
                else:
                    start_display = f"{downsample_kb_i}bp"
                    
                logger.info(f"Downsampled to one data point every {downsample_kb} kb starting from {start_display}. {len(topology_weights_downsampled)} remain.")
                downsample_mode = f"every {downsample_kb} kb starting from {start_display}"
            else:
                logger.warning("KB downsampling requested but no 'position' column found in topology weights. Skipping KB downsampling.")
                
        n_after_downsample = len(topology_weights_downsampled)
        if downsample_mode:
            logger.info(f"Number of data points after downsampling ({downsample_mode}): {n_after_downsample}")
        else:
            logger.info(f"No downsampling applied. Data points used: {n_after_downsample}")

        # Use downsampled data for all further analysis
        topology_weights = topology_weights_downsampled

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
        
        # Add tree count information for both modes
        n_total_trees = len(topology_weights)
        n_right = fundamental_results[0] 
        n_left = fundamental_results[1]
        n_trees_used = n_right + n_left
        
        if mode == "chromosome":
            logger.info(f"Total trees computed: {n_total_trees} (n_right + n_left = {n_trees_used})")
            
            # Always explain T2=T3 filtering (even when 0)
            n_filtered = n_total_trees - n_trees_used
            logger.info(f"Note: {n_filtered} trees were filtered out (trees where T2 = T3, which fall exactly on the y-axis in ternary space)")
                
        elif mode == "locus":
            logger.info(f"n_loci simulated: {config.simulation.n_loci}, total trees analyzed: {n_total_trees} (n_right + n_left = {n_trees_used})")
            
            # Always explain T2=T3 filtering (even when 0)
            n_filtered = n_total_trees - n_trees_used
            logger.info(f"Note: {n_filtered} trees were filtered out (trees where T2 = T3, which fall exactly on the y-axis in ternary space)")
            
            # Also show if loci were filtered during processing
            if n_total_trees < config.simulation.n_loci:
                n_loci_filtered = config.simulation.n_loci - n_total_trees
                logger.info(f"Additional note: {n_loci_filtered} loci were filtered out during processing.")
        
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
        # Ternary heatmap (count, no grid) - always uses 0.02 granularity
        plot_ternary_heatmap_data(topology_weights, 0.02, output_prefix, heatmap_colormap=heatmap_colormap)
        # Density radcount plot - mirrors legacy colormap behaviour
        density_colormap = "viridis"
        if hasattr(config, "visualization") and hasattr(config.visualization, "density"):
            density_colormap = config.visualization.density.colormap
        plot_density_colored_radcount(topology_weights, output_prefix, colormap=density_colormap)
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
        # Convert granularity to float for filename
        if granularity == "superfine":
            alpha = 0.05
        elif granularity == "fine":
            alpha = 0.1
        elif granularity == "coarse":
            alpha = 0.25
        else:
            alpha = float(granularity)
        
        triangles_csv = output_dir / f"{mode}_triangle_analysis_{alpha}.csv"
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
        logger.exception(error_msg)
        results["errors"].append(str(e))
        raise RuntimeError(error_msg) from e


def _generate_summary_report(results: Dict[str, Any], output_file: str) -> None:
    """Generate a summary report of the pipeline run."""

    with open(output_file, "w") as f:
        f.write("TWISSTNTERN_SIMULATE PIPELINE SUMMARY\n")
        f.write("=" * 50 + "\n\n")

        f.write(f"Configuration source: {results['config_source']}\n")
        f.write(f"Output directory: {results['output_dir']}\n\n")

        if "config" in results:
            config = results["config"]
            f.write("CONFIGURATION:\n")
            f.write(f"  Simulation mode: {config.simulation.mode}\n")
            f.write(f"  Random seed: {config.seed}\n")
            f.write(f"  Number of populations: {len(config.simulation.populations)}\n")
            if hasattr(config.simulation, "n_loci"):
                f.write(f"  Number of loci: {config.simulation.n_loci}\n")
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
