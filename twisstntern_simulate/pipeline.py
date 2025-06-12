"""
Pipeline orchestrator for twisstntern_simulate package.

This module coordinates the complete workflow:
1. Configuration loading and validation
2. Demographic simulation with msprime
3. Twisst availability checking/downloading  
4. Tree sequence processing to generate topology weights
5. Twisstntern analysis and visualization

The pipeline supports flexible workflow control allowing users to run
specific steps or skip certain components based on their needs.
"""

import os
import logging
from pathlib import Path
from typing import Optional, Dict, Any

from twisstntern_simulate.config import Config
from twisstntern_simulate.simulation import run_simulation
from twisstntern_simulate.download_twisst import ensure_twisst_available
from twisstntern_simulate.ts_processing import ts_to_twisst_weights
# Remove: from twisstntern_simulate.analysis import run_analysis

# Add these imports for direct analysis integration (like main twisstntern)
from twisstntern.analysis import triangles_analysis, fundamental_asymmetry
from twisstntern.visualization import (
    plot,
    plot_results,
    plotting_triangle_index,
    plot_fundamental_asymmetry,
)

# Get logger (logging configured in __init__.py)
logger = logging.getLogger(__name__)


def run_pipeline(
    config_path: str,
    output_dir: str,
    skip_twisst_check: bool = False,
    force_download: bool = False,
    seed_override: Optional[int] = None,
    mode_override: Optional[str] = None,
    granularity: float = 0.1,
    verbose: bool = False,
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

        # Apply command-line overrides
        if seed_override is not None:
            logger.info(f"Overriding seed: {config.seed} -> {seed_override}")
            config.seed = seed_override

        if mode_override is not None:
            logger.info(f"Overriding simulation mode: {config.simulation_mode} -> {mode_override}")
            config.simulation_mode = mode_override

        logger.info(f"Configuration loaded - Mode: {config.simulation_mode}, Seed: {config.seed}")
        results["config"] = config

        # ====================================================================
        # STEP 2: Run simulation
        # ====================================================================
        logger.info("Running msprime simulation...")
            # Ensure Results directory exists
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True)

        simulation_results = run_simulation(config_path, output_dir=output_dir)
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
           # results["twisst_path"] = twisst_path

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
            ts_data, output_file=weights_file, verbose=verbose
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
        
        # Generate output prefix based on simulation mode (like main twisstntern)
        output_prefix = str(output_dir / mode)
        
        logger.info("Generating visualizations...")
        # Generate all plots (exactly like main twisstntern - always generate)
        plot_fundamental_asymmetry(topology_weights, output_prefix)
        plot(topology_weights, granularity, output_prefix)
        plot_results(triangles_results, granularity, output_prefix)
        plotting_triangle_index(granularity, output_prefix)
        
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
        
        # Print summary like main twisstntern
        print("----------------------------------------------------------")
        print("Summary of the analysis:")
        print(f"Data file used: {weights_file}")
        print("\nFundamental asymmetry results:")
        print(f"n_right: {fundamental_results[0]}")
        print(f"n_left: {fundamental_results[1]}")
        print(f"D_LR: {fundamental_results[2]:.4f}")
        print(f"G-test: {fundamental_results[3]:.4f}")
        print(f"p-value: {fundamental_results[4]:.4e}")
        print(f"\nResults and plots have been saved to the '{output_dir}' directory.")

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
