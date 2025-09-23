#!/usr/bin/env python
# coding: utf-8

from pathlib import Path
import pandas as pd
from omegaconf import DictConfig

from ..core.utils import dump_data
from ..core.analysis import triangles_analysis, fundamental_asymmetry
from . import visualization as viz
from .tree_processing import (
    detect_and_read_trees,
    trees_to_twisst_weights_unified,
)
from ..core.logger import get_logger


def detect_file_type(file_path):
    """
    Detect whether the input file is a tree file or CSV file based on extension.

    Args:
        file_path (str): Path to the input file

    Returns:
        str: "tree" for tree files, "csv" for CSV files
    """
    file_path = Path(file_path)
    tree_extensions = {".trees", ".ts", ".newick", ".nwk", ".tree", ".nexus"}
    csv_extensions = {".csv"}

    if file_path.suffix.lower() in tree_extensions:
        return "tree"
    elif file_path.suffix.lower() in csv_extensions:
        return "csv"
    else:
        raise ValueError(
            f"Unsupported file format: {file_path.suffix}. "
            f"Supported formats: {tree_extensions | csv_extensions}"
        )


def ensure_twisst_available():
    """Ensure twisst is available for use."""
    try:
        from .tree_processing import weightTrees
        return True
    except ImportError:
        print("✗ twisst not available in the package")
        return False


def process_tree_file(
    tree_file: str,
    cfg: DictConfig,
    logger=None
):
    """
    Process a tree file to generate topology weights CSV file.

    Args:
        tree_file (str): Path to the tree file
        cfg (DictConfig): Hydra configuration
        logger: Optional logger instance

    Returns:
        str: Path to the generated CSV file with topology weights
    """
    # Ensure output directory exists
    output_dir = Path(cfg.output.output_dir)
    output_dir.mkdir(exist_ok=True)

    # Ensure twisst is available before processing
    if not ensure_twisst_available():
        raise RuntimeError(
            "twisst is required for tree file processing but could not be made available. "
            "Please install it manually or check your internet connection."
        )

    # Generate output filename
    input_name = Path(tree_file).stem
    csv_output = output_dir / f"{input_name}_topology_weights.csv"

    if logger:
        logger.info(f"Processing tree file: {tree_file}")
        logger.info(f"Output CSV will be saved to: {csv_output}")

    # Process trees using the unified function
    topology_weights_df = trees_to_twisst_weights_unified(
        file_path=tree_file,
        taxon_names=cfg.tree_processing.taxon_names,
        outgroup=cfg.tree_processing.outgroup,
        output_file=str(csv_output),
        verbose=cfg.output.verbose,
        topology_mapping=cfg.tree_processing.topology_mapping,
    )

    if logger:
        logger.info(f"✓ Successfully generated topology weights CSV: {csv_output}")
        logger.info(f"  - Shape: {topology_weights_df.shape}")
        logger.info(f"  - Columns: {list(topology_weights_df.columns)}")

    return str(csv_output)


def run_analysis(cfg: DictConfig):
    """
    Orchestrates the full analysis and visualization pipeline for both tree files and CSV files.

    Args:
        cfg (DictConfig): Hydra configuration containing all parameters

    Returns:
        tuple: (results, fundamental_results, csv_file_used)
    """
    logger = get_logger(__name__) if cfg.output.log_file else None

    # Sync visualization module styling with configuration to mirror legacy defaults
    viz.style = cfg.visualization.style
    viz.style_heatmap = cfg.visualization.style_heatmap
    viz.T1_color = cfg.visualization.t1_color
    viz.T2_color = cfg.visualization.t2_color
    viz.T3_color = cfg.visualization.t3_color
    viz.T1_color_data = cfg.visualization.t1_color_data
    viz.T2_color_data = cfg.visualization.t2_color_data
    viz.T3_color_data = cfg.visualization.t3_color_data
    
    # Ensure Results directory exists
    results_dir = Path(cfg.output.output_dir)
    results_dir.mkdir(exist_ok=True)
    if logger:
        logger.info(f"Output directory: {results_dir}")

    # Detect file type and process accordingly
    file_type = detect_file_type(cfg.file)
    if logger:
        logger.info(f"Detected file type: {file_type}")

    if file_type == "tree":
        if logger:
            logger.info("Processing tree file to generate topology weights...")
        print("Detected tree file format. Processing trees to generate topology weights...")

        # Validate tree file parameters
        tree_data, tree_type = detect_and_read_trees(cfg.file)
        if logger:
            logger.info(f"Detected tree format: {tree_type}")

        # Validate Newick file requirements
        if tree_type == "newick":
            if logger:
                logger.debug("Validating Newick file parameters...")
            if cfg.tree_processing.taxon_names is None:
                raise ValueError(
                    "❌ Taxon names are required for Newick files!\n"
                    "   Configure in config file or command line"
                )
            if cfg.tree_processing.outgroup is None:
                raise ValueError(
                    "❌ Outgroup is required for Newick files!\n"
                    "   Configure in config file or command line"
                )
            if logger:
                logger.info(f"Using taxon names: {cfg.tree_processing.taxon_names}")
                logger.info(f"Using outgroup: {cfg.tree_processing.outgroup}")

        # Process tree file to generate CSV
        if logger:
            logger.info("Converting trees to topology weights...")
        csv_file = process_tree_file(
            tree_file=cfg.file,
            cfg=cfg,
            logger=logger
        )
        if logger:
            logger.info(f"Generated topology weights CSV: {csv_file}")

        print("Tree processing complete.")

    elif file_type == "csv":
        if logger:
            logger.info("Using CSV file directly for analysis")
        print("Detected CSV file format. Using file directly for analysis...")
        csv_file = cfg.file

    else:
        error_msg = f"Unsupported file type: {file_type}"
        if logger:
            logger.error(error_msg)
        raise ValueError(error_msg)

    # Load and process the CSV data
    if logger:
        logger.info(f"Loading data from: {csv_file}")
    print(f"Loading data from: {csv_file}")
    axis_order = cfg.processing.axis_order if cfg.processing.axis_order else ["T1", "T2", "T3"]
    data = dump_data(
        csv_file,
        logger=logger,
        axis_order=axis_order,
        normalize=cfg.processing.normalize_data,
        remove_equal_t2_t3=cfg.processing.remove_equal_t2_t3,
    )
    n_before_trim = len(data)
    if logger:
        logger.info(f"Loaded data shape: {data.shape}")
        logger.debug(f"Data columns: {list(data.columns)}")

    # Apply downsampling if configured
    if cfg.processing.downsample_n is not None and cfg.processing.downsample_n > 1:
        downsample_i = cfg.processing.downsample_i if cfg.processing.downsample_i is not None else 0
        
        if logger:
            logger.info(f"Downsampling: keeping every {cfg.processing.downsample_n}th row starting from index {downsample_i}.")
        print(f"Downsampling: keeping every {cfg.processing.downsample_n}th row starting from index {downsample_i}.")
        
        # Create downsampled indices
        indices = list(range(downsample_i, len(data), cfg.processing.downsample_n))
        data_trimmed = data.iloc[indices, :].reset_index(drop=True)
        
        trimmed_csv_file = str(Path(csv_file).with_name(Path(csv_file).stem + "_trimmed.csv"))
        data_trimmed.to_csv(trimmed_csv_file, index=False)
        if logger:
            logger.info(f"Trimmed topology weights saved to: {trimmed_csv_file}")
        n_after_trim = len(data_trimmed)
        if logger:
            logger.info(f"Number of data points after downsampling: {n_after_trim}")
        data = data_trimmed
        csv_file = trimmed_csv_file
    else:
        n_after_trim = n_before_trim

    # Run triangle analysis
    if logger:
        logger.info("Running triangle analysis...")
    print("Running triangle analysis...")
    results = triangles_analysis(data, cfg.analysis.granularity)
    if logger:
        logger.info(f"Triangle analysis completed. Results shape: {results.shape}")

    # Run fundamental asymmetry analysis
    if logger:
        logger.info("Running fundamental asymmetry analysis...")
    print("Running fundamental asymmetry analysis...")
    fundamental_results = fundamental_asymmetry(data)
    n_right = fundamental_results[0]
    n_left = fundamental_results[1]
    n_used = n_right + n_left
    n_filtered = n_after_trim - n_used

    # Log fundamental asymmetry results
    if logger:
        logger.info("="*60)
        logger.info("FUNDAMENTAL ASYMMETRY RESULTS")
        logger.info("="*60)
        logger.info(f"Data file used: {csv_file}")
        logger.info(f"Total data points before downsampling: {n_before_trim}")
        if cfg.processing.downsample_n is not None and cfg.processing.downsample_n > 1:
            logger.info(f"Total data points after downsampling: {n_after_trim}")
        logger.info(f"Total data points used in symmetry analysis: {n_used} (n_right + n_left = {n_used})")
        logger.info(f"Note: {n_filtered} data points were filtered out (where T2 = T3)")
        logger.info(f"n_right: {n_right}")
        logger.info(f"n_left: {n_left}")
        logger.info(f"D_LR: {fundamental_results[2]:.4f}")
        logger.info(f"G-test: {fundamental_results[3]:.4f}")
        logger.info(f"p-value: {fundamental_results[4]:.4e}")
        logger.info("="*60)

    # Generate output prefix
    output_prefix = str(results_dir / Path(cfg.file).stem)
    if logger:
        logger.debug(f"Output prefix: {output_prefix}")

    # Generate all visualizations
    if logger:
        logger.info("Generating visualizations...")
    print("Generating visualizations...")
    
    # Update visualization functions to accept config parameters
    viz.plot_fundamental_asymmetry(data, output_prefix)
    if logger:
        logger.debug("Generated fundamental asymmetry plot")
        
    viz.plot(data, cfg.analysis.granularity, output_prefix)
    if logger:
        logger.debug("Generated ternary plot")

    # Ternary heatmap with fixed granularity and configurable colormap
    viz.plot_ternary_heatmap_data(
        data,
        cfg.analysis.heatmap_granularity,
        output_prefix,
        heatmap_colormap=cfg.visualization.style_heatmap,
    )
    if logger:
        logger.debug("Generated ternary heatmap")

    # Density radcount plot
    viz.plot_density_colored_radcount(
        data,
        output_prefix,
        colormap=cfg.visualization.style_heatmap,
    )
    if logger:
        logger.debug("Generated density radcount plot")
    
    viz.plot_results(results, cfg.analysis.granularity, output_prefix)
    if logger:
        logger.debug("Generated results plot")
        
    viz.plotting_triangle_index(cfg.analysis.granularity, output_prefix)
    if logger:
        logger.debug("Generated triangle index plot")

    # Add main sub-triangle results to the dataframe
    new_row = pd.DataFrame(
        [["main subtriangle", fundamental_results[0], fundamental_results[1], fundamental_results[2], fundamental_results[3], fundamental_results[4], "NA"]],
        columns=["coord. (T1, T2, T3)", "n_right", "n_left", "D-LR", "g-test", "p-value(g-test)", "index"],
        index=["full dataset"],
    )
    results = pd.concat([new_row, results])

    # Save results
    if logger:
        logger.info("Saving results...")
    
    # Convert granularity to float for filename
    granularity_names = {
        "superfine": cfg.analysis.superfine_granularity,
        "fine": cfg.analysis.fine_granularity,
        "coarse": cfg.analysis.coarse_granularity
    }
    
    if isinstance(cfg.analysis.granularity, str) and cfg.analysis.granularity in granularity_names:
        alpha = granularity_names[cfg.analysis.granularity]
    else:
        alpha = float(cfg.analysis.granularity)
    
    results_csv = results_dir / f"{Path(cfg.file).stem}_triangle_analysis_{alpha}.csv"
    results.to_csv(results_csv, index=False, float_format="%.3f")
    if logger:
        logger.info(f"Saved triangle analysis results to: {results_csv}")
    print(f"Saved triangle analysis results to: {results_csv}")

    if logger:
        logger.info("Analysis pipeline completed successfully!")
    print("Analysis pipeline completed successfully!")

    return results, fundamental_results, csv_file
