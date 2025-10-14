#!/usr/bin/env python
# coding: utf-8

import shutil
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from omegaconf import DictConfig

from ..core.utils import dump_data
from ..core.analysis import triangles_analysis, fundamental_asymmetry
from . import visualization as viz
from .tree_processing import (
    detect_and_read_trees,
    trees_to_twisst_weights_unified,
    ts_to_twisst_weights,
    newick_to_twisst_weights,
)
from ..core.hydra_utils import build_run_suffix
from ..core.logger import get_logger
from ..core.tree_export import save_ts_CHROM_as_newick, save_newick_strings


def export_tree_archive(
    tree_data,
    tree_type: str,
    destination: Union[Path, str],
    source_path: Optional[str] = None,
    logger=None,
) -> Optional[str]:
    """Persist the tree input for reproducibility and sweep audit."""

    dest_path = Path(destination)
    dest_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        if tree_type == "ts":
            saved_path = save_ts_CHROM_as_newick(tree_data, dest_path)
        else:
            if source_path is not None:
                source = Path(source_path)
                saved_path = str(dest_path)
                try:
                    if source.resolve() != dest_path.resolve():
                        shutil.copy2(source, dest_path)
                except FileNotFoundError:
                    # Fall back to writing from in-memory strings when the source cannot be copied.
                    saved_path = save_newick_strings(tree_data, dest_path)
            else:
                saved_path = save_newick_strings(tree_data, dest_path)
    except Exception as export_error:
        if logger:
            logger.warning(f"Failed to export tree archive: {export_error}")
        else:
            print(f"Warning: could not export tree archive: {export_error}")
        return None

    if logger:
        logger.info(f"Saved tree archive to: {saved_path}")

    return saved_path


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
    logger=None,
    tree_data=None,
    tree_type: Optional[str] = None,
    newick_export_path: Optional[str] = None,
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

    if tree_data is None or tree_type is None:
        tree_data, tree_type = detect_and_read_trees(tree_file)

    exported_archive = None
    if newick_export_path:
        exported_archive = export_tree_archive(
            tree_data,
            tree_type,
            newick_export_path,
            source_path=tree_file,
            logger=logger,
        )

    if tree_type == "ts":
        topology_weights_df = ts_to_twisst_weights(
            tree_data,
            outgroup=cfg.tree_processing.outgroup,
            output_file=str(csv_output),
            verbose=cfg.output.verbose,
            twisst_verbose=cfg.output.verbose,
            topology_mapping=cfg.tree_processing.topology_mapping,
        )
    elif tree_type == "newick":
        topology_weights_df = newick_to_twisst_weights(
            tree_data,
            taxon_names=cfg.tree_processing.taxon_names,
            outgroup=cfg.tree_processing.outgroup,
            output_file=str(csv_output),
            verbose=cfg.output.verbose,
            twisst_verbose=cfg.output.verbose,
            topology_mapping=cfg.tree_processing.topology_mapping,
        )
    else:
        topology_weights_df = trees_to_twisst_weights_unified(
            file_path=tree_file,
            taxon_names=cfg.tree_processing.taxon_names,
            outgroup=cfg.tree_processing.outgroup,
            output_file=str(csv_output),
            verbose=cfg.output.verbose,
            topology_mapping=cfg.tree_processing.topology_mapping,
        )
        if exported_archive is None and newick_export_path:
            export_tree_archive(
                *detect_and_read_trees(tree_file),
                destination=newick_export_path,
                source_path=tree_file,
                logger=logger,
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

        suffix = build_run_suffix(cfg)
        archive_name = f"{Path(cfg.file).stem}.newick" if suffix == "run" else f"{Path(cfg.file).stem}_{suffix}.newick"
        newick_export_path = results_dir / archive_name

        # Process tree file to generate CSV
        if logger:
            logger.info("Converting trees to topology weights...")
        csv_file = process_tree_file(
            tree_file=cfg.file,
            cfg=cfg,
            logger=logger,
            tree_data=tree_data,
            tree_type=tree_type,
            newick_export_path=str(newick_export_path),
        )
        if logger:
            logger.info(f"Generated topology weights CSV: {csv_file}")

        print(f"Tree archive saved to: {newick_export_path}")

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
