#!/usr/bin/env python
# coding: utf-8

from pathlib import Path
import pandas as pd

from twisstntern.utils import dump_data
from twisstntern.analysis import triangles_analysis, fundamental_asymmetry
from twisstntern.visualization import (
    plot,
    plot_results,
    plotting_triangle_index,
    plot_fundamental_asymmetry,
)
from twisstntern.tree_processing import (
    detect_and_read_trees,
    trees_to_twisst_weights_unified,
)
from twisstntern.logger import get_logger


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
    """
    Ensure twisst is available for tree processing. Download if necessary.

    Returns:
        bool: True if twisst is available, False otherwise
    """
    try:
        # Try importing twisst to check if it's available
        import sys
        from pathlib import Path

        # Add external directory to path
        external_dir = Path(__file__).parent / "external"
        sys.path.append(str(external_dir))

        from twisst import weightTrees  # twisst might be in users direct directory

        print("✓ twisst is already available")
        return True

    except ImportError:
        print("⚠️  twisst not found. Downloading automatically...")

        try:
            # Import and run the download function
            from twisstntern.download_twisst import (
                ensure_twisst_available as download_twisst,
            )

            success = download_twisst()  # download twisst to external directory

            if success:
                print("✓ twisst downloaded successfully")
                return True
            else:
                print("✗ Failed to download twisst automatically")
                print("Please run: python -m twisstntern.download_twisst")
                return False

        except Exception as e:
            print(f"✗ Error downloading twisst: {e}")
            print("Please run: python -m twisstntern.download_twisst")
            return False


def process_tree_file(
    tree_file,
    taxon_names=None,
    outgroup=None,
    output_dir="Results",
    verbose=True,
    topology_mapping=None,
):
    """
    Process a tree file to generate topology weights CSV file.

    Args:
        tree_file (str): Path to the tree file
        taxon_names (list, optional): List of taxon names for Newick files
        outgroup (str, optional): Outgroup taxon name
        output_dir (str): Directory to save the output CSV file
        verbose (bool): Whether to print verbose output
        topology_mapping (str, optional): User-defined topology mapping string

    Returns:
        str: Path to the generated CSV file with topology weights
    """
    # Ensure output directory exists
    output_dir = Path(output_dir)
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

    # print(f"Processing tree file: {tree_file}")
    print(f"Output CSV will be saved to: {csv_output}")

    # Process trees using the unified function, which automatically saves csv topology weights
    topology_weights_df = trees_to_twisst_weights_unified(
        file_path=tree_file,
        taxon_names=taxon_names,
        outgroup=outgroup,
        output_file=str(csv_output),
        verbose=verbose,
        topology_mapping=topology_mapping,
    )

    print(f"✓ Successfully generated topology weights CSV: {csv_output}")
    print(f"  - Shape: {topology_weights_df.shape}")
    print(f"  - Columns: {list(topology_weights_df.columns)}")

    return str(csv_output)


# default granularity is 0.1
def run_analysis(
    file,
    granularity=0.1,
    taxon_names=None,
    outgroup=None,
    output_dir="Results",
    topology_mapping=None,
):
    """
    Orchestrates the full analysis and visualization pipeline for both tree files and CSV files.

    Args:
        file (str): Path to the input file (tree file or CSV file).
        granularity (str or float): Granularity level ("superfine", "fine", "coarse", or a float).
        taxon_names (list, optional): List of taxon names for Newick tree files.
                                     Ignored for TreeSequence files and CSV files.
        outgroup (str, optional): Outgroup taxon name for tree files.
                                 Ignored for CSV files.
        output_dir (str): Path to the output directory where the results should be saved
        topology_mapping (str, optional): User-defined topology mapping for custom topology ordering.
                                         Format: 'T1="(0,(3,(1,2)))"; T2="(0,(1,(2,3)))"; T3="(0,(2,(1,3)))";'
                                         Ignored for CSV files.

    Returns:
        tuple: (results, fundamental_results, csv_file_used)
            - results: Triangle analysis results
            - fundamental_results: Fundamental asymmetry results
            - csv_file_used: Path to the CSV file that was analyzed (original or generated)
    """
    logger = get_logger(__name__)
    
    # Ensure Results directory exists
    results_dir = Path(output_dir)
    results_dir.mkdir(exist_ok=True)
    logger.info(f"Output directory: {results_dir}")

    # Detect file type and process accordingly
    file_type = detect_file_type(file)
    logger.info(f"Detected file type: {file_type}")

    if file_type == "tree":
        logger.info("Processing tree file to generate topology weights...")
        print(
            "Detected tree file format. Processing trees to generate topology weights..."
        )

        # if the tree is in Newick format, the users had to provide taxon names
        tree_data, tree_type = detect_and_read_trees(file)
        logger.info(f"Detected tree format: {tree_type}")

        # Add ploidy detection here
        if tree_type == "newick":
            logger.debug(f"Validating Newick file parameters...")
            if taxon_names is None:
                logger.error("Taxon names are required for Newick files")
                raise ValueError(
                    "❌ Taxon names are required for Newick files!\n"
                    "   Use --taxon-names to specify population names, e.g.:\n"
                    "   --taxon-names O,P1,P2,P3\n"
                    "   \n"
                    "   For Newick files, you must explicitly specify which population names\n"
                    "   correspond to your tree samples (e.g., O_1, P1_5 → populations O, P1)"
                )
            if outgroup is None:
                logger.error("Outgroup is required for Newick files")
                raise ValueError(
                    "❌ Outgroup is required for Newick files!\n"
                    "   Use --outgroup to specify the outgroup population, e.g.:\n"
                    "   --outgroup O\n"
                    "   \n"
                    "   The outgroup should be one of your taxon names and represents\n"
                    "   the ancestral/reference population for topology analysis."
                )
            logger.info(f"Using taxon names: {taxon_names}")
            logger.info(f"Using outgroup: {outgroup}")

            # Infer ploidy from sample names
            # ploidy = infer_ploidy_from_samples(tree_data[0])  # Pass first tree
            # print(f"Detected ploidy: {ploidy}x")

        # Process tree file to generate CSV (this will handle twisst installation)
        logger.info("Converting trees to topology weights...")
        csv_file = process_tree_file(
            tree_file=file,
            taxon_names=taxon_names,
            outgroup=outgroup,
            output_dir=results_dir,
            verbose=False,  # Disable verbose output here since we handle it in run_analysis
            topology_mapping=topology_mapping,
        )
        logger.info(f"Generated topology weights CSV: {csv_file}")

        print("Tree processing complete.")
        # print(f"Using generated CSV: {csv_file}")

    elif file_type == "csv":
        logger.info("Using CSV file directly for analysis")
        print("Detected CSV file format. Using file directly for analysis...")
        csv_file = file

    else:
        logger.error(f"Unsupported file type: {file_type}")
        raise ValueError(f"Unsupported file type: {file_type}")

    # Load and process the CSV data (either original or generated from trees)
    logger.info(f"Loading data from: {csv_file}")
    print(f"Loading data from: {csv_file}")
    data = dump_data(csv_file)
    logger.info(f"Loaded data shape: {data.shape}")
    logger.debug(f"Data columns: {list(data.columns)}")

    # Run standard twisstntern analyses
    logger.info("Running triangle analysis...")
    print("Running triangle analysis...")
    results = triangles_analysis(data, granularity)
    logger.info(f"Triangle analysis completed. Results shape: {results.shape}")

    logger.info("Running fundamental asymmetry analysis...")
    print("Running fundamental asymmetry analysis...")
    fundamental_results = fundamental_asymmetry(data)
    logger.info(f"Fundamental asymmetry: D_LR={fundamental_results[2]:.4f}, p-value={fundamental_results[4]:.4e}")

    # Generate output prefix based on original file name
    output_prefix = str(results_dir / Path(file).stem)
    logger.debug(f"Output prefix: {output_prefix}")

    # Generate all plots
    logger.info("Generating visualizations...")
    print("Generating visualizations...")
    plot_fundamental_asymmetry(data, output_prefix)
    logger.debug("Generated fundamental asymmetry plot")
    plot(data, granularity, output_prefix)
    logger.debug("Generated ternary plot")
    plot_results(results, granularity, output_prefix)
    logger.debug("Generated results plot")
    plotting_triangle_index(granularity, output_prefix)
    logger.debug("Generated triangle index plot")

    # ====================================================================
    # Add main sub-triangle results to the dataframe

    # logger.info("Adding main sub-triangle results to the dataframe...")
    # # adding the main sub-triangle results to the dataframe
    new_row = pd.DataFrame(
        [["main subtriangle", *fundamental_results[:5], "NA"]],
        columns=results.columns,
        index=["full dataset"],
    )
    results = pd.concat([new_row, results])

    ##################################
    # save plots and results
    ##################################
    logger.info("Saving results...")
    results_csv = results_dir / f"{Path(file).stem}_triangle_analysis.csv"
    results.to_csv(results_csv, index=False, float_format="%.3f")
    logger.info(f"Saved triangle analysis results to: {results_csv}")
    print(f"Saved triangle analysis results to: {results_csv}")

    logger.info("Analysis pipeline completed successfully!")
    print("Analysis pipeline completed successfully!")

    return results, fundamental_results, csv_file
