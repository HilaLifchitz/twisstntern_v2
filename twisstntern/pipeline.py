#!/usr/bin/env python
# coding: utf-8

from pathlib import Path
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


def detect_file_type(file_path):
    """
    Detect whether the input file is a tree file or CSV file based on extension.

    Args:
        file_path (str): Path to the input file

    Returns:
        str: "tree" for tree files, "csv" for CSV files
    """
    file_path = Path(file_path)
    tree_extensions = {".trees", ".ts", ".newick", ".nwk", ".nexus"}
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


def process_tree_file(tree_file, taxon_names=None, outgroup=None, output_dir="Results"):
    """
    Process a tree file to generate topology weights CSV file.

    Args:
        tree_file (str): Path to the tree file
        taxon_names (list, optional): List of taxon names for Newick files
        outgroup (str, optional): Outgroup taxon name
        output_dir (str): Directory to save the output CSV file

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

    print(f"Processing tree file: {tree_file}")
    print(f"Output CSV will be saved to: {csv_output}")

    # Process trees using the unified function, which automatically saves csv topology weights
    topology_weights_df = trees_to_twisst_weights_unified(
        file_path=tree_file,
        taxon_names=taxon_names,
        outgroup=outgroup,
        verbose=True,
    )

    print(f"✓ Successfully generated topology weights CSV: {csv_output}")
    print(f"  - Shape: {topology_weights_df.shape}")
    print(f"  - Columns: {list(topology_weights_df.columns)}")

    return str(csv_output)


# default granularity is 0.1
def run_analysis(file, granularity=0.1, taxon_names=None, outgroup=None):
    """
    Orchestrates the full analysis and visualization pipeline for both tree files and CSV files.

    Args:
        file (str): Path to the input file (tree file or CSV file).
        granularity (str or float): Granularity level ("superfine", "fine", "coarse", or a float).
        taxon_names (list, optional): List of taxon names for Newick tree files.
                                     Ignored for TreeSequence files and CSV files.
        outgroup (str, optional): Outgroup taxon name for tree files.
                                 Ignored for CSV files.

    Returns:
        tuple: (results, fundamental_results, csv_file_used)
            - results: Triangle analysis results
            - fundamental_results: Fundamental asymmetry results
            - csv_file_used: Path to the CSV file that was analyzed (original or generated)
    """
    # Ensure Results directory exists
    results_dir = Path("Results")
    results_dir.mkdir(exist_ok=True)

    # Detect file type and process accordingly
    file_type = detect_file_type(file)

    if file_type == "tree":
        print(
            "Detected tree file format. Processing trees to generate topology weights..."
        )

        # if the tree is in Newick format, the users had to provide taxon names
        _, tree_type = detect_and_read_trees(file)
        if tree_type == "newick":
            if taxon_names is None:
                raise ValueError("Taxon names are required for Newick files")
            if outgroup is None:
                raise ValueError("Outgroup is required for Newick files")

        # Process tree file to generate CSV (this will handle twisst installation)
        csv_file = process_tree_file(
            tree_file=file,
            taxon_names=taxon_names,
            outgroup=outgroup,
            output_dir=results_dir,
        )

        print(f"Tree processing complete. Using generated CSV: {csv_file}")

    elif file_type == "csv":
        print("Detected CSV file format. Using file directly for analysis...")
        csv_file = file

    else:
        raise ValueError(f"Unsupported file type: {file_type}")

    # Load and process the CSV data (either original or generated from trees)
    print(f"Loading data from: {csv_file}")
    data = dump_data(csv_file)

    # Run standard twisstntern analyses
    print("Running triangle analysis...")
    results = triangles_analysis(data, granularity)

    print("Running fundamental asymmetry analysis...")
    fundamental_results = fundamental_asymmetry(data)

    # Generate output prefix based on original file name
    output_prefix = str(results_dir / Path(file).stem)

    # Generate all plots
    print("Generating visualizations...")
    plot_fundamental_asymmetry(data, output_prefix)
    plot(data, granularity, output_prefix)
    plot_results(results, granularity, output_prefix)
    plotting_triangle_index(granularity, output_prefix)

    # Save results as CSV
    results_csv = results_dir / f"{Path(file).stem}_triangle_analysis.csv"
    results.to_csv(results_csv, index=False)
    print(f"Saved triangle analysis results to: {results_csv}")

    print("Analysis pipeline completed successfully!")

    return results, fundamental_results, csv_file
