"""
Tree Processing Module for TWISSTNTERN

This module provides comprehensive functionality for processing phylogenetic trees from various formats
and converting them to topology weights using the twisst method. Key capabilities include:

1. **Format Detection & Reading**: Automatically detects and reads TreeSequence (.ts), Newick (.newick), 
   and Nexus (.nexus) format files
2. **Topology Weight Generation**: Uses twisst to generate topology weights for phylogenetic inference
3. **User-Defined Topology Mapping**: Allows users to specify custom topology ordering (T1, T2, T3)
4. **Generator Support**: Efficiently handles both single TreeSequence objects and generators of multiple TreeSequences
5. **Beautiful Visualization**: Displays topology trees using ASCII art via twisst's built-in rendering

The module supports both TreeSequence data (from msprime/tskit) and Newick format trees, with automatic
format detection and appropriate processing pipelines for each.

Main Functions:
- trees_to_twisst_weights_unified(): Main entry point for processing any tree file format
- ts_to_twisst_weights(): Process TreeSequence objects (single or generator)
- newick_to_twisst_weights(): Process Newick format trees
- parse_topology_mapping(): Parse user topology preferences
- reorder_weights_by_topology_preference(): Apply custom topology ordering
"""

import os
import sys
import tempfile
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Union, List, Tuple, Optional
import tskit
import ete3
import msprime
import re
import logging


# Import twisst functions from core external directory
from ..core.external.twisst import weightTrees, summary

# Import logger
from ..core.logger import get_logger, log_topologies


# to know if the file is a tree sequence or a newick file
def detect_and_read_trees(
    file_path: str,
) -> Tuple[Optional[Union[tskit.TreeSequence, List[str]]], str]:
    """
    Automatically detect tree file format and load the appropriate data structure.

    This function serves as the main entry point for reading tree data from files.
    It supports multiple formats and handles format-specific parsing requirements.

    Supported Formats:
    - TreeSequence (.trees, .ts): Binary TSKit tree sequence files from msprime/tskit
    - Newick (.newick, .nwk, .tree): Standard Newick format phylogenetic trees
    - Nexus (.nexus): Nexus format files containing tree definitions

    Args:
        file_path (str): Absolute or relative path to the tree file

    Returns:
        Tuple containing:
        - Tree data: Either a tskit.TreeSequence object or List[str] of Newick strings
        - Format label: String indicating detected format ('ts' or 'newick')

    Raises:
        RuntimeError: If TreeSequence file cannot be loaded
        ValueError: If Newick/Nexus parsing fails or tree format is invalid
        IOError: If file cannot be read
    """
    path = Path(file_path)

    # Case 1: TreeSequence files (.trees or .ts)
    if path.suffix in [".trees", ".ts"]:
        try:
            ts = tskit.load(file_path)
            print("✅ Detected format: TreeSequence (.ts/.trees)")
            return ts, "ts"
        except Exception as e:
            raise RuntimeError(f"Failed to load TreeSequence: {e}")

    # Case 2: Newick files by extension
    if path.suffix in [".newick", ".nwk", ".tree"]:
        try:
            with open(file_path, "r") as f:
                newicks = [line.strip() for line in f if line.strip()]
                print(f"✅ Detected format: Newick ({path.suffix})")
                return newicks, "newick"
        except Exception as e:
            raise ValueError(
                f"Failed to read Newick trees from {path.suffix} file: {e}"
            )

    # Case 3: Content-based detection for files without clear extensions
    try:
        with open(file_path, "r") as f:
            first_line = f.readline()
            if first_line is not None:
                first_line = first_line.strip()
            else:
                first_line = ""
            rest = f.read()
    except Exception as e:
        raise IOError(f"Could not read file: {e}")

    # Case 3a: Nexus format detection
    if first_line.upper().startswith("#NEXUS"):
        try:
            with open(file_path, "r") as f:
                lines = f.readlines()

            # Extract tree definitions from Nexus file (lines starting with "tree")
            newicks = []
            for line in lines:
                if line.strip().lower().startswith("tree"):
                    # Parse tree definition: "tree name = (newick_string);"
                    newick = line.split("=", 1)[-1].strip()
                    try:
                        # Validate Newick string by parsing with ete3
                        t = ete3.Tree(newick)
                        newicks.append(t.write(format=0).strip())
                    except Exception as e:
                        raise ValueError(f"Invalid Newick tree in Nexus file: {e}")

            print("✅ Detected format: Nexus → Converted to Newick")
            return newicks, "newick"

        except Exception as e:
            raise ValueError(f"Failed to parse Nexus file: {e}")

    # Case 3b: Default to Newick format for unknown extensions
    try:
        with open(file_path, "r") as f:
            newicks = [line.strip() for line in f if line.strip()]
            print("✅ Detected format: Newick (unknown extension, treating as Newick)")
            return newicks, "newick"
    except Exception as e:
        raise ValueError(f"Failed to read Newick trees: {e}")


def simplify_topologies(weightsData):
    """
    Convert twisst topology trees to simplified Newick strings for comparison and display.

    This function processes the topology trees returned by twisst.weightTrees and converts
    them to standardized Newick strings suitable for:
    - CSV file headers
    - Topology comparison and matching
    - User-friendly display

    The simplification process:
    1. Removes all branch lengths (sets to 0.0)
    2. Removes internal node names/labels
    3. Preserves leaf names exactly as they appear (population IDs or names)
    4. Uses minimal Newick format (format=9 in ete3)

    Args:
        weightsData (dict): Dictionary returned by twisst.weightTrees containing 'topos' key
                           with ete3.Tree objects representing each topology

    Returns:
        List[str]: Simplified Newick strings, one for each topology in weightsData['topos']
                  Example: ["(0,(1,(2,3)));", "(0,(2,(1,3)));", "(0,(3,(1,2)));"]

    Note:
        This function preserves the exact leaf names from the original trees, which may be
        population IDs (like "0", "1", "2") or population names (like "P1", "P2", "P3").
    """
    simplified_topos = []
    topos = weightsData.get("topos", [])

    for i, tree in enumerate(topos):
        # Clone tree to avoid modifying original
        newick = tree.write(format=1)  # includes branch lengths
        t = ete3.Tree(newick, format=1)
        for node in t.traverse():
            node.dist = 0.0  # remove branch lengths
            if not node.is_leaf():
                node.name = ""  # remove internal node names
            # Keep leaf names (population names like p1, p2, p3) unchanged
        simplified = t.write(format=9)  # minimal Newick, preserves leaf names
        simplified_topos.append(simplified)

    return simplified_topos


def ts_chromosome_to_twisst_weights(
    ts,
    outgroup=None,
    output_file=None,
    verbose=False,
    twisst_verbose=False,
    topology_mapping=None,
):
    """
    Extract topology weights from any TreeSequence object using twisst.

    Args:
        ts (tskit.TreeSequence): Input TreeSequence object
        outgroup (str, optional): Population ID to use as outgroup. If None, uses the first population.
        output_file (str, optional): Path to save CSV file. If None, returns DataFrame.
        verbose (bool): Whether to print verbose output

    Returns:
        pd.DataFrame: Normalized weights ready for dump_data function, first row is the simplified topologies, the rest are the weights
    """

    # Debug info about the TreeSequence
    if verbose:
        print("=== TreeSequence Info ===")
        print(f"Number of populations: {ts.num_populations}")
        print(f"Number of samples: {ts.num_samples}")
        print(f"Number of trees: {ts.num_trees}")

        # Show samples by population
        print("\nSamples by population:")
        for pop_id in range(ts.num_populations):
            samples = [s for s in ts.samples() if ts.node(s).population == pop_id]
            print(f"  Population {pop_id}: {len(samples)} samples")

    # Get populations that actually have samples - THIS IS THE KEY FIX
    populations_with_samples = []
    for pop_id in range(ts.num_populations):
        samples = [s for s in ts.samples() if ts.node(s).population == pop_id]
        if len(samples) > 0:
            populations_with_samples.append(str(pop_id))

    if len(populations_with_samples) < 3:
        raise ValueError(
            f"Need at least 3 populations with samples for topology analysis. Found {len(populations_with_samples)}"
        )

    # Set default outgroup to first population if not specified
    if outgroup is None:
        outgroup = populations_with_samples[0]
        if verbose:
            print(f"Using population {outgroup} as outgroup")

    # Validate outgroup
    if outgroup not in populations_with_samples:
        raise ValueError(
            f"Outgroup population {outgroup} has no samples. Available populations: {populations_with_samples}"
        )

    if verbose:
        print(
            f"Analyzing {len(populations_with_samples)} populations: {populations_with_samples}"
        )
        print(f"Using outgroup: {outgroup}")

    # CRITICAL: We MUST specify taxonNames to only include populations with samples
    # Otherwise twisst will try to analyze all populations, including empty ancestral ones
    weightsData = weightTrees(
        ts,
        treeFormat="ts",
        taxonNames=populations_with_samples,  # EXPLICITLY specify which populations to use
        outgroup=outgroup,
        verbose=twisst_verbose,
    )

    # Apply topology reordering if user preference is provided
    if topology_mapping is not None:
        if verbose:
            print("\nApplying user topology preferences...")
        (
            weights_norm,
            simplified_topos,
            columns,
        ) = reorder_weights_by_topology_preference(weightsData, topology_mapping)
    else:  # if no topology mapping
        # we just print the default topologies order by twisst
        print("No topology mapping was provided; displaying the default topology axis")
        topos = weightsData["topos"]
        for i, topo in enumerate(topos):
            print(f"Topology {i+1}")
            print(topo)

        # Extract normalized weights
        if "weights_norm" in weightsData:
            weights_norm = weightsData["weights_norm"]
        else:
            # Normalize manually if not already normalized
            weights = weightsData["weights"]
            # Check if all weights are zero (which would cause division by zero)
            row_sums = weights.sum(axis=1)
            if np.all(row_sums == 0):
                raise ValueError(
                    "All topology weights are zero - this suggests a problem with the analysis"
                )
            weights_norm = weights / row_sums[:, np.newaxis]

        # Get simplified topologies for CSV header
        simplified_topos = simplify_topologies(weightsData)

        # Create column names based on number of topologies
        n_topos = weights_norm.shape[1]
        if n_topos == 3:
            columns = ["T1", "T2", "T3"]  # Standard 3-topology case (4 populations)
        else:
            columns = [f"Topo{i+1}" for i in range(n_topos)]

        # Log topologies to file
        logger = get_logger(__name__)
        log_topologies(topos, simplified_topos, columns, logger, "TreeSequence topologies (default order)")

    # Get number of topologies for reporting (works for both cases)
    n_topos = len(columns)

    # Create DataFrame
    df = pd.DataFrame(weights_norm, columns=columns)

    # Remove any rows with NaN values (trees where no valid topology was found)
    df = df.dropna()

    if verbose:
        print(f"\nExtracted {len(df)} valid trees with {n_topos} topologies")
        if len(df) > 0:
            print(f"Weight summary:")
            print(df.describe())
        else:
            print("WARNING: No valid topology weights found!")

    # Save to file if requested
    if output_file:
        # Create a new DataFrame with the simplified topologies as header row
        header_row = pd.DataFrame([simplified_topos], columns=columns)

        # Concatenate the header and the weights
        full_df = pd.concat([header_row, df], ignore_index=True)

        # Save
        full_df.to_csv(output_file, index=False, float_format="%.3f")

        if verbose:
            print(f"Saved results to: {output_file}")

    return df


def ts_to_twisst_weights(
    input_data,
    outgroup=None,
    output_file=None,
    verbose=False,
    twisst_verbose=False,
    topology_mapping=None,
):
    """
    Enhanced version of ts_to_twisst_weights that can handle both single TreeSequence objects
    and generators that yield multiple TreeSequence objects.

    Args:
        input_data: Either a single tskit.TreeSequence object or a generator yielding TreeSequence objects
        outgroup (str, optional): Population ID to use as outgroup. If None, uses the first population.
        output_file (str, optional): Path to save CSV file. If None, returns DataFrame.
        verbose (bool): Whether to print verbose output
        twisst_verbose (bool): Whether to print verbose twisst output

    Returns:
        pd.DataFrame: Normalized weights ready for analysis, with consistent topology ordering
    """

    # Check if input is a generator or single TreeSequence
    is_generator = hasattr(input_data, "__iter__") and not isinstance(
        input_data, tskit.TreeSequence
    )

    if not is_generator:
        # Single TreeSequence - use original function
        if verbose:
            print("Processing single TreeSequence...")
        return ts_chromosome_to_twisst_weights(
            input_data,
            outgroup=outgroup,
            output_file=output_file,
            verbose=verbose,
            twisst_verbose=twisst_verbose,
            topology_mapping=topology_mapping,
        )
    #########################################################
    # Generator case - process multiple TreeSequences
    input_data = list(input_data)

    all_weights = []
    canonical_topologies = None
    canonical_simplified_topos = None
    columns = None
    total_processed = 0

    if verbose:
        print(
            f"processing generator of TreeSequences with {len(input_data)} TreeSequences, this may take a while..."
        )

    for i in range(len(input_data)):
        ts = input_data[i]
        # if verbose:
        #     print(f"\nProcessing TreeSequence {i+1}...")

        # # Debug info about the TreeSequence
        # if verbose:
        #     print(f"  Number of populations: {ts.num_populations}")
        #     print(f"  Number of samples: {ts.num_samples}")
        #     print(f"  Number of trees: {ts.num_trees}")

        # Get populations that actually have samples
        populations_with_samples = []
        for pop_id in range(ts.num_populations):
            samples = [s for s in ts.samples() if ts.node(s).population == pop_id]
            if len(samples) > 0:
                populations_with_samples.append(str(pop_id))

        # if len(populations_with_samples) < 3:
            # if verbose:
            #     print(
            #         f"  Skipping: only {len(populations_with_samples)} populations with samples"
            #     )
            # continue

        # Set default outgroup to first population if not specified
        current_outgroup = outgroup
        if current_outgroup is None:
            current_outgroup = populations_with_samples[0]

        # # Validate outgroup
        # if current_outgroup not in populations_with_samples:
        #     if verbose:
        #         print(
        #             f"  Skipping: outgroup {current_outgroup} not in populations with samples"
        #         )
        #     continue

        # if verbose:
        #     print(
        #         f"  Analyzing {len(populations_with_samples)} populations: {populations_with_samples}"
        #     )
        #     print(f"  Using outgroup: {current_outgroup}")

        # Run twisst analysis
        weightsData = weightTrees(
            ts,
            treeFormat="ts",
            taxonNames=populations_with_samples,
            outgroup=current_outgroup,
            verbose=twisst_verbose,
        )

        # Extract simplified topologies for consistency checking
        current_simplified_topos = simplify_topologies(weightsData)

        # special treetments- because its many trees- making ure the topologies order is consistent

        # On first iteration, establish canonical topology order
        if canonical_topologies is None:
            canonical_topologies = weightsData["topos"]
            canonical_simplified_topos = current_simplified_topos
            n_topos = len(canonical_topologies)

            # Create column names
            if n_topos == 3:
                columns = ["T1", "T2", "T3"]
            else:
                columns = [f"Topo{i+1}" for i in range(n_topos)]

        # Map current topologies to canonical order
        topo_mapping = []
        for canonical_topo in canonical_simplified_topos:
            try:
                current_idx = current_simplified_topos.index(canonical_topo)
                topo_mapping.append(current_idx)
            except ValueError:
                raise ValueError(
                    f"TreeSequence {i+1} has different topologies than expected. "
                    f"Expected: {canonical_simplified_topos}, Got: {current_simplified_topos}"
                )

        if verbose and len(topo_mapping) != len(canonical_simplified_topos):
            print(f"  Warning: topology mapping incomplete: {topo_mapping}")

        # Extract and reorder weights according to canonical topology order
        if "weights_norm" in weightsData:
            weights_norm = weightsData["weights_norm"]
        else:
            weights = weightsData["weights"]
            row_sums = weights.sum(axis=1)
            if np.all(row_sums == 0):
                if verbose:
                    print(f"  Skipping: all topology weights are zero")
                continue
            weights_norm = weights / row_sums[:, np.newaxis]

        # Reorder weights to match canonical topology order
        reordered_weights = weights_norm[:, topo_mapping]

        # Create DataFrame and remove NaN rows
        df_current = pd.DataFrame(reordered_weights, columns=columns)
        df_current = df_current.dropna()

        if len(df_current) > 0:
            all_weights.append(df_current)
            total_processed += len(df_current)

        #     if verbose:
        #         print(f"  Extracted {len(df_current)} valid trees")
        # else:
        #     if verbose:
        #         print(f"  No valid trees found")

    # Combine all weights into single DataFrame
    if not all_weights:
        raise ValueError("No valid topology weights found across all TreeSequences")

    combined_df = pd.concat(all_weights, ignore_index=True)
    combined_df = combined_df.round(3)
    combined_df = combined_df.loc[combined_df.iloc[:, 1] != combined_df.iloc[:, 2]]

    print(f"\nWeights with T2 == T3 removed. Remaining rows: {len(combined_df)}")

    # Apply user topology mapping if provided
    if topology_mapping is not None:
        if verbose:
            print("\nApplying user topology preferences to combined data...")

        # Create a temporary weightsData structure for reordering
        temp_weightsData = {
            "topos": canonical_topologies,
            "weights_norm": combined_df.values,
        }

        # Apply reordering (this function also prints the topologies)
        (
            reordered_weights,
            reordered_simplified_topos,
            new_columns,
        ) = reorder_weights_by_topology_preference(temp_weightsData, topology_mapping)

        # Update combined_df with reordered data
        combined_df = pd.DataFrame(reordered_weights, columns=new_columns)
        canonical_simplified_topos = reordered_simplified_topos
        columns = new_columns

    else:  # if no topology mapping, we just print the topologies
        print("No topology mapping was provided; displaying the default topology axis")
        for i, topo in enumerate(canonical_topologies):
            print(f"Topology {i+1}")
            print(topo)
        
        # Log topologies to file
        logger = get_logger(__name__)
        log_topologies(canonical_topologies, canonical_simplified_topos, columns, logger, "Multi-TreeSequence canonical topologies (default order)")

    if verbose:
        print(f"\nCombined Results:")
        print(f"  Total TreeSequences processed: {len(all_weights)}")
        print(f"  Total valid trees: {total_processed}")
        print(f"  Final DataFrame shape: {combined_df.shape}")
        print(f"  Weight summary:")
        print(combined_df.describe())

    # Save to file if requested
    if output_file:
        # Create header row with canonical simplified topologies
        header_row = pd.DataFrame([canonical_simplified_topos], columns=columns)

        # Concatenate header and weights
        full_df = pd.concat([header_row, combined_df], ignore_index=True)

        # Save to CSV
        full_df.to_csv(output_file, index=False, float_format="%.3f")

        if verbose:
            print(f"  Saved results to: {output_file}")

    return combined_df


def newick_to_twisst_weights(
    newick_trees,
    taxon_names,
    outgroup=None,
    output_file=None,
    verbose=False,
    twisst_verbose=False,
    topology_mapping=None,
):
    """
    Extract topology weights from Newick tree strings using twisst.

    Args:
        newick_trees (List[str] or str): List of Newick tree strings or single Newick string
        taxon_names (List[str]): List of population names. .
        outgroup (str, optional): Population name to use as outgroup. If None, uses the first population.
        output_file (str, optional): Path to save CSV file. If None, returns DataFrame.
        verbose (bool): Whether to print verbose output
        twisst_verbose (bool): Whether to print verbose twisst output
        topology_mapping (dict or str, optional): User preference for topology ordering

    Returns:
        pd.DataFrame: Normalized weights ready for analysis
    """

    # Handle single tree string
    if isinstance(newick_trees, str):
        newick_trees = [newick_trees]

    if len(newick_trees) == 0:
        raise ValueError("No Newick trees provided")

    if verbose:
        print("=== Newick Trees Info ===")
        print(f"Number of trees: {len(newick_trees)}")

    # Extract all sample names from first tree
    try:
        first_tree = ete3.Tree(newick_trees[0])
        all_sample_names = [leaf.name for leaf in first_tree.get_leaves()]
        if verbose:
            print(f"Found {len(all_sample_names)} samples in first tree")
            print(
                f"Sample names: {all_sample_names[:10]}{'...' if len(all_sample_names) > 10 else ''}"
            )
    except Exception as e:
        raise ValueError(f"Failed to parse first Newick tree: {e}")

    # Group samples by population (infer from sample names like 'P1_1', 'P2_5', etc.)
    population_samples = {}
    population_names = []

    for sample in all_sample_names:
        # Extract population name (part before underscore)
        if "_" in sample:
            pop_name = sample.split("_")[0]
        else:
            # If no underscore, treat each sample as its own population
            pop_name = sample

        if pop_name not in population_samples:
            population_samples[pop_name] = []
            population_names.append(pop_name)

        population_samples[pop_name].append(sample)

    if verbose:
        print(f"Detected {len(population_names)} populations: {population_names}")
        for pop_name in population_names:
            sample_count = len(population_samples[pop_name])
            print(f"  {pop_name}: {sample_count} samples")

    # Override with user-provided taxon_names if provided
    if taxon_names is not None:
        # Validate that user-provided taxon_names match detected populations
        if set(taxon_names) != set(population_names):
            if verbose:
                print(
                    f"Warning: User-provided taxon_names {taxon_names} don't match detected populations {population_names}"
                )
                print("Using user-provided taxon_names anyway...")
            population_names = taxon_names
        else:
            population_names = taxon_names

    if len(population_names) < 3:
        raise ValueError(
            f"Need at least 3 populations for topology analysis. Found {len(population_names)}: {population_names}"
        )

    # Set default outgroup if not specified
    if outgroup is None:
        outgroup = population_names[0]
        if verbose:
            print(f"Using {outgroup} as outgroup")

    # Validate outgroup
    if outgroup not in population_names:
        raise ValueError(
            f"Outgroup '{outgroup}' not found in populations: {population_names}"
        )

    if verbose:
        print(f"Analyzing {len(population_names)} populations: {population_names}")
        print(f"Using outgroup: {outgroup}")

    # Validate all trees have the same sample set
    expected_samples = set(all_sample_names)
    for i, newick in enumerate(newick_trees):
        try:
            tree = ete3.Tree(newick)
            tree_samples = set(leaf.name for leaf in tree.get_leaves())

            if tree_samples != expected_samples:
                raise ValueError(
                    f"Tree {i+1} has different samples than expected. "
                    f"Expected {len(expected_samples)} samples, got {len(tree_samples)} samples"
                )
        except Exception as e:
            raise ValueError(f"Invalid Newick tree at position {i+1}: {e}")

    # Convert Newick strings to ETE3 tree objects
    try:
        ete_trees = [ete3.Tree(newick) for newick in newick_trees]
        if verbose:
            print(f"Successfully parsed {len(ete_trees)} Newick trees to ETE3 objects")
    except Exception as e:
        raise ValueError(f"Failed to parse Newick trees to ETE3 objects: {e}")

    # Create taxa structure for twisst (list of lists, one per population)
    taxa = []
    for pop_name in population_names:
        taxa.append(population_samples[pop_name])

    # Run twisst analysis
    if verbose:
        print("Running twisst analysis on Newick trees...")
        print(f"Taxa structure for twisst:")
        for i, (pop_name, samples) in enumerate(zip(population_names, taxa)):
            print(f"  {pop_name}: {len(samples)} samples")

    weightsData = weightTrees(
        ete_trees,  # Pass ETE3 tree objects
        treeFormat="ete3",  # Use ete3 format
        taxa=taxa,  # Dynamic taxa structure
        taxonNames=taxon_names,  # Dynamic population names
        outgroup=outgroup,
        verbose=twisst_verbose,
    )

    # explicatly printing the topologies--COME BACK TO THIS- WE WANNA PRINT THIS AFTER THE REORDERING!!
    # topos = weightsData["topos"]
    # for i, topo in enumerate(topos):
    # print(f"Topology {i+1} (T{i+1})")
    # print(topo)

    # Apply topology reordering if user preference is provided
    if topology_mapping is not None:
        print("\nApplying user topology preferences...")
        (
            weights_norm,
            simplified_topos,
            columns,
        ) = reorder_weights_by_topology_preference(weightsData, topology_mapping)
    else:
        # Use default ordering - always print topologies for Newick case
        print("No topology mapping was provided; displaying the default topology axis")
        topos = weightsData["topos"]
        for i, topo in enumerate(topos):
            print(f"Topology {i+1}")
            print(topo)

        # Extract normalized weights
        if "weights_norm" in weightsData:
            weights_norm = weightsData["weights_norm"]
        else:
            # Normalize manually if not already normalized
            weights = weightsData["weights"]
            # Check if all weights are zero (which would cause division by zero)
            row_sums = weights.sum(axis=1)
            if np.all(row_sums == 0):
                raise ValueError(
                    "All topology weights are zero - this suggests a problem with the analysis"
                )
            weights_norm = weights / row_sums[:, np.newaxis]

        # Get simplified topologies for CSV header
        simplified_topos = simplify_topologies(weightsData)
        
        # Create column names based on number of topologies
        n_topos = weights_norm.shape[1]
        if n_topos == 3:
            columns = ["T1", "T2", "T3"]  # Standard 3-topology case (4 populations)
        else:
            columns = [f"Topo{i+1}" for i in range(n_topos)]
        
        # Log topologies to file
        logger = get_logger(__name__)
        log_topologies(topos, simplified_topos, columns, logger, "Newick file topologies (default order)")

    # Create DataFrame
    df = pd.DataFrame(weights_norm, columns=columns)

    # Remove any rows with NaN values (trees where no valid topology was found)
    df = df.dropna()

    # Get number of topologies for reporting
    n_topos = len(columns)

    if verbose:
        print(f"\nExtracted {len(df)} valid trees with {n_topos} topologies")
        if len(df) > 0:
            print(f"Weight summary:")
            print(df.describe())
        else:
            print("WARNING: No valid topology weights found!")

    # Save to file if requested
    if output_file:
        # Use the simplified_topos that was already computed above
        # Create a new DataFrame with that single row
        header_row = pd.DataFrame([simplified_topos], columns=columns)

        # Concatenate the header and the weights
        full_df = pd.concat([header_row, df], ignore_index=True)

        # Save
        full_df.to_csv(output_file, index=False, float_format="%.3f")

        if verbose:
            print(f"Saved results to: {output_file}")

    return df


# the main function(!!!): detecting the format and routing to the appropriate function
def trees_to_twisst_weights_unified(
    file_path,
    taxon_names=None,
    outgroup=None,
    output_file=None,
    verbose=False,
    twisst_verbose=False,
    topology_mapping=None,
):
    """
    Unified function to extract topology weights from any tree file format.
    """
    if verbose:
        print(f"Processing tree file: {file_path}")

    # If no output file specified, create default one in Results directory
    if output_file is None:
        results_dir = Path("Results")
        results_dir.mkdir(exist_ok=True)
        input_filename = Path(file_path).stem  # filename without extension
        output_file = results_dir / f"{input_filename}_topology_weights.csv"

    if verbose:
        print(f"Will save results to: {output_file}")

    # Detect format and load trees
    try:
        tree_data, format_type = detect_and_read_trees(file_path)
    except Exception as e:
        raise ValueError(f"Failed to read tree file '{file_path}': {e}")

    # Route to appropriate processing function based on detected format
    if format_type == "ts":
        return ts_to_twisst_weights(
            tree_data,
            outgroup=outgroup,
            output_file=str(output_file),
            verbose=verbose,
            twisst_verbose=twisst_verbose,
            topology_mapping=topology_mapping,
        )
    elif format_type == "newick":
        return newick_to_twisst_weights(
            tree_data,
            taxon_names=taxon_names,
            outgroup=outgroup,
            output_file=str(output_file),
            verbose=verbose,
            twisst_verbose=twisst_verbose,
            topology_mapping=topology_mapping,
        )
    else:
        raise ValueError(f"Unsupported tree format: {format_type}")


# Not used directly, but still good to have:
# taking in a ts object- not in generator or Newick format- and returning the populations with samples
# Utility functions for debugging and backward compatibility


###########################################################
# reordering topologies- detmining which is T1, T2 and T3
###########################################################


# turns the input string of topologies order from the user into a dictionary
def parse_topology_mapping(mapping_string):
    """
    Parse user-specified topology mapping from a string like:
    "T1=(0,(1,(2,3))); T2=(0,(2,(1,3))); T3=(0,(3,(1,2)));"  
    or (for Newick files)
    "T1=(O,(P1,(P2,P3))); T2=(O,(P2,(P1,P3))); T3=(O,(P3,(P1,P2)));"
    
    Args:
        mapping_string (str): String containing topology assignments
        
    Returns:
        dict: Mapping from T1/T2/T3 to simplified topology strings
    """
    mapping = {}
    
    # Remove extra whitespace and split by semicolon
    assignments = [part.strip() for part in mapping_string.split(";") if part.strip()]
    
    for assignment in assignments:
        if "=" not in assignment:
            continue
            
        # Split on first '=' to handle cases where topology contains '='
        key, value = assignment.split("=", 1)
        key = key.strip()
        value = value.strip()
        
        # Remove surrounding quotes if present
        if value.startswith('"') and value.endswith('"'):
            value = value[1:-1]
        elif value.startswith("'") and value.endswith("'"):
            value = value[1:-1]
            
        # Validate key is T1, T2, or T3
        if key not in ["T1", "T2", "T3"]:
            raise ValueError(f"Invalid topology key: {key}. Must be T1, T2, or T3")
            
        mapping[key] = value
    
    # Ensure all three topologies are specified
    if len(mapping) != 3:
        missing = set(["T1", "T2", "T3"]) - set(mapping.keys())
        raise ValueError(f"Missing topology assignments: {missing}")
    
    # Validate that all three topologies are distinct
    keys = ["T1", "T2", "T3"]
    normalized_values = [normalize_topology_string(mapping[key]) for key in keys]
    
    if len(set(normalized_values)) != 3:
        # Find duplicates
        duplicates = []
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                if normalized_values[i] == normalized_values[j]:
                    duplicates.append(f"{keys[i]}={mapping[keys[i]]} and {keys[j]}={mapping[keys[j]]}")
        
        raise ValueError(
            f"All three topologies must be distinct. Found duplicate topologies: {'; '.join(duplicates)}"
        )
    
    return mapping


# SHOULDNT BE NECESSARY!
# Legacy functions - kept for backward compatibility but not actively used


# so we can compare topologies
def normalize_topology_string(topo_str):
    """
    Normalize topology string for comparison by standardizing population names to lowercase.
    Converts P1 → p1, P2 → p2, etc. (case-insensitive comparison)
    Does NOT convert p1 to 1 - preserves the letter format.
    Also removes trailing semicolons, quotes, and whitespace for consistent comparison.
    
    Args:
        topo_str (str): Topology string like "(0,(p1,(p2,p3)))" or "(0,(P1,(P2,P3)));"
        
    Returns:
        str: Normalized topology string like "(0,(p1,(p2,p3)))"
    """
    import re
    
    # Remove leading/trailing quotes, semicolons and whitespace
    normalized = topo_str.strip().rstrip(";").strip('"').strip("'")
    
    # Convert P1/P2/P3 to p1/p2/p3 (lowercase) for consistent comparison
    # This handles case-insensitivity without converting to digits
    normalized = re.sub(r"P(\d+)", r"p\1", normalized)
    
    return normalized


# needed for reordering topologies
def compare_topologies(topo1, topo2):
    """
    Compare two topology strings, accounting for different population naming conventions.
    
    Args:
        topo1 (str): First topology string
        topo2 (str): Second topology string
        
    Returns:
        bool: True if topologies are equivalent
    """
    return normalize_topology_string(topo1) == normalize_topology_string(topo2)


# The main one here: reordering the weights by the topology mapping and printing the trees
def reorder_weights_by_topology_preference(weightsData, topology_mapping):
    """
    Reorder topology weights according to user preference.
    This should be called BEFORE creating the DataFrame.
    
    Args:
        weightsData (dict): Output from twisst.weightTrees containing 'topos' and 'weights'
        topology_mapping (dict): Mapping from T1/T2/T3 to desired topology strings
        
    Returns:
        tuple: (reordered_weights, reordered_simplified_topos, column_names)
    """
    # Get current topologies and weights
    original_topos = weightsData["topos"]
    # we could print these if we want to see the original topologies
    # for i, topo in enumerate(original_topos):
    #     print(f"Topology {i+1} (T{i+1})")
    #     print(topo)

    original_simplified_topos = simplify_topologies(weightsData)
    
    # Parse topology mapping if it's a string (do this FIRST)
    if isinstance(topology_mapping, str):
        topology_mapping = parse_topology_mapping(topology_mapping)

    # Get weights (normalized or raw)
    if "weights_norm" in weightsData:
        weights = weightsData["weights_norm"]
    else:
        weights = weightsData["weights"]
        # Normalize if needed
        row_sums = weights.sum(axis=1)
        if not np.all(row_sums == 0):
            weights = weights / row_sums[:, np.newaxis]
    
    # Validate that we have exactly matching number of topologies
    num_expected_topologies = len(topology_mapping.keys())
    if len(original_simplified_topos) != num_expected_topologies:
        raise ValueError(
            f"Expected {num_expected_topologies} topologies, found {len(original_simplified_topos)}"
        )
    
    # Find mapping from user preferences to current column positions
    column_mapping = {}  # {target_column: source_column_index}
    
    for target_col in ["T1", "T2", "T3"]:
        desired_topo = topology_mapping[target_col]
        
        # Find which current column contains this topology
        found = False
        for i, current_topo in enumerate(original_simplified_topos):
            if compare_topologies(current_topo, desired_topo):
                column_mapping[target_col] = i
                found = True
                break
        
        if not found:
            print(f"Available topologies:")
            for i, topo in enumerate(original_simplified_topos):
                print(f"  {i}: {topo} (normalized: {normalize_topology_string(topo)})")
            print(
                f"User requested: {desired_topo} (normalized: {normalize_topology_string(desired_topo)})"
            )
            raise ValueError(
                f"Topology '{desired_topo}' not found in computed topologies."
            )
    
    # Create new column order
    new_column_order = [
        column_mapping["T1"],
        column_mapping["T2"],
        column_mapping["T3"],
    ]
    
    # Reorder weights and topologies
    reordered_weights = weights[:, new_column_order]
    reordered_simplified_topos = [
        original_simplified_topos[i] for i in new_column_order
    ]

    # Print beautiful tree representation
    print_topology_mapping_with_trees(weightsData, topology_mapping)

    return reordered_weights, reordered_simplified_topos, ["T1", "T2", "T3"]


# nice printing of the topologies to the terminal
def print_topology_mapping_with_trees(weightsData, topology_mapping):
    """
    Print topology mapping with beautiful ASCII tree representations.
    Uses twisst's built-in tree rendering to show the actual tree structure.

    Args:
        weightsData (dict): Output from twisst.weightTrees containing 'topos'
        topology_mapping (dict): Mapping from T1/T2/T3 to desired topology strings
    """
    # Get original topologies and simplified versions
    original_topos = weightsData[
        "topos"
    ]  # These are the twisst tree objects with .get_ascii()
    original_simplified_topos = simplify_topologies(weightsData)

    # Parse topology mapping if it's a string
    if isinstance(topology_mapping, str):
        topology_mapping = parse_topology_mapping(topology_mapping)

    print("Applied topology mapping:")
    print("=" * 50)

    # Prepare data for logging
    mapped_topos = []
    mapped_simplified = []
    mapped_columns = []
    
    # For each T1, T2, T3, find the corresponding original topology index
    for target_label in ["T1", "T2", "T3"]:
        desired_topo_string = topology_mapping[target_label]

        # Find which original topology matches this string
        original_index = None
        for i, simplified_topo in enumerate(original_simplified_topos):
            if compare_topologies(simplified_topo, desired_topo_string):
                original_index = i
                break

        if original_index is not None:
            print(f"\n{target_label}:")
            print(original_topos[original_index].get_ascii())
            print(f"String: {original_simplified_topos[original_index]}")
            
            # Collect data for logging
            mapped_topos.append(original_topos[original_index])
            mapped_simplified.append(original_simplified_topos[original_index])
            mapped_columns.append(target_label)
        else:
            print(f"\n{target_label}: ERROR - topology not found!")

    print("=" * 50)
    
    # Log the reordered topologies
    if mapped_topos:  # Only log if we have valid mappings
        logger = get_logger(__name__)
        log_topologies(mapped_topos, mapped_simplified, mapped_columns, logger, "Applied topology mapping")
