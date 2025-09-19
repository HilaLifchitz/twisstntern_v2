"""
Module for processing tree sequences and generating topology weights.
This module provides functionality to:
1. Read tree sequences from various formats
2. Generate topology weights using twisst
3. Convert topology weights to our analysis format
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Optional
import tskit
import ete3


# Add the external directory to the Python path
EXTERNAL_DIR = Path(__file__).parent / "external"
sys.path.append(str(EXTERNAL_DIR))

# Import twisst functions directly from local external directory
from twisst import weightTrees

# Import logging from twisstntern
from twisstntern.logger import get_logger
from twisstntern.tree_processing import log_topologies


# exrtacting the toplogies as easy strings
def simplify_topologies(weightsData):
    """
    Takes a twisst weightsData dictionary and returns simplified Newick strings
    for each topology, removing branch lengths and internal node names.

    Args:
        weightsData (dict): Output dictionary from twisst.weightTrees, must contain 'topos'

    Returns:
        List[str]: Simplified Newick strings for each topology
    """
    simplified_topos = []
    topos = weightsData.get("topos", [])

    # Handle case where topos might be None
    if topos is None:
        return []

    for i, tree in enumerate(topos):
        # Clone tree to avoid modifying original
        newick = tree.write(format=1)  # includes branch lengths
        t = ete3.Tree(newick, format=1)
        for node in t.traverse():
            node.dist = 0.0  # remove branch lengths
            if not node.is_leaf():
                node.name = ""  # remove internal node names
        simplified = t.write(format=9)  # minimal Newick
        simplified_topos.append(simplified)

    return simplified_topos


# this is the main function for the chromosome mode
def ts_chromosome_to_twisst_weights(
    ts,
    outgroup=None,
    output_file=None,
    verbose=False,
    twisst_verbose=False,
    topology_mapping=None,
    population_labels=None,
):
    """
    Extract topology weights from any TreeSequence object using twisst.

    Args:
        ts (tskit.TreeSequence): Input TreeSequence object
        outgroup (str, optional): Population ID to use as outgroup. If None, uses the first population.
        output_file (str, optional): Path to save CSV file. If None, returns DataFrame.
        verbose (bool): Whether to print verbose output
        topology_mapping (str or dict, optional): Custom topology mapping for T1/T2/T3

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
        ) = reorder_weights_by_topology_preference(
            weightsData, topology_mapping, population_labels
        )
    else:  # if no topology mapping
        # we just print the default topologies order by twisst
        # print("No topology mapping was provided; displaying the default topology axis")
        # topos = weightsData["topos"]
        # for i, topo in enumerate(topos):
        #     print(f"Topology {i+1}")
        #     print(topo)

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
        topos = weightsData["topos"]
        log_topologies(
            topos,
            simplified_topos,
            columns,
            logger,
            "TreeSequence topologies (default order)",
            population_labels,
        )

    # Get number of topologies for reporting (works for both cases)
    n_topos = len(columns)

    # Create DataFrame
    df = pd.DataFrame(weights_norm, columns=columns)

    # Add position column: start position of each tree
    positions = [tree.interval[0] for tree in ts.trees()]
    df.insert(0, "position", positions)

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


##################################################################################
# main topology weights generating function
#################################################################################


def ts_to_twisst_weights(
    input_data,
    outgroup=None,
    output_file=None,
    verbose=False,
    twisst_verbose=False,
    topology_mapping=None,
    population_labels=None,
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
        topology_mapping (str or dict, optional): Custom topology mapping for T1/T2/T3

    Returns:
        pd.DataFrame: Normalized weights ready for analysis, with consistent topology ordering
    """

    # Check if input is a generator or single TreeSequence
    is_generator = hasattr(input_data, "__iter__") and not isinstance(
        input_data, tskit.TreeSequence
    )

    if not is_generator:
        # Single TreeSequence case (chromosome mode)
        df = ts_chromosome_to_twisst_weights(
            input_data,
            outgroup=outgroup,
            output_file=output_file,
            verbose=verbose,
            twisst_verbose=twisst_verbose,
            topology_mapping=topology_mapping,
            population_labels=population_labels,
        )
        # Ensure 'position' column is present (should already be added by ts_chromosome_to_twisst_weights)
        if "position" not in df.columns:
            if isinstance(input_data, tskit.TreeSequence):
                positions = [tree.interval[0] for tree in input_data.trees()]
                df.insert(0, "position", positions)
        return df
    #########################################################
    # Generator case - process multiple TreeSequences
    # Convert generator to list for processing - but check if it's actually iterable first
    if hasattr(input_data, "__iter__"):
        input_data = list(input_data)
    else:
        # This shouldn't happen given our earlier check, but handle it gracefully
        raise TypeError(
            "Input data is neither a TreeSequence nor an iterable of TreeSequences"
        )

    if verbose:
        print("Processing generator of TreeSequences...")

    all_weights = []
    canonical_topologies = None
    canonical_simplified_topos: Optional[List[str]] = None
    columns: Optional[List[str]] = None
    total_processed = 0

    if verbose:
        print(
            f"processing generator of TreeSequences with {len(input_data)} TreeSequences, this may take a while..."
        )

    for i in range(len(input_data)):
        ts = input_data[i]
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

        # Validate outgroup
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

            # Log canonical topologies
            logger = get_logger(__name__)
            log_topologies(
                canonical_topologies,
                canonical_simplified_topos,
                columns,
                logger,
                "Multi-TreeSequence canonical topologies (default order)",
                population_labels,
            )

        # Type assertion for Pylance
        assert canonical_simplified_topos is not None
        assert columns is not None

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

    # Type assertion for Pylance
    assert canonical_simplified_topos is not None
    assert columns is not None

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
        ) = reorder_weights_by_topology_preference(
            temp_weightsData, topology_mapping, population_labels
        )

        # Update combined_df with reordered data
        combined_df = pd.DataFrame(reordered_weights, columns=new_columns)
        canonical_simplified_topos = reordered_simplified_topos
        columns = new_columns

    else:  # if no topology mapping, we just print the topologies
        print("No topology mapping was provided; displaying the default topology axis")
        for i, topo in enumerate(canonical_topologies):
            print(f"Topology {i+1}")
            print(topo)

        # Topologies already logged during first iteration, no need to log again

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


# Not used directly, but still good to have:
# taking in a ts object- not in generator or Newick format- and returning the populations with samples
def debug_ts_populations(ts):
    """
    Debug function to inspect TreeSequence population structure.
    Useful for understanding what populations are available.
    """
    print("=== TreeSequence Population Debug ===")
    print(f"Total populations: {ts.num_populations}")
    print(f"Total samples: {ts.num_samples}")
    print(f"Total trees: {ts.num_trees}")
    print()

    # Show population details
    print("Population details:")
    for i in range(ts.num_populations):
        pop = ts.population(i)
        metadata = pop.metadata if pop.metadata else {}
        name = metadata.get("name", f"Pop{i}") if metadata else f"Pop{i}"

        # Count samples in this population
        samples = [int(s) for s in ts.samples() if ts.node(s).population == i]

        print(f"  Population {i}: name='{name}', {len(samples)} samples")
        if len(samples) > 0:
            print(f"    Sample IDs: {samples[:5]}{'...' if len(samples) > 5 else ''}")

    print()

    # Show what twisst would extract
    populations_with_samples = []
    for pop_id in range(ts.num_populations):
        samples = [s for s in ts.samples() if ts.node(s).population == pop_id]
        if len(samples) > 0:
            populations_with_samples.append(str(pop_id))

    print(f"Populations with samples (for twisst): {populations_with_samples}")

    if len(populations_with_samples) >= 3:
        print(
            f"✓ Ready for topology analysis ({len(populations_with_samples)} populations)"
        )
    else:
        print(
            f"✗ Need at least 3 populations with samples (found {len(populations_with_samples)})"
        )

    return populations_with_samples


######################################################################################
# TOPOLOGY MAPPING FUNCTIONS (copied from twisstntern)
######################################################################################


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
                    duplicates.append(
                        f"{keys[i]}={mapping[keys[i]]} and {keys[j]}={mapping[keys[j]]}"
                    )

        raise ValueError(
            f"All three topologies must be distinct. Found duplicate topologies: {'; '.join(duplicates)}"
        )

    return mapping


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


def reorder_weights_by_topology_preference(
    weightsData, topology_mapping, population_labels=None
):
    """
    Reorder topology weights according to user preference.
    This should be called BEFORE creating the DataFrame.

    Args:
        weightsData (dict): Output from twisst.weightTrees containing 'topos' and 'weights'
        topology_mapping (dict): Mapping from T1/T2/T3 to desired topology strings
        population_labels (dict, optional): Mapping from population IDs to descriptive labels

    Returns:
        tuple: (reordered_weights, reordered_simplified_topos, column_names)
    """
    # Get current topologies and weights
    original_topos = weightsData["topos"]
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

    # Print beautiful tree representation (this function also handles logging)
    print_topology_mapping_with_trees(weightsData, topology_mapping, population_labels)

    return reordered_weights, reordered_simplified_topos, ["T1", "T2", "T3"]


def print_topology_mapping_with_trees(
    weightsData, topology_mapping, population_labels=None
):
    """
    Print topology mapping with beautiful ASCII tree representations.
    Uses twisst's built-in tree rendering to show the actual tree structure.

    Args:
        weightsData (dict): Output from twisst.weightTrees containing 'topos'
        topology_mapping (dict): Mapping from T1/T2/T3 to desired topology strings
        population_labels (dict, optional): Mapping from population IDs to descriptive labels
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

    # Collect topology information for logging
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

            # Collect for logging
            mapped_topos.append(original_topos[original_index])
            mapped_simplified.append(original_simplified_topos[original_index])
            mapped_columns.append(target_label)
        else:
            print(f"\n{target_label}: ERROR - topology not found!")

    print("=" * 50)

    # Log the applied topology mapping
    if mapped_topos:
        logger = get_logger(__name__)
        log_topologies(
            mapped_topos,
            mapped_simplified,
            mapped_columns,
            logger,
            "Applied topology mapping",
            population_labels,
        )
