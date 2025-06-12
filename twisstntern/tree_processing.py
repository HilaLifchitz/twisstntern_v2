"""
Module for processing tree sequences and generating topology weights.
This module provides functionality to:
1. Read tree sequences from various formats
2. Generate topology weights using twisst
3. Convert topology weights to our analysis format
"""

import os
import sys
import tempfile
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Union, List, Tuple
import tskit
import ete3
import msprime

####################### ts on the spot ##########################

# # KEEP THIS CONSTANT -- SHOULD BE KEPT SMALL FOR TESTING
# nP1 = 10  # number of samples in population 1
# nP2 = 10  # number of samples in population 2
# nP3 = 10  # number of samples in population 3
# n0 = 10  # number of samples in the outgroup

# # Provide divergence times
# # WE KEEP CONSTANT FOR THIS RUN
# t1 = 100  # time of split between population 1 and population 2
# t2 = 200  # time of split between population (1,2) and population 3
# t3 = 300  # time of split between population (1,2,3) and the outgroup 0

# # Provide migration rate
# m = 0.0001  # migration rate between population 2 and population 3
# mig_rate = m

# Ne = 1000  # population size we keep constant for all populations for simplicity
# # Provide population sizes
# NeP1 = Ne
# NeP2 = Ne
# NeP3 = Ne
# NeO = Ne
# NeP12 = Ne
# NeP123 = Ne
# NeANC = Ne


# # A demography where (1,2) first coalesce, with our given Ne

# demography = msprime.Demography()

# # initializing populations
# demography.add_population(name="O", initial_size=NeO)
# demography.add_population(name="P1", initial_size=NeP1)
# demography.add_population(name="P2", initial_size=NeP2)
# demography.add_population(name="P3", initial_size=NeP3)
# demography.add_population(name="P13", initial_size=NeP12)
# demography.add_population(name="P123", initial_size=NeP123)
# demography.add_population(name="ANC", initial_size=NeANC)

# # adding split times
# demography.add_population_split(time=t1, derived=["P1", "P2"], ancestral="P13")
# demography.add_population_split(time=t2, derived=["P13", "P3"], ancestral="P123")
# demography.add_population_split(time=t3, derived=["P123", "O"], ancestral="ANC")

# # setting up gene flow
# demography.set_migration_rate("P2", "P3", mig_rate)

# ploidy=2
# # from collections import defaultdict
# num_replicates = 20
# COMMENTED OUT: These were executing at module import time, causing unwanted simulations
# # For the locus usage -- ts1 is a tskit.generator object !!
# ts1 = msprime.sim_ancestry(
#     samples={"O": n0, "P1": nP1, "P2": nP2, "P3": nP3},
#     demography=demography,
#     num_replicates=num_replicates,
#     ploidy=ploidy,
# )
# # For the Chromosome usage -- ts is a tskit.TreeSequence object  
# ts = msprime.sim_ancestry(
#     samples={"O": n0, "P1": nP1, "P2": nP2, "P3": nP3},
#     demography=demography,
#     sequence_length=100000000,
#     recombination_rate=0.000000001,
#     ploidy= ploidy,
# )

########################################################################################

# Add the external directory to the Python path
EXTERNAL_DIR = Path(__file__).parent / "external"
sys.path.append(str(EXTERNAL_DIR))

# Import twisst functions
try:
    from twisst import weightTrees, summary
except ImportError:
    # This will be handled by the pipeline calling download_twisst
    weightTrees = None
    summary = None


def detect_and_read_trees(
    file_path: str,
) -> Tuple[Union[tskit.TreeSequence, List[str]], str]:
    """
    Detects the format of the given tree file (ts, Newick or Nexus) and returns both:
      - if the file is a tskit.TreeSequence, it returns the tree sequence object
      - if the file is a Newick or Nexus file, it returns a list of Newick strings
      - a string indicating the format: 'ts' or 'newick'

    Supported formats:
      - TreeSequence (.trees, .ts): TSKit tree sequence files
      - Newick (.newick, .nwk, .tree): Single or multiple Newick format trees
      - Nexus (.nexus): Nexus format files

    Args:
        file_path (str): Path to the tree file.

    Returns:
        Tuple[Union[tskit.TreeSequence, List[str]], str]: (Tree data, format label)
    """
    path = Path(file_path)

    # Case 1: TreeSequence (.trees or .ts)
    if path.suffix in [".trees", ".ts"]:
        try:
            ts = tskit.load(file_path)
            print("✅ Detected format: TreeSequence (.ts/.trees)")
            return ts, "ts"
        except Exception as e:
            raise RuntimeError(f"Failed to load TreeSequence: {e}")

    # Case 2: Check for Newick files by extension first
    if path.suffix in [".newick", ".nwk", ".tree"]:
        try:
            with open(file_path, "r") as f:
                newicks = [line.strip() for line in f if line.strip()]
                print(f"✅ Detected format: Newick ({path.suffix})")
                return newicks, "newick"
        except Exception as e:
            raise ValueError(f"Failed to read Newick trees from {path.suffix} file: {e}")

    # Case 3: Try reading first line to detect Nexus or other formats
    try:
        with open(file_path, "r") as f:
            first_line = f.readline().strip()
            rest = f.read()
    except Exception as e:
        raise IOError(f"Could not read file: {e}")

    # Case 3a: Nexus format
    if first_line.upper().startswith("#NEXUS"):
        try:
            with open(file_path, "r") as f:
                lines = f.readlines()

            # Extract lines that define trees (e.g., "tree t1 = ...")
            newicks = []
            for line in lines:
                if line.strip().lower().startswith("tree"):
                    newick = line.split("=", 1)[-1].strip()
                    try:
                        t = ete3.Tree(newick)
                        newicks.append(t.write(format=0).strip())
                    except Exception as e:
                        raise ValueError(f"Invalid Newick tree in Nexus file: {e}")

            print("✅ Detected format: Nexus → Converted to Newick")
            return newicks, "newick"

        except Exception as e:
            raise ValueError(f"Failed to parse Nexus file: {e}")

    # Case 3b: Default to Newick (one tree per line) for unknown extensions
    try:
        with open(file_path, "r") as f:
            newicks = [line.strip() for line in f if line.strip()]
            print("✅ Detected format: Newick (unknown extension, treating as Newick)")
            return newicks, "newick"
    except Exception as e:
        raise ValueError(f"Failed to read Newick trees: {e}")


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


def ts_chromosome_to_twisst_weights(
    ts, outgroup=None, output_file=None, verbose=False, twisst_verbose=False
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

    # Check if twisst is available
    if weightTrees is None:
        raise ImportError(
            "twisst is not available. Please ensure twisst.py is in the "
            "twisstntern/external directory. You can download it from: "
            "https://github.com/simonhmartin/twisst"
        )

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
    # explicatly printing the topologies
    if verbose:
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

    # Create DataFrame
    n_topos = weights_norm.shape[1]

    # Create column names based on number of topologies
    if n_topos == 3:
        columns = ["T1", "T2", "T3"]  # Standard 3-topology case (4 populations)
    else:
        columns = [f"Topo{i+1}" for i in range(n_topos)]

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
        # Get simplified topologies
        simplified_topos = simplify_topologies(
            weightsData
        )  # e.g. ['(0,(1,(2,3)))', '(0,(2,(1,3)))', '(0,(3,(1,2)))']

        # Create a new DataFrame with that single row
        header_row = pd.DataFrame([simplified_topos], columns=columns)

        # Concatenate the header and the weights
        full_df = pd.concat([header_row, df], ignore_index=True)

        # Save
        full_df.to_csv(output_file, index=False)

    return df


def ts_to_twisst_weights(
    input_data, outgroup=None, output_file=None, verbose=False, twisst_verbose=False
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
        )
    #########################################################
    # Generator case - process multiple TreeSequences
    input_data = list(input_data)
    if verbose:
        print("Processing generator of TreeSequences...")

    all_weights = []
    canonical_topologies = None
    canonical_simplified_topos = None
    columns = None
    total_processed = 0

    for i in range(len(input_data)):
        ts = input_data[i]
        if verbose:
            print(f"\nProcessing TreeSequence {i+1}...")

        # Debug info about the TreeSequence
        if verbose:
            print(f"  Number of populations: {ts.num_populations}")
            print(f"  Number of samples: {ts.num_samples}")
            print(f"  Number of trees: {ts.num_trees}")

        # Get populations that actually have samples
        populations_with_samples = []
        for pop_id in range(ts.num_populations):
            samples = [s for s in ts.samples() if ts.node(s).population == pop_id]
            if len(samples) > 0:
                populations_with_samples.append(str(pop_id))

        if len(populations_with_samples) < 3:
            if verbose:
                print(
                    f"  Skipping: only {len(populations_with_samples)} populations with samples"
                )
            continue

        # Set default outgroup to first population if not specified
        current_outgroup = outgroup
        if current_outgroup is None:
            current_outgroup = populations_with_samples[0]

        # Validate outgroup
        if current_outgroup not in populations_with_samples:
            if verbose:
                print(
                    f"  Skipping: outgroup {current_outgroup} not in populations with samples"
                )
            continue

        if verbose:
            print(
                f"  Analyzing {len(populations_with_samples)} populations: {populations_with_samples}"
            )
            print(f"  Using outgroup: {current_outgroup}")

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

            if verbose:
                print(f"  Established canonical topology order:")
                for j, topo in enumerate(canonical_simplified_topos):
                    print(f"    {columns[j]}: {topo}")

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

            if verbose:
                print(f"  Extracted {len(df_current)} valid trees")
        else:
            if verbose:
                print(f"  No valid trees found")

    # Combine all weights into single DataFrame
    if not all_weights:
        raise ValueError("No valid topology weights found across all TreeSequences")

    combined_df = pd.concat(all_weights, ignore_index=True)

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
        full_df.to_csv(output_file, index=False, float_format='%.3f')

        if verbose:
            print(f"  Saved results to: {output_file}")

    return combined_df


def newick_to_twisst_weights(
    newick_trees,
    taxon_names=None,
    outgroup=None,
    output_file=None,
    verbose=False,
    twisst_verbose=False,
):
    """
    Extract topology weights from Newick tree strings using twisst.

    Args:
        newick_trees (List[str] or str): List of Newick tree strings or single Newick string
        taxon_names (List[str], optional): List of population names. If None, will be inferred from sample names.
        outgroup (str, optional): Population name to use as outgroup. If None, uses the first population.
        output_file (str, optional): Path to save CSV file. If None, returns DataFrame.
        verbose (bool): Whether to print verbose output
        twisst_verbose (bool): Whether to print verbose twisst output

    Returns:
        pd.DataFrame: Normalized weights ready for analysis
    """

    # Check if twisst is available
    if weightTrees is None:
        raise ImportError(
            "twisst is not available. Please ensure twisst.py is in the "
            "twisstntern/external directory. You can download it from: "
            "https://github.com/simonhmartin/twisst"
        )

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
                print("Using detected populations...")
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
        taxonNames=population_names,  # Dynamic population names
        outgroup=outgroup,
        verbose=twisst_verbose,
    )

    # Print topologies if verbose
    if verbose:
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

    # Create DataFrame
    n_topos = weights_norm.shape[1]

    # Create column names based on number of topologies
    if n_topos == 3:
        columns = ["T1", "T2", "T3"]  # Standard 3-topology case (4 populations)
    else:
        columns = [f"Topo{i+1}" for i in range(n_topos)]

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
        # Get simplified topologies
        simplified_topos = simplify_topologies(weightsData)

        # Create a new DataFrame with that single row
        header_row = pd.DataFrame([simplified_topos], columns=columns)

        # Concatenate the header and the weights
        full_df = pd.concat([header_row, df], ignore_index=True)

        # Save
        full_df.to_csv(output_file, index=False, float_format='%.3f')

        if verbose:
            print(f"Saved results to: {output_file}")

    return df


def trees_to_twisst_weights_unified(
    file_path,
    taxon_names=None,
    outgroup=None,
    output_file=None,
    verbose=False,
    twisst_verbose=False,
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

    #if verbose:
       # print(f"Detected format: {format_type}")

    # Route to appropriate processing function based on detected format
    if format_type == "ts":
        return ts_to_twisst_weights(
            tree_data,
            outgroup=outgroup,
            output_file=str(output_file),
            verbose=verbose,
            twisst_verbose=twisst_verbose,
        )
    elif format_type == "newick":
        return newick_to_twisst_weights(
            tree_data,
            taxon_names=taxon_names,
            outgroup=outgroup,
            output_file=str(output_file),
            verbose=verbose,
            twisst_verbose=twisst_verbose,
        )
    else:
        raise ValueError(f"Unsupported tree format: {format_type}")


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


# A thin wrapper
def process_trees_from_file(file_path, **kwargs):
    """
    Convenience wrapper for trees_to_twisst_weights_unified.
    Maintains backward compatibility with existing code.
    """
    return trees_to_twisst_weights_unified(file_path, **kwargs)


# def infer_ploidy_from_samples(newick_tree):
#     """
#     Infer ploidy from sample names in a Newick tree.
    
#     Args:
#         newick_tree (str): Newick tree string
        
#     Returns:
#         int: Inferred ploidy (1 for haploid, 2 for diploid)
#     """
#     try:
#         tree = ete3.Tree(newick_tree)
#         all_sample_names = [leaf.name for leaf in tree.get_leaves()]
        
#         # Group samples by population
#         population_samples = {}
#         for sample in all_sample_names:
#             if "_" in sample:
#                 pop_name = sample.split("_")[0]
#                 if pop_name not in population_samples:
#                     population_samples[pop_name] = []
#                 population_samples[pop_name].append(sample)
        
#         # Check sample counts per population
#         sample_counts = [len(samples) for samples in population_samples.values()]
        
#         # If all populations have even number of samples, likely diploid
#         # If all populations have odd number of samples, likely haploid
#         # If mixed, default to haploid
#         if all(count % 2 == 0 for count in sample_counts):
#             return 2
#         else:
#             return 1
            
#     except Exception as e:
#         print(f"Warning: Could not infer ploidy: {e}")
#         return 1  # Default to haploid if inference fails
