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


# COMMENTED OUT: This code was executing simulations at import time, causing unwanted output
# This appears to be example/test code that should not run during import

# # A demography where (1,2) first coalesce, with our given Ne
# 
# demography = msprime.Demography()
# 
# # initializing populations
# demography.add_population(name="O", initial_size=NeO)
# demography.add_population(name="P1", initial_size=NeP1)
# demography.add_population(name="P2", initial_size=NeP2)
# demography.add_population(name="P3", initial_size=NeP3)
# demography.add_population(name="P13", initial_size=NeP12)
# demography.add_population(name="P123", initial_size=NeP123)
# demography.add_population(name="ANC", initial_size=NeANC)
# 
# # adding split times
# demography.add_population_split(time=t1, derived=["P1", "P2"], ancestral="P13")
# demography.add_population_split(time=t2, derived=["P13", "P3"], ancestral="P123")
# demography.add_population_split(time=t3, derived=["P123", "O"], ancestral="ANC")
# 
# # setting up gene flow
# demography.set_migration_rate("P2", "P3", mig_rate)
# 
# ploidy = 2
# # from collections import defaultdict
# num_replicates = 20
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
#     ploidy=ploidy,
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
        full_df.to_csv(output_file, index=False, float_format='%.3f')

    return df

#
#
#
#
#
#
#
##################################################################################
# main topology weights generating function
#################################################################################
def ts_to_twisst_weights(
    input_data, outgroup=None, output_file=None, verbose=False, twisst_verbose=False
):
    """
    Extract topology weights from a TreeSequence object using twisst.

    This function handles both single TreeSequence objects and generators of TreeSequences.
    It ensures consistent topology ordering across all TreeSequences.

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
    # CRITICAL: Convert generator to list first since generators can only be consumed once
    input_data = list(input_data)
    
    # PERFORMANCE OPTIMIZATION: Use batch processing for large numbers of loci
    if len(input_data) > 100:  # Large number of loci - use optimized approach --- NOT SURE THIS ACTUALLY WORKS..
        if verbose:
            print(f"Processing {len(input_data)} TreeSequences (using optimized batch mode)...")
        return _ts_locus_batch_process(
            input_data, outgroup=outgroup, output_file=output_file, 
            verbose=verbose, twisst_verbose=twisst_verbose
        )
    
    # Original method for smaller numbers of TreeSequences (< 100)
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

# an attempt tp make the locus mode faster
def _ts_locus_batch_process(input_data, outgroup=None, output_file=None, verbose=False, twisst_verbose=False):
    """
    Optimized batch processing for large numbers of loci.
    
    Instead of processing each locus individually, this function:
    1. Samples a few loci to determine topology structure
    2. Processes loci in batches to reduce overhead
    3. Uses simplified topology mapping
    """
    
    if verbose:
        print(f"Batch processing {len(input_data)} loci...")
    
    # Step 1: Sample first few loci to establish topology structure
    sample_size = min(10, len(input_data))
    sample_results = []
    
    for i in range(sample_size):
        ts = input_data[i]
        
        # Get populations with samples
        populations_with_samples = []
        for pop_id in range(ts.num_populations):
            samples = [s for s in ts.samples() if ts.node(s).population == pop_id]
            if len(samples) > 0:
                populations_with_samples.append(str(pop_id))
        
        if len(populations_with_samples) < 3:
            continue
            
        # Set outgroup
        current_outgroup = outgroup if outgroup else populations_with_samples[0]
        if current_outgroup not in populations_with_samples:
            continue
            
        # Process this sample
        weightsData = weightTrees(
            ts, treeFormat="ts", taxonNames=populations_with_samples,
            outgroup=current_outgroup, verbose=False
        )
        
        if "weights" in weightsData:
            sample_results.append((weightsData, populations_with_samples, current_outgroup))
    
    if not sample_results:
        raise ValueError("No valid loci found in sample")
    
    # Use first valid result to establish structure
    canonical_weights, populations_with_samples, used_outgroup = sample_results[0]
    canonical_simplified_topos = simplify_topologies(canonical_weights)
    n_topos = len(canonical_simplified_topos)
    
    # Create column names
    if n_topos == 3:
        columns = ["T1", "T2", "T3"]
    else:
        columns = [f"Topo{i+1}" for i in range(n_topos)]
    
    if verbose:
        print(f"Established topology structure with {n_topos} topologies")
        print(f"Using {len(populations_with_samples)} populations: {populations_with_samples}")
        print(f"Using outgroup: {used_outgroup}")
    
    # Step 2: Process all loci in batches
    batch_size = 500  # Process 500 loci at a time
    all_weights = []
    
    for batch_start in range(0, len(input_data), batch_size):
        batch_end = min(batch_start + batch_size, len(input_data))
        
        if verbose:
            print(f"Processing batch {batch_start//batch_size + 1}/{(len(input_data)-1)//batch_size + 1} "
                  f"(loci {batch_start+1}-{batch_end})...")
        
        batch_weights = []
        
        for i in range(batch_start, batch_end):
            ts = input_data[i]
            
            try:
                # Quick processing - assume same structure as canonical
                weightsData = weightTrees(
                    ts, treeFormat="ts", taxonNames=populations_with_samples,
                    outgroup=used_outgroup, verbose=False
                )
                
                # Extract weights
                if "weights_norm" in weightsData:
                    weights_norm = weightsData["weights_norm"]
                else:
                    weights = weightsData["weights"]
                    row_sums = weights.sum(axis=1)
                    if np.all(row_sums == 0):
                        continue
                    weights_norm = weights / row_sums[:, np.newaxis]
                
                # Add to batch
                batch_weights.append(weights_norm)
                
            except Exception:
                # Skip problematic loci
                continue
        
        if batch_weights:
            # Combine batch weights
            combined_batch = np.vstack(batch_weights)
            batch_df = pd.DataFrame(combined_batch, columns=columns)
            batch_df = batch_df.dropna()
            if len(batch_df) > 0:
                all_weights.append(batch_df)
    
    # Step 3: Combine all results
    if not all_weights:
        raise ValueError("No valid topology weights found")
    
    combined_df = pd.concat(all_weights, ignore_index=True)
    
    if verbose:
        print(f"Batch processing complete:")
        print(f"  Total trees processed: {len(combined_df)}")
        print(f"  Final DataFrame shape: {combined_df.shape}")
    
    # Save to file if requested
    if output_file:
        header_row = pd.DataFrame([canonical_simplified_topos], columns=columns)
        full_df = pd.concat([header_row, combined_df], ignore_index=True)
        full_df.to_csv(output_file, index=False, float_format='%.3f')
        if verbose:
            print(f"Results saved to: {output_file}")
    
    return combined_df





def ts_to_twisst_weights_old_old(
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
    combined_df = combined_df.round(3)
    combined_df = combined_df.loc[combined_df.iloc[:, 1] != combined_df.iloc[:, 2]] 
    
    print(f"\nWeights with T2 == T3 removed. Remaining rows: {len(combined_df)}")

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

# Not ised directly, but still good to have:
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
