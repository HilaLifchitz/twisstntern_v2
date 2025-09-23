"""
Module for simulating demographic scenarios and generating tree sequences.
This module provides functionality to:
1. Simulate demographic scenarios using msprime
2. Generate tree sequences for different simulation modes:
   - Locus mode: independent non-recombining loci
   - Chromosome mode: recombining chromosome

The simulation module supports two main modes:

1. Locus Mode:
   - Simulates independent non-recombining loci
   - Each locus is simulated independently
   - Useful for studying population structure without recombination
   - Parameters:
     * n_loci: Number of loci to simulate
     * locus_length: Length of each locus in base pairs - can be set to 1 for a single locus

2. Chromosome Mode:
   - Simulates a recombining chromosome
   - Models recombination along the chromosome
   - Useful for studying linkage and recombination
   - Parameters:
     * chromosome_length: Total length of chromosome
     * rec_rate: Recombination rate per base per generation

Example usage:
    from twisstntern_simulate.config import Config
    from twisstntern_simulate.simulation import run_simulation

    # Load configuration from YAML file
    config = Config("config.yaml")

    # Run simulation (trees are automatically saved)
    results = run_simulation(config, output_dir="results")

    # Access results
    if config.simulation_mode == 'locus':
        ts_locus = results['locus']
    elif config.simulation_mode == 'chromosome':
        ts_chrom = results['chromosome']
"""

import msprime
import tskit
from typing import Optional
from pathlib import Path
import logging
from .config import Config
import random

# Get logger (logging configured in __init__.py)
logger = logging.getLogger(__name__)
######################################################################################
# SIMULATION FUNCTIONS
######################################################################################


def simulate_locus(config: Config):
    """
    Simulates independent non-recombining loci.

    This function simulates multiple independent loci without recombination. Each locus
    is simulated using the demographic model specified in the configuration.

    Args:
        config: Configuration object containing demographic parameters

    Returns:
        Generator of tskit.TreeSequence: Generator yielding TreeSequence objects for each locus

    Note:
        The recombination rate is set to 0 for locus mode simulations.
    """
    # Creating demographic model from the yaml parameters
    demography = msprime.Demography()

    # Add populations
    for pop in config.populations:
        demography.add_population(
            name=pop.name, initial_size=pop.Ne, growth_rate=pop.growth_rate
        )

    # Add population splits
    for split in config.splits:
        demography.add_population_split(
            time=split.time,
            derived=[split.derived_pop1, split.derived_pop2],
            ancestral=split.ancestral_pop,
        )

    # Add migration rates
    if hasattr(config, "migration") and config.migration:
        for migration_route, rate in config.migration.items():
            if rate > 0:  # Only add non-zero migration rates
                # Parse migration route like "p1>p2" -> source="p1", dest="p2"
                source, dest = migration_route.split(">")
                demography.set_migration_rate(source=source, dest=dest, rate=rate)

    # Default values, in case the user hasn't specified them
    if config.locus_length:
        locus_length = config.locus_length
    else:
        locus_length = 1  # we don't need a locus length, all we care about are the trees themselves -> it is 1 by default

    # Default is haploid
    if config.ploidy:
        ploidy = config.ploidy
    else:
        ploidy = 1

    # Default if no random seed is specified, we specify a random seed and log it
    if config.seed:
        seed = config.seed
    else:
        seed = random.randint(0, 2**32 - 1)
        logger.info(f"Using random seed: {seed}")
        print(f"Using random seed: {seed}")  # for the user to see

    # Simulate tree sequence
    ts = msprime.sim_ancestry(
        samples={
            pop.name: pop.sample_size for pop in config.populations if pop.sample_size
        },
        demography=demography,
        num_replicates=config.n_loci,
        sequence_length=locus_length,
        ploidy=ploidy,
        random_seed=seed,
        recombination_rate=0,  # No recombination for locus mode
    )

    return ts


def simulate_chromosome(config: Config) -> tskit.TreeSequence:
    """
    Simulates a chromosome with recombination.

    This function simulates a chromosome with recombination using the
    demographic model specified in the configuration. The recombination
    rate is applied along the entire chromosome length.

    Args:
        config: Configuration object containing demographic parameters

    Returns:
        tskit.TreeSequence: Tree sequence for the simulated chromosome

    Note:
        The recombination rate is applied per base pair per generation.
    """
    # Create demographic model
    demography = msprime.Demography()

    # Add populations
    for pop in config.populations:
        demography.add_population(
            name=pop.name, initial_size=pop.Ne, growth_rate=pop.growth_rate
        )

    # Add population splits
    for split in config.splits:
        demography.add_population_split(
            time=split.time,
            derived=[split.derived_pop1, split.derived_pop2],
            ancestral=split.ancestral_pop,
        )

    # Add migration rates
    if hasattr(config, "migration") and config.migration:
        for migration_route, rate in config.migration.items():
            if rate > 0:  # Only add non-zero migration rates
                # Parse migration route like "p1>p2" -> source="p1", dest="p2"
                source, dest = migration_route.split(">")
                demography.set_migration_rate(source=source, dest=dest, rate=rate)

    # Default is haploid
    if config.ploidy:
        ploidy = config.ploidy
    else:
        ploidy = 1

    # Default if no random seed is specified, we specify a random seed and log it
    if config.seed:
        seed = config.seed
    else:
        seed = random.randint(0, 2**32 - 1)
        logger.info(f"Using random seed: {seed}")
        print(f"Using random seed: {seed}")  # for the user to see

    # Simulate tree sequence
    ts = msprime.sim_ancestry(
        samples={
            pop.name: pop.sample_size for pop in config.populations if pop.sample_size
        },
        demography=demography,
        sequence_length=config.chromosome_length,
        recombination_rate=config.rec_rate,
        ploidy=ploidy,
        random_seed=seed,
    )
    return ts


######################################################################################
# MAIN SIMULATION FUNCTION
######################################################################################


def run_simulation(
    config: Config, output_dir: str, mode_override: Optional[str] = None
) -> dict:
    """
    Runs simulation based on the specified mode in config.
    This function acts as a dispatcher for different simulation modes,
    running the appropriate simulation mode based on the configuration.
    Trees are automatically saved in Newick format.

    Args:
        config: Configuration object containing simulation parameters (with any overrides already applied)
        output_dir: Directory to save tree files
        mode_override: Optional override for simulation mode (overrides config file)

    Returns:
        dict: Dictionary containing simulation results:
            - 'locus': Tree sequence for locus mode (if requested)
            - 'chromosome': Tree sequence for chromosome mode (if requested)
            - 'newick_file': Path to saved Newick file

    Note:
        Only one mode (locus OR chromosome) will be run based on config.simulation_mode.
    """
    results = {}

    # Apply mode override if provided
    if mode_override is not None:
        # Type assertion for Pylance - we know mode_override is a string here
        assert isinstance(mode_override, str)
        config.simulation_mode = mode_override

    # Run locus simulation if requested
    if config.simulation_mode == "locus":
        print("Simulating independent non-recombining loci...")
        ts_locus = simulate_locus(config)
        print(f"Generated {config.n_loci} non-recombining loci")

        # Convert generator to list once for both saving and processing
        ts_list = list(ts_locus)

        # Store the list (not the exhausted generator) for pipeline processing
        results["locus"] = ts_list

        # Ensure output directory exists
        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # Save trees as Newick format
        newick_path = Path(output_dir) / f"{config.simulation_mode}_trees.newick"
        newick_file = save_ts_LOCUS_as_plain_newick(ts_list, newick_path)

        print(f"✅ Saved trees: {newick_file}")

        results["newick_file"] = newick_file

    # Run chromosome simulation if requested
    elif config.simulation_mode == "chromosome":
        print("Simulating recombining chromosome...")
        ts_chrom = simulate_chromosome(config)
        results["chromosome"] = ts_chrom
        print(
            f"Generated chromosome of length {config.chromosome_length:.1e} with reocmbination rate of {config.rec_rate:.1e}"
        )

        # Always save trees
        # Ensure output directory exists
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        newick_path = Path(output_dir) / f"{config.simulation_mode}_trees.newick"
        newick_file = save_ts_CHROM_as_newick(ts_chrom, newick_path)
        results["newick_file"] = newick_file

    else:
        raise ValueError(
            f"Unknown simulation mode: {config.simulation_mode}. Must be 'locus' or 'chromosome'"
        )

    return results


######################################################################################
# HELPER FUNCTIONS FOR SAVING TREE SEQUENCES AS NEWICK FILES
######################################################################################


def get_population_map(ts):
    """
    Create a mapping from sample ID to unique population sample name for a TreeSequence.

    Args:
        ts: msprime TreeSequence object

    Returns:
        dict: Mapping from sample_id (int) to population_sample_name (str) like "P1_1", "P1_2", etc.
    """
    pop_map = {}

    # Create mapping from population ID to population name
    pop_id_to_name = {}
    for pop_id in range(ts.num_populations):
        pop = ts.population(pop_id)
        if pop.metadata and "name" in pop.metadata:
            pop_name = pop.metadata["name"]
        else:
            pop_name = str(pop_id)  # fallback to ID if no name
        pop_id_to_name[pop_id] = pop_name

    # Count samples within each population to create unique identifiers
    pop_sample_counts = {}

    # Map each sample to its unique population sample name
    for sample_id in ts.samples():
        node = ts.node(sample_id)
        pop_name = pop_id_to_name[node.population]

        # Increment counter for this population
        if pop_name not in pop_sample_counts:
            pop_sample_counts[pop_name] = 0
        pop_sample_counts[pop_name] += 1

        # Create unique sample identifier like "P1_1", "P1_2", etc.
        unique_sample_name = f"{pop_name}_{pop_sample_counts[pop_name]}"
        pop_map[sample_id] = unique_sample_name

    return pop_map


def create_newick_with_sample_labels(tree, pop_map):
    """
    Create a Newick string from a tskit tree with proper sample labels.

    Args:
        tree: tskit Tree object
        pop_map: dict mapping sample_id to population_name

    Returns:
        str: Newick string with population labels for samples
    """

    def _get_newick_recursive(node):
        """Recursively build Newick string"""
        if tree.is_sample(node):
            # This is a sample - use population label
            pop_label = pop_map.get(node, f"Sample_{node}")
            return pop_label
        else:
            # This is an internal node - recurse on children
            children = tree.children(node)
            if len(children) == 0:
                return ""

            child_strings = []
            for child in children:
                child_str = _get_newick_recursive(child)
                if child_str:  # Only add if not empty
                    # Add branch length if available
                    branch_length = tree.branch_length(child)
                    if branch_length > 0:
                        child_str += f":{branch_length}"
                    child_strings.append(child_str)

            if len(child_strings) == 1:
                return child_strings[0]
            else:
                return f"({','.join(child_strings)})"

    # Start from root
    root = tree.root
    newick = _get_newick_recursive(root)

    # Add semicolon at the end
    if not newick.endswith(";"):
        newick += ";"

    return newick


#######################################################################################
# FOR THE CHROMOSOME MODE:
def save_ts_CHROM_as_newick(ts, output_path):
    """
    Saves marginal trees as Newick format with population labels, no interval annotations.
    This format is ideal for twisst.
    """
    pop_map = get_population_map(ts)

    with open(output_path, "w") as f:
        for tree in ts.trees():
            labeled_newick = create_newick_with_sample_labels(tree, pop_map)
            f.write(labeled_newick + "\n")

    return str(output_path)


# EXAMPLE:
# save_ts_chromosome_as_newick(ts, os.path.join(output_dir, "CHROM_pop_plain.newick"))
##################################################################################
# LOCUS MODE:
##################################################################################


def save_ts_LOCUS_as_plain_newick(ts_list, output_path):
    """
    Save TreeSequence genealogies as Newick format using plain format.
    This creates files that work properly with twisst.
    """
    with open(output_path, "w") as file:
        for replicate_index, ts in enumerate(ts_list):
            pop_map = get_population_map(ts)
            for tree in ts.trees():
                labeled_newick = create_newick_with_sample_labels(tree, pop_map)
                file.write(labeled_newick + "\n")

    return str(output_path)
