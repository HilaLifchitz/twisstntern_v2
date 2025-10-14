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
    from omegaconf import DictConfig
    from twisstntern_simulate.simulation import run_simulation

    # Load configuration from YAML file
    # config is now passed as a Hydra DictConfig object

    # Run simulation (trees are automatically saved)
    results = run_simulation(config, output_dir="results")

    # Access results
    if config.simulation.mode == 'locus':
        ts_locus = results['locus']
    elif config.simulation.mode == 'chromosome':
        ts_chrom = results['chromosome']
"""

import os
import logging
import random
from pathlib import Path
from typing import Dict, List, Union, Tuple, Optional, Literal

import msprime
import tskit
from omegaconf import DictConfig 

from ..core.hydra_utils import build_run_suffix
from ..core.tree_export import (
    save_ts_CHROM_as_newick,
    save_ts_LOCUS_as_plain_newick,
)

# Get logger (logging configured in __init__.py)
logger = logging.getLogger(__name__)
######################################################################################
# SIMULATION FUNCTIONS
######################################################################################

def simulate_locus(config: DictConfig):
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
    for pop in config.simulation.populations:
        demography.add_population(
            name=pop.name, initial_size=pop.Ne, growth_rate=pop.growth_rate
        )

    # Add population splits
    for split in config.simulation.splits:
        demography.add_population_split(
            time=split.time,
            derived=[split.derived_pop1, split.derived_pop2],
            ancestral=split.ancestral_pop,
        )

    # Add migration rates
    if hasattr(config.simulation, 'migration') and config.simulation.migration:
        for migration_route, rate in config.simulation.migration.items():
            if rate > 0:  # Only add non-zero migration rates
                # Parse migration route like "p1>p2" -> source="p1", dest="p2"
                source, dest = migration_route.split('>')
                demography.set_migration_rate(source=source, dest=dest, rate=rate)

    # Default values, in case the user hasn't specified them
    if config.simulation.locus_length:
        locus_length = config.simulation.locus_length
    else:
        locus_length = 1  # we don't need a locus length, all we care about are the trees themselves -> it is 1 by default

    # Default is haploid
    if config.simulation.ploidy:
        ploidy = config.simulation.ploidy
    else:
        ploidy = 1

    # Default if no random seed is specified, we specify a random seed and log it
    if config.seed:
        seed = config.seed
    else:
        seed = random.randint(0, 2**32 - 1)
        logger.info(f"Using random seed: {seed}")
        print(f"Using random seed: {seed}")  # for the user to see
        config.seed = seed

    # Simulate tree sequence
    samples = {}
    for pop in config.simulation.populations:
        sample_size = getattr(pop, 'sample_size', 0)
        if sample_size is None:
            sample_size = 0
        if sample_size > 0:
            samples[pop.name] = sample_size

    ts = msprime.sim_ancestry(
        samples=samples,
        demography=demography,
        num_replicates=config.simulation.n_loci,
        sequence_length=locus_length,
        ploidy=ploidy,
        random_seed=seed,
        recombination_rate=0,  # No recombination for locus mode
    )

    return ts


def simulate_chromosome(config: DictConfig) -> tskit.TreeSequence:
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
    for pop in config.simulation.populations:
        demography.add_population(
            name=pop.name, initial_size=pop.Ne, growth_rate=pop.growth_rate
        )

    # Add population splits
    for split in config.simulation.splits:
        demography.add_population_split(
            time=split.time,
            derived=[split.derived_pop1, split.derived_pop2],
            ancestral=split.ancestral_pop,
        )

    # Add migration rates
    if hasattr(config.simulation, 'migration') and config.simulation.migration:
        for migration_route, rate in config.simulation.migration.items():
            if rate > 0:  # Only add non-zero migration rates
                # Parse migration route like "p1>p2" -> source="p1", dest="p2"
                source, dest = migration_route.split('>')
                demography.set_migration_rate(source=source, dest=dest, rate=rate)

    # Default is haploid
    if config.simulation.ploidy:
        ploidy = config.simulation.ploidy
    else:
        ploidy = 1

    # Default if no random seed is specified, we specify a random seed and log it
    if config.seed:
        seed = config.seed
    else:
        seed = random.randint(0, 2**32 - 1)
        logger.info(f"Using random seed: {seed}")
        print(f"Using random seed: {seed}")  # for the user to see
        config.seed = seed

    # Simulate tree sequence
    samples = {}
    for pop in config.simulation.populations:
        sample_size = getattr(pop, 'sample_size', 0)
        if sample_size is None:
            sample_size = 0
        if sample_size > 0:
            samples[pop.name] = sample_size

    ts = msprime.sim_ancestry(
        samples=samples,
        demography=demography,
        sequence_length=config.simulation.chromosome_length,
        recombination_rate=config.simulation.rec_rate,
        ploidy=ploidy,
        random_seed=seed,
    )
    return ts

######################################################################################
# MAIN SIMULATION FUNCTION
######################################################################################


def run_simulation(config: DictConfig, output_dir: str, mode_override: Optional[str] = None) -> dict:
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
        Only one mode (locus OR chromosome) will be run based on config.simulation.mode.
    """
    results = {}
    
    # Apply mode override if provided
    if mode_override is not None:
        # Type assertion for Pylance - we know mode_override is a string here
        assert isinstance(mode_override, str)
        config.simulation.mode = mode_override
    
    # Run locus simulation if requested
    if config.simulation.mode == "locus":
        print("Simulating independent non-recombining loci...")
        ts_locus = simulate_locus(config)
        print(f"Generated {config.simulation.n_loci} non-recombining loci")

        # Convert generator to list once for both saving and processing
        ts_list = list(ts_locus)
        
        # Store the list (not the exhausted generator) for pipeline processing
        results["locus"] = ts_list
        
        # Ensure output directory exists
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # Save trees as Newick format with sweep-aware naming
        seed_value = getattr(config, "seed", None)
        extra_parts = [f"seed={seed_value}"] if seed_value is not None else None
        suffix = build_run_suffix(config, extra_parts=extra_parts)
        base_name = f"{config.simulation.mode}_trees"
        if suffix != "run":
            base_name = f"{base_name}_{suffix}"
        newick_path = Path(output_dir) / f"{base_name}.newick"
        newick_file = save_ts_LOCUS_as_plain_newick(ts_list, newick_path)
        
        print(f"✅ Saved trees: {newick_file}")
        
        results["newick_file"] = newick_file

    # Run chromosome simulation if requested
    elif config.simulation.mode == "chromosome":
        print("Simulating recombining chromosome...")
        ts_chrom = simulate_chromosome(config)
        results["chromosome"] = ts_chrom
        print(f"Generated chromosome of length {config.simulation.chromosome_length:.1e} with reocmbination rate of {config.simulation.rec_rate:.1e}")
        
        # Always save trees
        # Ensure output directory exists
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        seed_value = getattr(config, "seed", None)
        extra_parts = [f"seed={seed_value}"] if seed_value is not None else None
        suffix = build_run_suffix(config, extra_parts=extra_parts)
        base_name = f"{config.simulation.mode}_trees"
        if suffix != "run":
            base_name = f"{base_name}_{suffix}"
        newick_path = Path(output_dir) / f"{base_name}.newick"
        newick_file = save_ts_CHROM_as_newick(ts_chrom, newick_path)
        results["newick_file"] = newick_file
    
    else:
        raise ValueError(f"Unknown simulation mode: {config.simulation.mode}. Must be 'locus' or 'chromosome'")

    return results
