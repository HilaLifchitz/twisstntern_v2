"""
Module for simulating demographic scenarios and generating tree sequences.
This module provides functionality to:
1. Simulate demographic scenarios using msprime
2. Generate tree sequences for different simulation modes:
   - Locus mode: independent non-recombining loci
   - Chromosome mode: recombining chromosome
   - Both modes: generates both types of simulations

The simulation module supports two main modes:

1. Locus Mode:
   - Simulates independent non-recombining loci
   - Each locus is simulated independently
   - Useful for studying population structure without recombination
   - Parameters:
     * n_loci: Number of loci to simulate
     * locus_length: Length of each locus in base pairs

2. Chromosome Mode:
   - Simulates a recombining chromosome
   - Models recombination along the chromosome
   - Useful for studying linkage and recombination
   - Parameters:
     * chromosome_length: Total length of chromosome
     * rec_rate: Recombination rate per base per generation

3. Both Modes:
   - Runs both locus and chromosome simulations
   - Useful for comparing results between modes
   - Requires parameters for both modes

Example usage:
    from twisstntern.config import Config
    from twisstntern.simulation import run_simulation

    # Load configuration from YAML file
    config = Config("config.yaml")

    # Run simulation
    results = run_simulation(config)

    # Access results
    if 'locus' in results:
        ts_locus = results['locus']
    if 'chromosome' in results:
        ts_chrom = results['chromosome']
"""

import os
import msprime
import numpy as np
import pandas as pd
from typing import Dict, List, Union, Tuple
from pathlib import Path
import json
import logging
from .config import Config
from .core import dump_data

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class DemographicModel:
    """Class to handle demographic model configuration and simulation."""
    
    def __init__(self, config_file: str):
        """
        Initialize demographic model from configuration file.
        
        Args:
            config_file (str): Path to configuration file
        """
        self.config = self._load_config(config_file)
        self._validate_config()
        
    def _load_config(self, config_file: str) -> Dict:
        """Load and parse configuration file."""
        try:
            with open(config_file, 'r') as f:
                config = {}
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        if '=' in line:
                            key, value = line.split('=')
                            key = key.strip()
                            value = value.strip()
                            # Convert numeric values
                            try:
                                value = float(value)
                            except ValueError:
                                pass
                            config[key] = value
                return config
        except Exception as e:
            raise Exception(f"Error loading configuration file: {str(e)}")
    
    def _validate_config(self):
        """Validate configuration parameters."""
        required_params = [
            't1', 't2', 't3',
            'ne_p1', 'ne_p2', 'ne_p3', 'ne_p12', 'ne_p123', 'ne_O'
        ]
        for param in required_params:
            if param not in self.config:
                raise ValueError(f"Missing required parameter: {param}")
    
    def create_demography(self) -> msprime.Demography:
        """Create msprime demography object from configuration."""
        demography = msprime.Demography()
        
        # Add populations
        demography.add_population(name="p1", initial_size=self.config['ne_p1'])
        demography.add_population(name="p2", initial_size=self.config['ne_p2'])
        demography.add_population(name="p3", initial_size=self.config['ne_p3'])
        demography.add_population(name="p12", initial_size=self.config['ne_p12'])
        demography.add_population(name="p123", initial_size=self.config['ne_p123'])
        demography.add_population(name="O", initial_size=self.config['ne_O'])
        
        # Add population splits
        demography.add_population_split(
            time=self.config['t1'],
            derived=["p1", "p2"],
            ancestral="p12"
        )
        demography.add_population_split(
            time=self.config['t2'],
            derived=["p12", "p3"],
            ancestral="p123"
        )
        demography.add_population_split(
            time=self.config['t3'],
            derived=["p123", "O"],
            ancestral="O"
        )
        
        # Add migration events
        migration_rates = {k: v for k, v in self.config.items() if k.startswith('m_')}
        for rate_name, rate in migration_rates.items():
            if rate > 0:
                source, dest = rate_name[2:].split('>')
                demography.add_migration_rate_change(
                    time=0,
                    rate=rate,
                    source=source,
                    dest=dest
                )
        
        return demography

def simulate_locus(config: Config, locus_length: int = 10000) -> msprime.TreeSequence:
    """
    Simulates independent non-recombining loci.
    
    This function simulates a single locus without recombination. The locus
    is simulated using the demographic model specified in the configuration.
    
    Args:
        config: Configuration object containing demographic parameters
        locus_length: Length of each locus in base pairs (default: 10000)
    
    Returns:
        msprime.TreeSequence: Tree sequence for the simulated locus
    
    Note:
        The recombination rate is set to 0 for locus mode simulations.
    """
    # Create demographic model
    demography = msprime.Demography()
    
    # Add populations
    for pop in config.populations:
        demography.add_population(
            name=pop.name,
            initial_size=pop.Ne,
            growth_rate=pop.growth_rate
        )
    
    # Add population splits
    for split in config.splits:
        demography.add_population_split(
            time=split.time,
            derived=[split.derived_pop],
            ancestral=split.ancestral_pop
        )
    
    # Simulate tree sequence
    ts = msprime.sim_ancestry(
        samples={pop.name: pop.sample_size for pop in config.populations},
        demography=demography,
        sequence_length=locus_length,
        recombination_rate=0,  # No recombination for locus mode
        random_seed=config.seed
    )
    
    return ts

def simulate_chromosome(config: Config) -> msprime.TreeSequence:
    """
    Simulates a chromosome with recombination.
    
    This function simulates a chromosome with recombination using the
    demographic model specified in the configuration. The recombination
    rate is applied along the entire chromosome length.
    
    Args:
        config: Configuration object containing demographic parameters
    
    Returns:
        msprime.TreeSequence: Tree sequence for the simulated chromosome
    
    Note:
        The recombination rate is applied per base pair per generation.
    """
    # Create demographic model
    demography = msprime.Demography()
    
    # Add populations
    for pop in config.populations:
        demography.add_population(
            name=pop.name,
            initial_size=pop.Ne,
            growth_rate=pop.growth_rate
        )
    
    # Add population splits
    for split in config.splits:
        demography.add_population_split(
            time=split.time,
            derived=[split.derived_pop],
            ancestral=split.ancestral_pop
        )
    
    # Simulate tree sequence
    ts = msprime.sim_ancestry(
        samples={pop.name: pop.sample_size for pop in config.populations},
        demography=demography,
        sequence_length=config.chromosome_length,
        recombination_rate=config.rec_rate,
        random_seed=config.seed
    )
    
    return ts

def run_simulation(config: Config) -> dict:
    """
    Runs simulation based on the specified mode in config.
    
    This function is the main entry point for simulations. It handles
    running the appropriate simulation mode(s) based on the configuration.
    
    Args:
        config: Configuration object containing simulation parameters
    
    Returns:
        dict: Dictionary containing simulation results for each mode:
            - 'locus': Tree sequence for locus mode (if requested)
            - 'chromosome': Tree sequence for chromosome mode (if requested)
    
    Note:
        The function will run both modes if simulation_mode is 'both'.
    """
    results = {}
    
    # Run locus simulation if requested
    if config.simulation_mode in ['locus', 'both']:
        print("Simulating independent non-recombining loci...")
        ts_locus = simulate_locus(config, config.locus_length)
        results['locus'] = ts_locus
        print(f"Generated {config.n_loci} non-recombining loci")
    
    # Run chromosome simulation if requested
    if config.simulation_mode in ['chromosome', 'both']:
        print("Simulating recombining chromosome...")
        ts_chrom = simulate_chromosome(config)
        results['chromosome'] = ts_chrom
        print(f"Generated chromosome of length {config.chromosome_length}")
    
    return results

def save_tree_sequences(
    tree_sequences: Union[msprime.TreeSequence, List[msprime.TreeSequence]],
    output_dir: str,
    prefix: str = "simulation"
):
    """
    Save tree sequences to files.
    
    This function saves the simulated tree sequences to files in the
    specified output directory. For locus mode, it saves multiple files
    (one per locus), while for chromosome mode it saves a single file.
    
    Args:
        tree_sequences: Single tree sequence or list of tree sequences
        output_dir: Directory to save files
        prefix: Prefix for output files
    
    Note:
        The output files are saved in .ts format, which is the native
        format for msprime tree sequences.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    if isinstance(tree_sequences, list):
        for i, ts in enumerate(tree_sequences):
            output_file = os.path.join(output_dir, f"{prefix}_locus_{i}.ts")
            ts.dump(output_file)
    else:
        output_file = os.path.join(output_dir, f"{prefix}_chromosome.ts")
        tree_sequences.dump(output_file)
    
    logger.info(f"Saved tree sequences to {output_dir}")

def main():
    """
    Command-line interface for simulation.
    
    This function provides a command-line interface for running simulations.
    It parses command-line arguments and runs the appropriate simulation mode.
    
    Usage:
        python -m twisstntern.simulation -m MODE -c CONFIG -o OUTPUT [options]
    
    Arguments:
        -m, --mode: Simulation mode (locus, chromosome, or both)
        -c, --config: Path to configuration file
        -o, --output: Output directory
        --n-ind: Number of haploid individuals per population
        --n-loci: Number of loci for locus mode
        --rec-rate: Recombination rate for chromosome mode
        --seq-length: Sequence length for chromosome mode
        --mutation-rate: Mutation rate
        --seed: Random seed
    """
    import argparse
    
    parser = argparse.ArgumentParser(description="Simulate demographic scenarios")
    parser.add_argument("-m", "--mode", choices=["locus", "chromosome", "both"], required=True,
                      help="Simulation mode")
    parser.add_argument("-c", "--config", required=True,
                      help="Path to configuration file")
    parser.add_argument("-o", "--output", required=True,
                      help="Output directory")
    parser.add_argument("--n-ind", type=int, default=20,
                      help="Number of haploid individuals per population")
    parser.add_argument("--n-loci", type=int, default=10000,
                      help="Number of loci for locus mode")
    parser.add_argument("--rec-rate", type=float, default=1e-8,
                      help="Recombination rate for chromosome mode")
    parser.add_argument("--seq-length", type=int, default=1000000,
                      help="Sequence length for chromosome mode")
    parser.add_argument("--mutation-rate", type=float, default=1e-8,
                      help="Mutation rate")
    parser.add_argument("--seed", type=int,
                      help="Random seed")
    
    args = parser.parse_args()
    
    try:
        config = Config(args.config)
        config.simulation_mode = args.mode
        config.n_loci = args.n_loci
        config.rec_rate = args.rec_rate
        config.seq_length = args.seq_length
        config.seed = args.seed
        
        results = run_simulation(config)
        
        for mode, ts in results.items():
            save_tree_sequences(ts, args.output, mode)
        
    except Exception as e:
        logger.error(f"Error during simulation: {str(e)}")
        raise

if __name__ == "__main__":
    main() 