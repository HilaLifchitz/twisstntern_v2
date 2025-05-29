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
from typing import List, Dict, Union, Tuple
import ete3
from pathlib import Path

# Add the external directory to the Python path
EXTERNAL_DIR = Path(__file__).parent / 'external'
sys.path.append(str(EXTERNAL_DIR))

# Import twisst functions
try:
    from twisst import (
        load_trees,
        get_taxa,
        get_topologies,
        get_weights
    )
except ImportError:
    print("Warning: twisst.py not found in external directory.")
    print("Please download twisst.py from https://github.com/simonhmartin/twisst")
    print("and place it in the twisstntern/external directory.")

def read_tree_sequence(file_path: str, format: str = 'newick') -> List[str]:
    """
    Read a tree sequence file and return a list of trees.
    
    Args:
        file_path (str): Path to the tree sequence file
        format (str): Format of the tree file ('newick', 'nexus', etc.)
    
    Returns:
        List[str]: List of trees in newick format
    """
    try:
        with open(file_path, 'r') as f:
            if format == 'newick':
                # Read newick format (one tree per line)
                trees = [line.strip() for line in f if line.strip()]
            else:
                raise ValueError(f"Unsupported tree format: {format}")
        return trees
    except Exception as e:
        raise Exception(f"Error reading tree sequence file: {str(e)}")

def process_trees_with_twisst(trees: List[str],
                            taxa_groups: Dict[str, List[str]],
                            window_size: int = 1000,
                            window_step: int = 100) -> pd.DataFrame:
    """
    Process trees using twisst to generate topology weights.
    
    Args:
        trees (List[str]): List of trees in newick format
        taxa_groups (Dict[str, List[str]]): Dictionary mapping taxon names to lists of tip labels
        window_size (int): Size of the sliding window
        window_step (int): Step size for the sliding window
    
    Returns:
        pd.DataFrame: DataFrame with topology weights
    """
    try:
        # Load trees
        tree_list = load_trees(trees)
        
        # Get unique topologies
        topologies = get_topologies(tree_list, taxa_groups)
        
        # Calculate weights
        weights = get_weights(tree_list, topologies, taxa_groups,
                            window_size=window_size,
                            window_step=window_step)
        
        # Convert to DataFrame
        weights_df = pd.DataFrame(weights)
        
        return weights_df
    except Exception as e:
        raise Exception(f"Error processing trees with twisst: {str(e)}")

def convert_weights_to_ternary(weights_df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert twisst weights to ternary coordinates format.
    
    Args:
        weights_df (pd.DataFrame): DataFrame with topology weights
    
    Returns:
        pd.DataFrame: DataFrame with ternary coordinates
    """
    try:
        # Convert weights to ternary coordinates
        # Assuming weights are in columns: w1, w2, w3
        ternary_df = pd.DataFrame({
            'T1': weights_df['w1'],
            'T2': weights_df['w2'],
            'T3': weights_df['w3']
        })
        
        return ternary_df
    except Exception as e:
        raise Exception(f"Error converting weights to ternary coordinates: {str(e)}")

def process_tree_sequence(file_path: str,
                         taxa_groups: Dict[str, List[str]],
                         format: str = 'newick',
                         window_size: int = 1000,
                         window_step: int = 100) -> pd.DataFrame:
    """
    Process a tree sequence file and return ternary coordinates.
    
    Args:
        file_path (str): Path to the tree sequence file
        taxa_groups (Dict[str, List[str]]): Dictionary mapping taxon names to lists of tip labels
        format (str): Format of the tree file
        window_size (int): Size of the sliding window for twisst
        window_step (int): Step size for the sliding window
    
    Returns:
        pd.DataFrame: DataFrame with ternary coordinates ready for analysis
    """
    # Read tree sequence
    trees = read_tree_sequence(file_path, format)
    
    # Process trees with twisst
    weights_df = process_trees_with_twisst(trees, taxa_groups, window_size, window_step)
    
    # Convert to ternary coordinates
    ternary_df = convert_weights_to_ternary(weights_df)
    
    return ternary_df 