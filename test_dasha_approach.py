#!/usr/bin/env python3
"""
Test script to compare our original tree processing approach with the DaSh-bash approach.
This will help us see which method produces more reasonable topology weights for twisst.
"""

import sys
import os
import msprime
import pandas as pd
import numpy as np

# Add the current directory to path so we can import our modules
sys.path.insert(0, '.')

from twisstntern.tree_processing import (
    ts_to_twisst_weights,
    dasha_ts_to_twisst_weights,
    newick_to_twisst_weights
)

# Set up parameters matching generate_trees.py
nP1 = 10  
nP2 = 10  
nP3 = 10  
n0 = 10   

t1 = 100  
t2 = 200  
t3 = 300  

mig_rate = 0.00001

Ne = 1000

def create_demography():
    """Create the same demography as in generate_trees.py"""
    demography = msprime.Demography()

    # initializing populations
    demography.add_population(name="O", initial_size=Ne)
    demography.add_population(name="P1", initial_size=Ne)
    demography.add_population(name="P2", initial_size=Ne)
    demography.add_population(name="P3", initial_size=Ne)
    demography.add_population(name="P13", initial_size=Ne)
    demography.add_population(name="P123", initial_size=Ne)
    demography.add_population(name="ANC", initial_size=Ne)

    # adding split times
    demography.add_population_split(time=t1, derived=["P1", "P2"], ancestral="P13")
    demography.add_population_split(time=t2, derived=["P13", "P3"], ancestral="P123")
    demography.add_population_split(time=t3, derived=["P123", "O"], ancestral="ANC")

    # setting up gene flow
    demography.set_migration_rate("P2", "P3", mig_rate)
    
    return demography

def test_approaches():
    """Test both approaches and compare results"""
    print("=" * 70)
    print("COMPARING ORIGINAL VS DASHA APPROACH")
    print("=" * 70)
    
    # Generate test data
    demography = create_demography()
    
    print("Generating test TreeSequences...")
    num_replicates = 10  # Small number for testing
    
    genealogies = msprime.sim_ancestry(
        samples={"O": n0, "P1": nP1, "P2": nP2, "P3": nP3},
        demography=demography,
        num_replicates=num_replicates,
        ploidy=1,  # Use haploid samples like DaSh-bash approach
    )
    
    # Convert to list so we can use it twice
    genealogies = list(genealogies)
    print(f"Generated {len(genealogies)} TreeSequences")
    
    # Test 1: Original approach
    print("\n" + "-" * 50)
    print("TESTING ORIGINAL APPROACH")
    print("-" * 50)
    
    try:
        weights_original = ts_to_twisst_weights(
            genealogies,
            outgroup="0",  # Population O has ID 0
            verbose=True
        )
        
        print(f"Original approach results:")
        print(f"  Shape: {weights_original.shape}")
        print(f"  Columns: {list(weights_original.columns)}")
        print(f"  Mean weights: {weights_original.mean()}")
        print(f"  Weight ranges:")
        for col in weights_original.columns:
            print(f"    {col}: {weights_original[col].min():.3f} - {weights_original[col].max():.3f}")
        
        # Check for suspicious patterns
        print(f"  Suspicious patterns:")
        print(f"    All zeros rows: {(weights_original.sum(axis=1) == 0).sum()}")
        print(f"    Extreme weights (>0.9): {(weights_original.max(axis=1) > 0.9).sum()}")
        
    except Exception as e:
        print(f"Original approach failed: {e}")
        weights_original = None
    
    # Test 2: DaSh-bash approach
    print("\n" + "-" * 50)
    print("TESTING DASHA APPROACH")
    print("-" * 50)
    
    try:
        weights_dasha = dasha_ts_to_twisst_weights(
            genealogies,
            nP1=nP1, nP2=nP2, nP3=nP3, n0=n0,
            verbose=True
        )
        
        print(f"DaSh-bash approach results:")
        print(f"  Shape: {weights_dasha.shape}")
        print(f"  Columns: {list(weights_dasha.columns)}")
        print(f"  Mean weights: {weights_dasha.mean()}")
        print(f"  Weight ranges:")
        for col in weights_dasha.columns:
            print(f"    {col}: {weights_dasha[col].min():.3f} - {weights_dasha[col].max():.3f}")
        
        # Check for suspicious patterns
        print(f"  Suspicious patterns:")
        print(f"    All zeros rows: {(weights_dasha.sum(axis=1) == 0).sum()}")
        print(f"    Extreme weights (>0.9): {(weights_dasha.max(axis=1) > 0.9).sum()}")
        
    except Exception as e:
        print(f"DaSh-bash approach failed: {e}")
        weights_dasha = None
    
    # Test 3: Test with the generated Newick file
    print("\n" + "-" * 50)
    print("TESTING WITH GENERATED NEWICK FILE")
    print("-" * 50)
    
    newick_file = "tree files/dasha_approach_for_twisst.newick"
    if os.path.exists(newick_file):
        try:
            # Read the newick file
            with open(newick_file, 'r') as f:
                newick_trees = [line.strip() for line in f if line.strip()]
            
            print(f"Loaded {len(newick_trees)} trees from {newick_file}")
            
            weights_newick = newick_to_twisst_weights(
                newick_trees=newick_trees,
                taxon_names=["P1", "P2", "P3", "O"],
                outgroup="O",
                verbose=True
            )
            
            print(f"Newick file approach results:")
            print(f"  Shape: {weights_newick.shape}")
            print(f"  Columns: {list(weights_newick.columns)}")
            print(f"  Mean weights: {weights_newick.mean()}")
            print(f"  Weight ranges:")
            for col in weights_newick.columns:
                print(f"    {col}: {weights_newick[col].min():.3f} - {weights_newick[col].max():.3f}")
            
            # Check for suspicious patterns
            print(f"  Suspicious patterns:")
            print(f"    All zeros rows: {(weights_newick.sum(axis=1) == 0).sum()}")
            print(f"    Extreme weights (>0.9): {(weights_newick.max(axis=1) > 0.9).sum()}")
            
        except Exception as e:
            print(f"Newick file approach failed: {e}")
            weights_newick = None
    else:
        print(f"Newick file not found: {newick_file}")
        weights_newick = None
    
    # Comparison
    print("\n" + "=" * 70)
    print("COMPARISON SUMMARY")
    print("=" * 70)
    
    if weights_original is not None:
        print(f"Original approach: {weights_original.shape[0]} trees processed")
        print(f"  Mean weight distribution: {weights_original.mean().values}")
        
    if weights_dasha is not None:
        print(f"DaSh-bash approach: {weights_dasha.shape[0]} trees processed")
        print(f"  Mean weight distribution: {weights_dasha.mean().values}")
        
    if weights_newick is not None:
        print(f"Newick file approach: {weights_newick.shape[0]} trees processed")
        print(f"  Mean weight distribution: {weights_newick.mean().values}")
    
    # Expected behavior for this demographic model
    print(f"\nExpected behavior for this demographic model:")
    print(f"  - Population topology: ((P1,P2),P3),O")
    print(f"  - With migration P2<->P3, we should see some weight on other topologies")
    print(f"  - No single topology should dominate completely")
    print(f"  - Weights should be reasonably distributed, not extreme")
    
    return weights_original, weights_dasha, weights_newick


if __name__ == "__main__":
    test_approaches() 