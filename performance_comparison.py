#!/usr/bin/env python3
"""
Performance comparison script: Simplified approach vs ts_to_twisst_weights_old_old

This script tests both approaches with varying numbers of TreeSequences
to see which is faster and more memory efficient.
"""

import time
import sys
import gc
import psutil
import os
import pandas as pd
import numpy as np
import msprime
import tskit
from pathlib import Path

# Add the simulate directory to path to import the functions
sys.path.append('twisstntern_simulate')
from ts_processing import ts_to_twisst_weights_old_old, ts_chromosome_to_twisst_weights

def get_memory_usage():
    """Get current memory usage in MB."""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024

def ts_to_twisst_weights_simplified(input_data, outgroup=None, output_file=None, verbose=False, **kwargs):
    """
    My simplified approach - true streaming without topology reconciliation.
    """
    if isinstance(input_data, tskit.TreeSequence):
        return ts_chromosome_to_twisst_weights(input_data, outgroup=outgroup, 
                                               output_file=output_file, verbose=verbose)
    
    all_dfs = []
    for i, ts in enumerate(input_data):  # True streaming - don't convert to list!
        if verbose:
            print(f"Processing TreeSequence {i+1}...")
        
        df = ts_chromosome_to_twisst_weights(ts, outgroup=outgroup, 
                                           output_file=None, verbose=False)
        all_dfs.append(df)
    
    combined_df = pd.concat(all_dfs, ignore_index=True)
    
    if output_file:
        combined_df.to_csv(output_file, index=False, float_format='%.3f')
        if verbose:
            print(f"Results saved to: {output_file}")
    
    return combined_df

def generate_test_tree_sequences(n_sequences, sequence_length=1000000, recombination_rate=1e-8):
    """
    Generate a list of TreeSequence objects for testing.
    """
    print(f"Generating {n_sequences} test TreeSequences...")
    
    # Set up a simple 4-population demography
    demography = msprime.Demography()
    demography.add_population(name="O", initial_size=1000)   # Outgroup
    demography.add_population(name="P1", initial_size=1000)
    demography.add_population(name="P2", initial_size=1000)
    demography.add_population(name="P3", initial_size=1000)
    demography.add_population(name="P12", initial_size=1000)
    demography.add_population(name="P123", initial_size=1000)
    demography.add_population(name="ANC", initial_size=1000)

    # Split times
    demography.add_population_split(time=100, derived=["P1", "P2"], ancestral="P12")
    demography.add_population_split(time=200, derived=["P12", "P3"], ancestral="P123")
    demography.add_population_split(time=300, derived=["P123", "O"], ancestral="ANC")

    tree_sequences = []
    for i in range(n_sequences):
        ts = msprime.sim_ancestry(
            samples={"O": 5, "P1": 5, "P2": 5, "P3": 5},
            demography=demography,
            sequence_length=sequence_length,
            recombination_rate=recombination_rate,
            ploidy=1,
            random_seed=42 + i  # Different seed for each
        )
        tree_sequences.append(ts)
    
    print(f"Generated {len(tree_sequences)} TreeSequences")
    return tree_sequences

def time_function(func, *args, **kwargs):
    """Time a function call and return (result, elapsed_time, memory_used)."""
    gc.collect()  # Clean up before measurement
    start_memory = get_memory_usage()
    start_time = time.time()
    
    result = func(*args, **kwargs)
    
    end_time = time.time()
    end_memory = get_memory_usage()
    
    elapsed_time = end_time - start_time
    memory_used = end_memory - start_memory
    
    return result, elapsed_time, memory_used

def run_performance_test(n_sequences_list=[5, 10, 20]):
    """
    Run performance comparison for different numbers of TreeSequences.
    """
    results = []
    
    print("="*60)
    print("PERFORMANCE COMPARISON: Simplified vs Old_Old")
    print("="*60)
    
    for n_seq in n_sequences_list:
        print(f"\n--- Testing with {n_seq} TreeSequences ---")
        
        # Generate test data
        tree_sequences = generate_test_tree_sequences(n_seq)
        
        # Test 1: Old_Old approach (converts generator to list)
        print("\n1. Testing old_old approach...")
        result_old, time_old, memory_old = time_function(
            ts_to_twisst_weights_old_old,
            tree_sequences,  # Pass as list
            outgroup="0",
            verbose=False
        )
        
        # Test 2: Simplified approach (true streaming)
        print("2. Testing simplified approach...")
        result_simplified, time_simplified, memory_simplified = time_function(
            ts_to_twisst_weights_simplified,
            iter(tree_sequences),  # Pass as generator
            outgroup="0",
            verbose=False
        )
        
        # Compare results
        print(f"\n--- Results for {n_seq} TreeSequences ---")
        print(f"Old_Old approach:")
        print(f"  Time: {time_old:.3f} seconds")
        print(f"  Memory: {memory_old:.1f} MB")
        print(f"  Result shape: {result_old.shape}")
        
        print(f"Simplified approach:")
        print(f"  Time: {time_simplified:.3f} seconds")
        print(f"  Memory: {memory_simplified:.1f} MB")
        print(f"  Result shape: {result_simplified.shape}")
        
        # Calculate improvements
        time_improvement = ((time_old - time_simplified) / time_old) * 100
        memory_improvement = ((memory_old - memory_simplified) / memory_old) * 100 if memory_old > 0 else 0
        
        print(f"\nImprovements:")
        print(f"  Time: {time_improvement:+.1f}% ({'faster' if time_improvement > 0 else 'slower'})")
        print(f"  Memory: {memory_improvement:+.1f}% ({'less' if memory_improvement > 0 else 'more'})")
        
        # Store results
        results.append({
            'n_sequences': n_seq,
            'old_time': time_old,
            'old_memory': memory_old,
            'old_shape': result_old.shape,
            'simplified_time': time_simplified,
            'simplified_memory': memory_simplified,
            'simplified_shape': result_simplified.shape,
            'time_improvement': time_improvement,
            'memory_improvement': memory_improvement
        })
        
        # Quick sanity check - do they produce similar number of rows?
        if abs(result_old.shape[0] - result_simplified.shape[0]) > n_seq * 10:  # Allow some variance
            print(f"⚠️  WARNING: Results have very different shapes!")
    
    return results

def generate_summary_report(results):
    """Generate a summary report of the performance comparison."""
    print("\n" + "="*60)
    print("SUMMARY REPORT")
    print("="*60)
    
    df = pd.DataFrame(results)
    
    print("\nPerformance Summary:")
    print(f"{'N_Seq':<8} {'Old(s)':<8} {'New(s)':<8} {'Time%':<8} {'Old(MB)':<8} {'New(MB)':<8} {'Mem%':<8}")
    print("-" * 60)
    
    for _, row in df.iterrows():
        print(f"{row['n_sequences']:<8} "
              f"{row['old_time']:<8.2f} "
              f"{row['simplified_time']:<8.2f} "
              f"{row['time_improvement']:<+7.1f}% "
              f"{row['old_memory']:<8.1f} "
              f"{row['simplified_memory']:<8.1f} "
              f"{row['memory_improvement']:<+7.1f}%")
    
    print(f"\nAverage improvements:")
    print(f"  Time: {df['time_improvement'].mean():+.1f}%")
    print(f"  Memory: {df['memory_improvement'].mean():+.1f}%")
    
    # Conclusions
    print(f"\n{'CONCLUSIONS:':<20}")
    if df['time_improvement'].mean() > 5:
        print("✅ Simplified approach is significantly faster")
    elif df['time_improvement'].mean() > 0:
        print("✅ Simplified approach is slightly faster")
    else:
        print("❌ Old_old approach is faster")
    
    if df['memory_improvement'].mean() > 5:
        print("✅ Simplified approach uses significantly less memory")
    elif df['memory_improvement'].mean() > 0:
        print("✅ Simplified approach uses slightly less memory")
    else:
        print("❌ Old_old approach uses less memory")

def main():
    """Main function to run the performance comparison."""
    print("Starting performance comparison...")
    print(f"Python version: {sys.version}")
    print(f"Initial memory usage: {get_memory_usage():.1f} MB")
    
    # Run tests with increasing numbers of TreeSequences
    test_sizes = [3, 5, 10]  # Start small for testing
    
    try:
        results = run_performance_test(test_sizes)
        generate_summary_report(results)
        
    except Exception as e:
        print(f"Error during testing: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 