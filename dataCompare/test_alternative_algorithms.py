#!/usr/bin/env python3
"""
Test alternative Wasserstein algorithms to find methods with less bias.
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path
import time
import matplotlib.pyplot as plt
import seaborn as sns

# Add parent directory for imports
sys.path.append(str(Path(__file__).parent.parent))

try:
    import ot
    from scipy.spatial.distance import pdist, squareform
    from scipy.stats import wasserstein_distance as scipy_wasserstein
    HAS_ADVANCED_OT = True
except ImportError as e:
    print(f"Warning: Some libraries not available: {e}")
    HAS_ADVANCED_OT = False

from fast_topology_metrics import (
    load_and_sample_data,
    compute_wasserstein_euclidean_cpu,
    compute_wasserstein_kl_cpu
)

def compute_alternative_wasserstein_methods(data1, data2, max_points=2000):
    """Test various alternative Wasserstein computation methods."""
    
    # Sample data if needed
    if len(data1) > max_points:
        data1 = data1.sample(n=max_points, random_state=42)
    if len(data2) > max_points:
        data2 = data2.sample(n=max_points, random_state=42)
    
    data1_arr = data1[['T1', 'T2', 'T3']].values
    data2_arr = data2[['T1', 'T2', 'T3']].values
    
    methods = {}
    
    # 1. Our exact method (baseline)
    try:
        start_time = time.time()
        exact_eucl = compute_wasserstein_euclidean_cpu(data1, data2, verbose=False)
        exact_time = time.time() - start_time
        methods['exact_euclidean'] = {'value': exact_eucl, 'time': exact_time}
    except Exception as e:
        print(f"Exact Euclidean failed: {e}")
        methods['exact_euclidean'] = {'value': np.nan, 'time': np.nan}
    
    # 2. Network Simplex (exact but different algorithm)
    try:
        start_time = time.time()
        a = np.ones(len(data1_arr)) / len(data1_arr)
        b = np.ones(len(data2_arr)) / len(data2_arr)
        M = ot.dist(data1_arr, data2_arr)
        
        # Use network simplex instead of EMD
        network_eucl = ot.emd2(a, b, M, solver='network_simplex')
        network_time = time.time() - start_time
        methods['network_simplex'] = {'value': network_eucl, 'time': network_time}
    except Exception as e:
        print(f"Network simplex failed: {e}")
        methods['network_simplex'] = {'value': np.nan, 'time': np.nan}
    
    # 3. Unbalanced OT (less constraints)
    try:
        start_time = time.time()
        a = np.ones(len(data1_arr)) / len(data1_arr)
        b = np.ones(len(data2_arr)) / len(data2_arr)
        M = ot.dist(data1_arr, data2_arr)
        
        # Unbalanced with different regularization
        unbalanced_eucl = ot.unbalanced.sinkhorn_unbalanced2(a, b, M, reg=0.01, reg_m=0.1)
        unbalanced_time = time.time() - start_time
        methods['unbalanced_sinkhorn'] = {'value': unbalanced_eucl, 'time': unbalanced_time}
    except Exception as e:
        print(f"Unbalanced Sinkhorn failed: {e}")
        methods['unbalanced_sinkhorn'] = {'value': np.nan, 'time': np.nan}
    
    # 4. Stabilized Sinkhorn (better numerical stability)
    try:
        start_time = time.time()
        a = np.ones(len(data1_arr)) / len(data1_arr)
        b = np.ones(len(data2_arr)) / len(data2_arr)
        M = ot.dist(data1_arr, data2_arr)
        
        stabilized_eucl = ot.sinkhorn2(a, b, M, reg=0.01, method='sinkhorn_stabilized')
        stabilized_time = time.time() - start_time
        methods['stabilized_sinkhorn'] = {'value': stabilized_eucl, 'time': stabilized_time}
    except Exception as e:
        print(f"Stabilized Sinkhorn failed: {e}")
        methods['stabilized_sinkhorn'] = {'value': np.nan, 'time': np.nan}
    
    # 5. Epsilon-scaling method
    try:
        start_time = time.time()
        a = np.ones(len(data1_arr)) / len(data1_arr)
        b = np.ones(len(data2_arr)) / len(data2_arr)
        M = ot.dist(data1_arr, data2_arr)
        
        epsilon_eucl = ot.sinkhorn2(a, b, M, reg=0.01, method='sinkhorn_epsilon_scaling')
        epsilon_time = time.time() - start_time
        methods['epsilon_scaling'] = {'value': epsilon_eucl, 'time': epsilon_time}
    except Exception as e:
        print(f"Epsilon scaling failed: {e}")
        methods['epsilon_scaling'] = {'value': np.nan, 'time': np.nan}
    
    # 6. Barycenters approach (approximate)
    try:
        start_time = time.time()
        
        # Compute barycenters and use them for distance approximation
        # This is experimental
        weights = [0.5, 0.5]
        measures = [data1_arr, data2_arr]
        
        # This is a heuristic approximation
        barycenter_eucl = np.linalg.norm(np.mean(data1_arr, axis=0) - np.mean(data2_arr, axis=0))
        barycenter_time = time.time() - start_time
        methods['barycenter_approx'] = {'value': barycenter_eucl, 'time': barycenter_time}
    except Exception as e:
        print(f"Barycenter approximation failed: {e}")
        methods['barycenter_approx'] = {'value': np.nan, 'time': np.nan}
    
    # 7. Sliced Wasserstein (1D projections)
    try:
        start_time = time.time()
        n_projections = 50
        
        sliced_distances = []
        for _ in range(n_projections):
            # Random direction
            theta = np.random.randn(3)
            theta = theta / np.linalg.norm(theta)
            
            # Project data
            proj1 = data1_arr @ theta
            proj2 = data2_arr @ theta
            
            # 1D Wasserstein (much faster)
            sliced_dist = scipy_wasserstein(proj1, proj2)
            sliced_distances.append(sliced_dist)
        
        sliced_eucl = np.mean(sliced_distances)
        sliced_time = time.time() - start_time
        methods['sliced_wasserstein'] = {'value': sliced_eucl, 'time': sliced_time}
    except Exception as e:
        print(f"Sliced Wasserstein failed: {e}")
        methods['sliced_wasserstein'] = {'value': np.nan, 'time': np.nan}
    
    # 8. Low-rank approximation Sinkhorn
    try:
        start_time = time.time()
        a = np.ones(len(data1_arr)) / len(data1_arr)
        b = np.ones(len(data2_arr)) / len(data2_arr)
        M = ot.dist(data1_arr, data2_arr)
        
        # Use lower rank approximation
        lowrank_eucl = ot.lowrank_sinkhorn(a, b, M, reg=0.01, rank=50)
        lowrank_time = time.time() - start_time
        methods['lowrank_sinkhorn'] = {'value': lowrank_eucl, 'time': lowrank_time}
    except Exception as e:
        print(f"Low-rank Sinkhorn failed: {e}")
        methods['lowrank_sinkhorn'] = {'value': np.nan, 'time': np.nan}
    
    return methods

def test_alternative_algorithms():
    """Test alternative Wasserstein algorithms for bias and speed."""
    
    print("="*80)
    print("ALTERNATIVE WASSERSTEIN ALGORITHMS TEST")
    print("="*80)
    
    # Load test datasets
    print("Loading test datasets...")
    data1 = load_and_sample_data("migration_topology_weights.csv", sample_size=1000, random_state=42)
    data2 = load_and_sample_data("NoMigration_topology_weights.csv", sample_size=1000, random_state=42)
    
    print(f"Dataset 1: {len(data1)} points")
    print(f"Dataset 2: {len(data2)} points")
    
    # Test different sample sizes
    sample_sizes = [500, 1000, 2000]
    all_results = []
    
    for sample_size in sample_sizes:
        print(f"\n{'='*60}")
        print(f"TESTING WITH {sample_size} POINTS")
        print(f"{'='*60}")
        
        # Get reference exact value
        data1_sample = data1.sample(n=min(sample_size, len(data1)), random_state=42)
        data2_sample = data2.sample(n=min(sample_size, len(data2)), random_state=42)
        
        exact_eucl = compute_wasserstein_euclidean_cpu(data1_sample, data2_sample, verbose=False)
        print(f"Reference exact Euclidean: {exact_eucl:.6f}")
        
        # Test all methods
        methods = compute_alternative_wasserstein_methods(data1_sample, data2_sample, max_points=sample_size)
        
        print(f"\nResults:")
        print(f"{'Method':<20} {'Value':<10} {'Error':<8} {'Rel Error':<10} {'Time':<8}")
        print("-" * 65)
        
        for method_name, result in methods.items():
            value = result['value']
            time_taken = result['time']
            
            if not np.isnan(value) and not np.isnan(exact_eucl):
                error = abs(value - exact_eucl)
                rel_error = (error / exact_eucl) * 100
                
                print(f"{method_name:<20} {value:<10.4f} {error:<8.4f} {rel_error:<10.1f}% {time_taken:<8.3f}s")
                
                all_results.append({
                    'sample_size': sample_size,
                    'method': method_name,
                    'value': value,
                    'exact_value': exact_eucl,
                    'error': error,
                    'rel_error': rel_error,
                    'time': time_taken
                })
            else:
                print(f"{method_name:<20} {'FAILED':<10} {'-':<8} {'-':<10} {time_taken:<8.3f}s")
    
    # Convert to DataFrame and analyze
    results_df = pd.DataFrame(all_results)
    
    if len(results_df) > 0:
        # Save results
        results_df.to_csv('alternative_algorithms_results.csv', index=False)
        
        # Analysis
        print(f"\n{'='*80}")
        print("ALTERNATIVE ALGORITHMS ANALYSIS")
        print(f"{'='*80}")
        
        # Best methods by accuracy
        print(f"\nBEST METHODS BY ACCURACY:")
        best_accuracy = results_df.loc[results_df['rel_error'].abs().idxmin()]
        print(f"  {best_accuracy['method']} - {best_accuracy['rel_error']:+.1f}% error, "
              f"{best_accuracy['time']:.3f}s")
        
        # Best methods by speed
        print(f"\nFASTEST METHODS:")
        fastest = results_df.loc[results_df['time'].idxmin()]
        print(f"  {fastest['method']} - {fastest['time']:.3f}s, "
              f"{fastest['rel_error']:+.1f}% error")
        
        # Good compromises
        good_methods = results_df[(results_df['rel_error'].abs() < 20) & (results_df['time'] < 1.0)]
        if len(good_methods) > 0:
            print(f"\nGOOD COMPROMISES (<20% error, <1s):")
            for _, row in good_methods.iterrows():
                print(f"  {row['method']} - {row['rel_error']:+.1f}% error, {row['time']:.3f}s")
        
        # Create visualizations
        create_algorithm_visualizations(results_df)
    
    return results_df

def create_algorithm_visualizations(results_df):
    """Create visualizations comparing alternative algorithms."""
    
    if len(results_df) == 0:
        print("No data to visualize!")
        return
    
    # Set up plotting
    plt.style.use('default')
    sns.set_palette("Set2")
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Alternative Wasserstein Algorithms Comparison', fontsize=16)
    
    # 1. Error comparison
    ax = axes[0, 0]
    for size in results_df['sample_size'].unique():
        subset = results_df[results_df['sample_size'] == size]
        ax.scatter(subset['method'], subset['rel_error'], label=f'{size} points', alpha=0.7, s=80)
    ax.set_ylabel('Relative Error (%)')
    ax.set_title('Accuracy Comparison')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='x', rotation=45)
    ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
    
    # 2. Time comparison
    ax = axes[0, 1]
    for size in results_df['sample_size'].unique():
        subset = results_df[results_df['sample_size'] == size]
        ax.scatter(subset['method'], subset['time'], label=f'{size} points', alpha=0.7, s=80)
    ax.set_ylabel('Computation Time (s)')
    ax.set_title('Speed Comparison')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='x', rotation=45)
    ax.set_yscale('log')
    
    # 3. Error vs Time tradeoff
    ax = axes[1, 0]
    methods = results_df['method'].unique()
    colors = plt.cm.Set3(np.linspace(0, 1, len(methods)))
    
    for method, color in zip(methods, colors):
        subset = results_df[results_df['method'] == method]
        ax.scatter(subset['time'], subset['rel_error'].abs(), 
                  label=method, alpha=0.7, s=80, color=color)
    
    ax.set_xlabel('Computation Time (s)')
    ax.set_ylabel('|Relative Error| (%)')
    ax.set_title('Error vs Time Tradeoff')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    
    # 4. Method performance heatmap
    ax = axes[1, 1]
    pivot_data = results_df.pivot_table(values='rel_error', 
                                       index='method', columns='sample_size', aggfunc='mean')
    
    sns.heatmap(pivot_data, annot=True, fmt='.1f', cmap='RdBu_r', center=0, ax=ax)
    ax.set_title('Error Heatmap by Method and Sample Size (%)')
    ax.set_xlabel('Sample Size')
    ax.set_ylabel('Method')
    
    plt.tight_layout()
    plt.savefig('alternative_algorithms_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Visualization saved to: alternative_algorithms_comparison.png")

if __name__ == "__main__":
    if not HAS_ADVANCED_OT:
        print("Error: Required libraries not available. Please install: pip install POT scipy")
        sys.exit(1)
    
    results = test_alternative_algorithms() 