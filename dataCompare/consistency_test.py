#!/usr/bin/env python3
"""
Consistency Test for Topology Metrics Algorithms

Test script to compare consistency between:
1. Normal exact Wasserstein algorithm
2. Sinkhorn approximation with 200 points
3. Sinkhorn approximation with 10,000 points

Creates 20 pairs of small datasets and compares results.
"""

import sys
import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Add parent directory for imports
sys.path.append(str(Path(__file__).parent.parent))
from twisstntern.utils import dump_data

# Import our fast metrics functions
from fast_topology_metrics import (
    compute_wasserstein_euclidean_cpu,
    compute_wasserstein_kl_cpu,
    compute_sinkhorn_wasserstein_fast,
    compute_sinkhorn_kl_fast,
    compute_l2_chi2_ultrafast
)

def load_and_sample_data(filepath, sample_size=1000, random_state=42):
    """Load data and create a random sample."""
    print(f"Loading data from {filepath}...")
    
    # Load full dataset
    data_array = dump_data(str(filepath))
    data = pd.DataFrame(data_array, columns=['T1', 'T2', 'T3'])
    
    # Sample the requested number of points
    if len(data) > sample_size:
        data = data.sample(n=sample_size, random_state=random_state)
    
    return data

def create_test_datasets(data1_path, data2_path, n_pairs=20, base_sample_size=1000):
    """Create multiple pairs of test datasets."""
    print(f"Creating {n_pairs} pairs of test datasets...")
    
    # Load the full datasets
    full_data1 = load_and_sample_data(data1_path, sample_size=base_sample_size*n_pairs, random_state=42)
    full_data2 = load_and_sample_data(data2_path, sample_size=base_sample_size*n_pairs, random_state=43)
    
    dataset_pairs = []
    
    for i in range(n_pairs):
        # Create different random samples for each pair
        start_idx = i * base_sample_size
        end_idx = (i + 1) * base_sample_size
        
        data1_sample = full_data1.iloc[start_idx:end_idx].copy().reset_index(drop=True)
        data2_sample = full_data2.iloc[start_idx:end_idx].copy().reset_index(drop=True)
        
        dataset_pairs.append((data1_sample, data2_sample, f"Pair_{i+1:02d}"))
    
    return dataset_pairs

def test_single_pair(data1, data2, pair_name, verbose=True):
    """Test a single pair of datasets with all three methods."""
    if verbose:
        print(f"\nTesting {pair_name} ({len(data1)} vs {len(data2)} points)...")
    
    results = {
        'pair_name': pair_name,
        'n_data1': len(data1),
        'n_data2': len(data2)
    }
    
    # 1. Exact Wasserstein (normal algorithm)
    if verbose:
        print("  1. Normal exact Wasserstein...")
    
    start_time = time.time()
    try:
        exact_eucl = compute_wasserstein_euclidean_cpu(data1, data2, max_points=None, verbose=False)
        exact_kl = compute_wasserstein_kl_cpu(data1, data2, max_points=None, verbose=False)
        exact_time = time.time() - start_time
        
        results.update({
            'exact_wasserstein_eucl': exact_eucl,
            'exact_wasserstein_kl': exact_kl,
            'exact_time': exact_time
        })
        
        if verbose:
            print(f"    Euclidean: {exact_eucl:.6f}, KL: {exact_kl:.6f}, Time: {exact_time:.3f}s")
    
    except Exception as e:
        print(f"    Error in exact method: {e}")
        results.update({
            'exact_wasserstein_eucl': np.nan,
            'exact_wasserstein_kl': np.nan,
            'exact_time': np.nan
        })
    
    # 2. Sinkhorn with 2000 points
    if verbose:
        print("  2. Sinkhorn approximation (2000 points)...")
    
    start_time = time.time()
    try:
        sinkhorn_200_eucl = compute_sinkhorn_wasserstein_fast(data1, data2, max_points=2000, reg=0.01, verbose=False)
        sinkhorn_200_kl = compute_sinkhorn_kl_fast(data1, data2, max_points=2000, reg=0.005, verbose=False)
        sinkhorn_200_time = time.time() - start_time
        
        results.update({
            'sinkhorn_200_wasserstein_eucl': sinkhorn_200_eucl,
            'sinkhorn_200_wasserstein_kl': sinkhorn_200_kl,
            'sinkhorn_200_time': sinkhorn_200_time
        })
        
        if verbose:
            print(f"    Euclidean: {sinkhorn_200_eucl:.6f}, KL: {sinkhorn_200_kl:.6f}, Time: {sinkhorn_200_time:.3f}s")
    
    except Exception as e:
        print(f"    Error in Sinkhorn 2000: {e}")
        results.update({
            'sinkhorn_200_wasserstein_eucl': np.nan,
            'sinkhorn_200_wasserstein_kl': np.nan,
            'sinkhorn_200_time': np.nan
        })
    
    # 3. Sinkhorn with 10,000 points
    if verbose:
        print("  3. Sinkhorn approximation (10,000 points)...")
    
    start_time = time.time()
    try:
        sinkhorn_10k_eucl = compute_sinkhorn_wasserstein_fast(data1, data2, max_points=10000, reg=0.01, verbose=False)
        sinkhorn_10k_kl = compute_sinkhorn_kl_fast(data1, data2, max_points=10000, reg=0.005, verbose=False)
        sinkhorn_10k_time = time.time() - start_time
        
        results.update({
            'sinkhorn_10k_wasserstein_eucl': sinkhorn_10k_eucl,
            'sinkhorn_10k_wasserstein_kl': sinkhorn_10k_kl,
            'sinkhorn_10k_time': sinkhorn_10k_time
        })
        
        if verbose:
            print(f"    Euclidean: {sinkhorn_10k_eucl:.6f}, KL: {sinkhorn_10k_kl:.6f}, Time: {sinkhorn_10k_time:.3f}s")
    
    except Exception as e:
        print(f"    Error in Sinkhorn 10k: {e}")
        results.update({
            'sinkhorn_10k_wasserstein_eucl': np.nan,
            'sinkhorn_10k_wasserstein_kl': np.nan,
            'sinkhorn_10k_time': np.nan
        })
    
    # Also compute L2 and chi2 for reference
    try:
        L2_dist, chi2_stat, p_val, dof, n_meaningful = compute_l2_chi2_ultrafast(data1, data2, 0.1, verbose=False)
        results.update({
            'L2_distance': L2_dist,
            'chi2_statistic': chi2_stat,
            'meaningful_triangles': n_meaningful
        })
    except Exception as e:
        print(f"    Error in L2/chi2: {e}")
        results.update({
            'L2_distance': np.nan,
            'chi2_statistic': np.nan,
            'meaningful_triangles': np.nan
        })
    
    return results

def analyze_consistency(results_df):
    """Analyze consistency between different methods."""
    print("\n" + "="*80)
    print("CONSISTENCY ANALYSIS")
    print("="*80)
    
    # Remove rows with NaN values for analysis
    clean_df = results_df.dropna()
    
    if len(clean_df) == 0:
        print("No valid results to analyze!")
        return
    
    print(f"Analyzing {len(clean_df)} successful test pairs...\n")
    
    # 1. Correlation analysis
    print("CORRELATION ANALYSIS:")
    print("-" * 40)
    
    methods = ['exact', 'sinkhorn_200', 'sinkhorn_10k']
    metrics = ['wasserstein_eucl', 'wasserstein_kl']
    
    for metric in metrics:
        print(f"\n{metric.upper()} correlations:")
        
        # Create correlation matrix
        metric_cols = [f"{method}_{metric}" for method in methods]
        available_cols = [col for col in metric_cols if col in clean_df.columns]
        
        if len(available_cols) >= 2:
            corr_matrix = clean_df[available_cols].corr()
            print(corr_matrix.round(4))
            
            # Specific correlations of interest
            if f'exact_{metric}' in available_cols and f'sinkhorn_10k_{metric}' in available_cols:
                corr_exact_10k = clean_df[f'exact_{metric}'].corr(clean_df[f'sinkhorn_10k_{metric}'])
                print(f"  Exact vs Sinkhorn 10k: {corr_exact_10k:.4f}")
            
            if f'sinkhorn_200_{metric}' in available_cols and f'sinkhorn_10k_{metric}' in available_cols:
                corr_200_10k = clean_df[f'sinkhorn_200_{metric}'].corr(clean_df[f'sinkhorn_10k_{metric}'])
                print(f"  Sinkhorn 2000 vs 10k: {corr_200_10k:.4f}")
    
    # 2. Error analysis
    print(f"\n\nERROR ANALYSIS:")
    print("-" * 40)
    
    for metric in metrics:
        exact_col = f'exact_{metric}'
        sinkhorn_10k_col = f'sinkhorn_10k_{metric}'
        sinkhorn_200_col = f'sinkhorn_200_{metric}'
        
        if exact_col in clean_df.columns and sinkhorn_10k_col in clean_df.columns:
            errors_10k = np.abs(clean_df[exact_col] - clean_df[sinkhorn_10k_col])
            rel_errors_10k = errors_10k / clean_df[exact_col] * 100
            
            print(f"\n{metric.upper()} - Sinkhorn 10k vs Exact:")
            print(f"  Mean absolute error: {errors_10k.mean():.6f}")
            print(f"  Max absolute error: {errors_10k.max():.6f}")
            print(f"  Mean relative error: {rel_errors_10k.mean():.2f}%")
            print(f"  Max relative error: {rel_errors_10k.max():.2f}%")
        
        if exact_col in clean_df.columns and sinkhorn_200_col in clean_df.columns:
            errors_200 = np.abs(clean_df[exact_col] - clean_df[sinkhorn_200_col])
            rel_errors_200 = errors_200 / clean_df[exact_col] * 100
            
            print(f"\n{metric.upper()} - Sinkhorn 2000 vs Exact:")
            print(f"  Mean absolute error: {errors_200.mean():.6f}")
            print(f"  Max absolute error: {errors_200.max():.6f}")
            print(f"  Mean relative error: {rel_errors_200.mean():.2f}%")
            print(f"  Max relative error: {rel_errors_200.max():.2f}%")
    
    # 3. Timing analysis
    print(f"\n\nTIMING ANALYSIS:")
    print("-" * 40)
    
    time_cols = ['exact_time', 'sinkhorn_200_time', 'sinkhorn_10k_time']
    available_time_cols = [col for col in time_cols if col in clean_df.columns]
    
    for col in available_time_cols:
        method_name = col.replace('_time', '').replace('_', ' ').title()
        mean_time = clean_df[col].mean()
        std_time = clean_df[col].std()
        print(f"  {method_name}: {mean_time:.3f} ± {std_time:.3f} seconds")
    
    if 'exact_time' in clean_df.columns and 'sinkhorn_10k_time' in clean_df.columns:
        speedup = clean_df['exact_time'].mean() / clean_df['sinkhorn_10k_time'].mean()
        print(f"  Speedup (Sinkhorn 10k vs Exact): {speedup:.1f}x")

def create_visualization(results_df, output_dir='consistency_test_results'):
    """Create visualizations of the consistency analysis."""
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    clean_df = results_df.dropna()
    
    if len(clean_df) == 0:
        print("No valid results for visualization!")
        return
    
    # Set up the plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # 1. Correlation scatter plots
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle('Wasserstein Distance Method Comparison', fontsize=16)
    
    metrics = ['wasserstein_eucl', 'wasserstein_kl']
    
    for i, metric in enumerate(metrics):
        exact_col = f'exact_{metric}'
        sinkhorn_10k_col = f'sinkhorn_10k_{metric}'
        sinkhorn_200_col = f'sinkhorn_200_{metric}'
        
        # Exact vs Sinkhorn 10k
        if exact_col in clean_df.columns and sinkhorn_10k_col in clean_df.columns:
            ax = axes[i, 0]
            ax.scatter(clean_df[exact_col], clean_df[sinkhorn_10k_col], alpha=0.7)
            ax.plot([clean_df[exact_col].min(), clean_df[exact_col].max()], 
                   [clean_df[exact_col].min(), clean_df[exact_col].max()], 'r--', alpha=0.8)
            ax.set_xlabel(f'Exact {metric}')
            ax.set_ylabel(f'Sinkhorn 10k {metric}')
            ax.set_title(f'{metric.replace("_", " ").title()}: Exact vs Sinkhorn 10k')
            
            # Add correlation coefficient
            corr = clean_df[exact_col].corr(clean_df[sinkhorn_10k_col])
            ax.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax.transAxes, 
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        # Sinkhorn 2000 vs Sinkhorn 10k
        if sinkhorn_200_col in clean_df.columns and sinkhorn_10k_col in clean_df.columns:
            ax = axes[i, 1]
            ax.scatter(clean_df[sinkhorn_200_col], clean_df[sinkhorn_10k_col], alpha=0.7, color='orange')
            ax.plot([clean_df[sinkhorn_200_col].min(), clean_df[sinkhorn_200_col].max()], 
                   [clean_df[sinkhorn_200_col].min(), clean_df[sinkhorn_200_col].max()], 'r--', alpha=0.8)
            ax.set_xlabel(f'Sinkhorn 2000 {metric}')
            ax.set_ylabel(f'Sinkhorn 10k {metric}')
            ax.set_title(f'{metric.replace("_", " ").title()}: Sinkhorn 2000 vs 10k')
            
            # Add correlation coefficient
            corr = clean_df[sinkhorn_200_col].corr(clean_df[sinkhorn_10k_col])
            ax.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax.transAxes,
                   bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(output_path / 'method_comparison_scatter.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # 2. Timing comparison
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    time_cols = ['exact_time', 'sinkhorn_200_time', 'sinkhorn_10k_time']
    time_labels = ['Exact', 'Sinkhorn 2000', 'Sinkhorn 10k']
    available_times = []
    available_labels = []
    
    for col, label in zip(time_cols, time_labels):
        if col in clean_df.columns:
            available_times.append(clean_df[col])
            available_labels.append(label)
    
    if available_times:
        ax.boxplot(available_times, labels=available_labels)
        ax.set_ylabel('Computation Time (seconds)')
        ax.set_title('Computation Time Comparison')
        ax.grid(True, alpha=0.3)
        
        plt.yscale('log')
        plt.tight_layout()
        plt.savefig(output_path / 'timing_comparison.png', dpi=300, bbox_inches='tight')
        plt.show()
    
    print(f"Visualizations saved to: {output_path}")

def main():
    """Main function to run the consistency test."""
    print("="*80)
    print("TOPOLOGY METRICS CONSISTENCY TEST")
    print("="*80)
    
    # File paths
    data1_path = Path("migration_topology_weights.csv")
    data2_path = Path("NoMigration_topology_weights.csv")
    
    # Check files exist
    if not data1_path.exists() or not data2_path.exists():
        print("Error: Required CSV files not found!")
        print(f"Looking for: {data1_path} and {data2_path}")
        return
    
    # Create test datasets
    dataset_pairs = create_test_datasets(data1_path, data2_path, n_pairs=20, base_sample_size=1000)
    
    # Run tests
    print(f"\nRunning consistency tests on {len(dataset_pairs)} dataset pairs...")
    results = []
    
    for i, (data1, data2, pair_name) in enumerate(dataset_pairs):
        try:
            result = test_single_pair(data1, data2, pair_name, verbose=(i < 3))  # Verbose for first 3
            results.append(result)
            
            if i >= 3:  # Print progress for remaining tests
                print(f"  Completed {pair_name} ({i+1}/{len(dataset_pairs)})")
                
        except Exception as e:
            print(f"Error testing {pair_name}: {e}")
            continue
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    # Save results
    output_path = Path("consistency_test_results")
    output_path.mkdir(exist_ok=True)
    results_df.to_csv(output_path / "consistency_test_results.csv", index=False)
    
    # Analyze results
    analyze_consistency(results_df)
    
    # Create visualizations
    create_visualization(results_df)
    
    print(f"\n\nFull results saved to: {output_path / 'consistency_test_results.csv'}")
    print("="*80)

if __name__ == "__main__":
    main() 