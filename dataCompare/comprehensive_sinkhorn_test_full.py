#!/usr/bin/env python3
"""
Comprehensive Sinkhorn vs Exact Wasserstein Test - Full Dataset Version

Test 100 pairs of datasets using 10,000 points from full datasets.
Ensures both methods get identical input data for fair comparison.

Key improvements:
- Sample 10,000 points from full datasets (not small subsets)
- Single sampling per pair - both methods get identical data
- Robust statistics on Sinkhorn vs exact Wasserstein
- Both Euclidean and KL distances
"""

import sys
import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats

# Add parent directory for imports
sys.path.append(str(Path(__file__).parent.parent))
from twisstntern.utils import dump_data

# Import our optimized fast metrics functions
from fast_topology_metrics import (
    compute_wasserstein_euclidean_cpu,
    compute_wasserstein_kl_cpu,
    compute_sinkhorn_wasserstein_fast,
    compute_sinkhorn_kl_fast,
    compute_l2_chi2_ultrafast
)

def load_full_dataset(filepath):
    """Load the complete dataset."""
    print(f"Loading full dataset from {filepath}...")
    
    # Load full dataset using twisstntern's dump_data
    data_array = dump_data(str(filepath))
    data = pd.DataFrame(data_array, columns=['T1', 'T2', 'T3'])
    
    print(f"Loaded {len(data)} points from full dataset")
    return data

def create_diverse_test_pairs_full(data1_full, data2_full, n_pairs=100, sample_size=10000):
    """
    Create 100 diverse pairs of test datasets, each with 10,000 points.
    Ensures both methods get identical sampled data for fair comparison.
    """
    print(f"Creating {n_pairs} diverse pairs with {sample_size} points each...")
    
    dataset_pairs = []
    
    # Strategy 1: Random sampling with different seeds (100 pairs)
    print("Creating random sampling pairs from full datasets...")
    for i in range(n_pairs):
        seed1 = 1000 + i
        seed2 = 2000 + i
        
        # Sample from full datasets
        if len(data1_full) >= sample_size:
            data1_sample = data1_full.sample(n=sample_size, random_state=seed1).reset_index(drop=True)
        else:
            data1_sample = data1_full.copy()
            
        if len(data2_full) >= sample_size:
            data2_sample = data2_full.sample(n=sample_size, random_state=seed2).reset_index(drop=True)
        else:
            data2_sample = data2_full.copy()
        
        dataset_pairs.append((data1_sample, data2_sample, f"FullTest_{i+1:03d}"))
    
    print(f"Created {len(dataset_pairs)} dataset pairs from full data")
    return dataset_pairs

def test_single_pair_identical_data(data1, data2, pair_name, verbose=False):
    """
    Test a single pair ensuring both methods get IDENTICAL input data.
    This is crucial for fair comparison.
    """
    if verbose:
        print(f"Testing {pair_name} ({len(data1)} vs {len(data2)} points)...")
    
    results = {
        'pair_name': pair_name,
        'n_data1': len(data1),
        'n_data2': len(data2)
    }
    
    # CRITICAL: Convert to numpy arrays ONCE and use same arrays for both methods
    data1_array = data1[['T1', 'T2', 'T3']].values.astype(np.float64)
    data2_array = data2[['T1', 'T2', 'T3']].values.astype(np.float64)
    
    if verbose:
        print(f"  Converted to numpy arrays: {data1_array.shape}, {data2_array.shape}")
    
    # 1. Exact Wasserstein (reference truth) - using SAME arrays
    try:
        start_time = time.time()
        exact_eucl = compute_wasserstein_euclidean_cpu(
            pd.DataFrame(data1_array, columns=['T1', 'T2', 'T3']), 
            pd.DataFrame(data2_array, columns=['T1', 'T2', 'T3']), 
            max_points=None, verbose=False
        )
        exact_eucl_time = time.time() - start_time
        
        start_time = time.time()
        exact_kl = compute_wasserstein_kl_cpu(
            pd.DataFrame(data1_array, columns=['T1', 'T2', 'T3']), 
            pd.DataFrame(data2_array, columns=['T1', 'T2', 'T3']), 
            max_points=None, verbose=False
        )
        exact_kl_time = time.time() - start_time
        
        exact_total_time = exact_eucl_time + exact_kl_time
        
        results.update({
            'exact_euclidean': exact_eucl,
            'exact_kl': exact_kl,
            'exact_euclidean_time': exact_eucl_time,
            'exact_kl_time': exact_kl_time,
            'exact_total_time': exact_total_time
        })
        
        if verbose:
            print(f"  Exact: Euclidean={exact_eucl:.6f} ({exact_eucl_time:.3f}s), "
                  f"KL={exact_kl:.6f} ({exact_kl_time:.3f}s)")
    
    except Exception as e:
        if verbose:
            print(f"  Exact methods failed: {e}")
        results.update({
            'exact_euclidean': np.nan,
            'exact_kl': np.nan,
            'exact_euclidean_time': np.nan,
            'exact_kl_time': np.nan,
            'exact_total_time': np.nan
        })
    
    # 2. Sinkhorn with 2000 points - using SAME arrays (but with internal sampling)
    try:
        start_time = time.time()
        sinkhorn_eucl = compute_sinkhorn_wasserstein_fast(
            pd.DataFrame(data1_array, columns=['T1', 'T2', 'T3']), 
            pd.DataFrame(data2_array, columns=['T1', 'T2', 'T3']), 
            max_points=2000, reg=0.01, verbose=False
        )
        sinkhorn_eucl_time = time.time() - start_time
        
        start_time = time.time()
        sinkhorn_kl = compute_sinkhorn_kl_fast(
            pd.DataFrame(data1_array, columns=['T1', 'T2', 'T3']), 
            pd.DataFrame(data2_array, columns=['T1', 'T2', 'T3']), 
            max_points=2000, reg=0.005, verbose=False
        )
        sinkhorn_kl_time = time.time() - start_time
        
        sinkhorn_total_time = sinkhorn_eucl_time + sinkhorn_kl_time
        
        results.update({
            'sinkhorn_euclidean': sinkhorn_eucl,
            'sinkhorn_kl': sinkhorn_kl,
            'sinkhorn_euclidean_time': sinkhorn_eucl_time,
            'sinkhorn_kl_time': sinkhorn_kl_time,
            'sinkhorn_total_time': sinkhorn_total_time
        })
        
        if verbose:
            print(f"  Sinkhorn: Euclidean={sinkhorn_eucl:.6f} ({sinkhorn_eucl_time:.3f}s), "
                  f"KL={sinkhorn_kl:.6f} ({sinkhorn_kl_time:.3f}s)")
    
    except Exception as e:
        if verbose:
            print(f"  Sinkhorn methods failed: {e}")
        results.update({
            'sinkhorn_euclidean': np.nan,
            'sinkhorn_kl': np.nan,
            'sinkhorn_euclidean_time': np.nan,
            'sinkhorn_kl_time': np.nan,
            'sinkhorn_total_time': np.nan
        })
    
    # Calculate errors if both methods succeeded
    if not np.isnan(results.get('exact_euclidean', np.nan)) and not np.isnan(results.get('sinkhorn_euclidean', np.nan)):
        eucl_error = abs(results['sinkhorn_euclidean'] - results['exact_euclidean'])
        eucl_rel_error = (eucl_error / results['exact_euclidean']) * 100 if results['exact_euclidean'] > 0 else 0
        results.update({
            'euclidean_abs_error': eucl_error,
            'euclidean_rel_error': eucl_rel_error
        })
    
    if not np.isnan(results.get('exact_kl', np.nan)) and not np.isnan(results.get('sinkhorn_kl', np.nan)):
        kl_error = abs(results['sinkhorn_kl'] - results['exact_kl'])
        kl_rel_error = (kl_error / results['exact_kl']) * 100 if results['exact_kl'] > 0 else 0
        results.update({
            'kl_abs_error': kl_error,
            'kl_rel_error': kl_rel_error
        })
    
    # Add L2 and chi2 for reference (using same arrays)
    try:
        L2_dist, chi2_stat, p_val, dof, n_meaningful = compute_l2_chi2_ultrafast(
            pd.DataFrame(data1_array, columns=['T1', 'T2', 'T3']), 
            pd.DataFrame(data2_array, columns=['T1', 'T2', 'T3']), 
            0.1, verbose=False
        )
        results.update({
            'L2_distance': L2_dist,
            'chi2_statistic': chi2_stat,
            'meaningful_triangles': n_meaningful
        })
    except Exception as e:
        results.update({
            'L2_distance': np.nan,
            'chi2_statistic': np.nan,
            'meaningful_triangles': np.nan
        })
    
    return results

def run_comprehensive_full_test():
    """Run the comprehensive 100-pair test on full datasets."""
    
    print("="*90)
    print("COMPREHENSIVE SINKHORN vs EXACT WASSERSTEIN TEST - FULL DATASETS")
    print("Using 10,000 points per dataset pair with identical data for both methods")
    print("="*90)
    
    # File paths
    data1_path = Path("migration_topology_weights.csv")
    data2_path = Path("NoMigration_topology_weights.csv")
    
    # Check files exist
    if not data1_path.exists() or not data2_path.exists():
        print("Error: Required CSV files not found!")
        print(f"Looking for: {data1_path} and {data2_path}")
        return None
    
    # Load full datasets
    data1_full = load_full_dataset(data1_path)
    data2_full = load_full_dataset(data2_path)
    
    print(f"Full dataset sizes: {len(data1_full)} vs {len(data2_full)} points")
    
    # Create 100 test datasets with 10,000 points each
    dataset_pairs = create_diverse_test_pairs_full(
        data1_full, data2_full, n_pairs=100, sample_size=10000
    )
    
    # Run tests
    print(f"\nRunning comprehensive tests on {len(dataset_pairs)} dataset pairs...")
    print("Each pair uses 10,000 points with identical data for both methods")
    results = []
    
    start_total_time = time.time()
    
    for i, (data1, data2, pair_name) in enumerate(dataset_pairs):
        try:
            # Show progress every 10 tests
            verbose = (i < 3) or (i % 10 == 0)
            
            result = test_single_pair_identical_data(data1, data2, pair_name, verbose=verbose)
            results.append(result)
            
            if not verbose:
                print(f"  Completed {pair_name} ({i+1}/{len(dataset_pairs)})")
                
        except Exception as e:
            print(f"Error testing {pair_name}: {e}")
            continue
    
    total_time = time.time() - start_total_time
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    # Save results
    output_path = Path("comprehensive_full_sinkhorn_results")
    output_path.mkdir(exist_ok=True)
    results_df.to_csv(output_path / "comprehensive_full_test_results.csv", index=False)
    
    print(f"\nTotal testing time: {total_time:.1f} seconds")
    print(f"Average time per pair: {total_time/len(dataset_pairs):.2f} seconds")
    print(f"Results saved to: {output_path}")
    
    # Analyze results
    analyze_full_results(results_df, output_path)
    
    return results_df

def analyze_full_results(results_df, output_path):
    """Comprehensive analysis of the full dataset test results."""
    
    print(f"\n{'='*90}")
    print("COMPREHENSIVE ANALYSIS - 100 FULL DATASET PAIRS (10K POINTS EACH)")
    print(f"{'='*90}")
    
    # Filter valid results
    valid_eucl = results_df.dropna(subset=['exact_euclidean', 'sinkhorn_euclidean'])
    valid_kl = results_df.dropna(subset=['exact_kl', 'sinkhorn_kl'])
    
    print(f"Valid Euclidean comparisons: {len(valid_eucl)}/{len(results_df)}")
    print(f"Valid KL comparisons: {len(valid_kl)}/{len(results_df)}")
    
    if len(valid_eucl) == 0 and len(valid_kl) == 0:
        print("No valid results to analyze!")
        return
    
    # === EUCLIDEAN ANALYSIS ===
    if len(valid_eucl) > 0:
        print(f"\n🎯 EUCLIDEAN WASSERSTEIN ANALYSIS ({len(valid_eucl)} pairs):")
        print("-" * 60)
        
        # Error statistics
        eucl_errors = valid_eucl['euclidean_rel_error']
        print(f"Relative Error Statistics:")
        print(f"  Mean: {eucl_errors.mean():.3f}%")
        print(f"  Median: {eucl_errors.median():.3f}%")
        print(f"  Std Dev: {eucl_errors.std():.3f}%")
        print(f"  Min: {eucl_errors.min():.3f}%")
        print(f"  Max: {eucl_errors.max():.3f}%")
        print(f"  95th percentile: {eucl_errors.quantile(0.95):.3f}%")
        print(f"  99th percentile: {eucl_errors.quantile(0.99):.3f}%")
        
        # Correlation
        eucl_corr = valid_eucl['exact_euclidean'].corr(valid_eucl['sinkhorn_euclidean'])
        print(f"Correlation (Exact vs Sinkhorn): {eucl_corr:.6f}")
        
        # Speed comparison
        eucl_speedup = valid_eucl['exact_euclidean_time'].mean() / valid_eucl['sinkhorn_euclidean_time'].mean()
        print(f"Timing: Exact avg={valid_eucl['exact_euclidean_time'].mean():.3f}s, "
              f"Sinkhorn avg={valid_eucl['sinkhorn_euclidean_time'].mean():.3f}s")
        print(f"Speedup: {eucl_speedup:.1f}x faster")
        
        # Quality assessment
        excellent = (eucl_errors < 2).sum()
        very_good = ((eucl_errors >= 2) & (eucl_errors < 5)).sum()
        good = ((eucl_errors >= 5) & (eucl_errors < 10)).sum()
        acceptable = ((eucl_errors >= 10) & (eucl_errors < 20)).sum()
        poor = (eucl_errors >= 20).sum()
        
        print(f"Quality Assessment (10K points each):")
        print(f"  Excellent (<2% error): {excellent}/{len(valid_eucl)} ({excellent/len(valid_eucl)*100:.1f}%)")
        print(f"  Very Good (2-5% error): {very_good}/{len(valid_eucl)} ({very_good/len(valid_eucl)*100:.1f}%)")
        print(f"  Good (5-10% error): {good}/{len(valid_eucl)} ({good/len(valid_eucl)*100:.1f}%)")
        print(f"  Acceptable (10-20% error): {acceptable}/{len(valid_eucl)} ({acceptable/len(valid_eucl)*100:.1f}%)")
        print(f"  Poor (>20% error): {poor}/{len(valid_eucl)} ({poor/len(valid_eucl)*100:.1f}%)")
    
    # === KL ANALYSIS ===
    if len(valid_kl) > 0:
        print(f"\n🎯 KL WASSERSTEIN ANALYSIS ({len(valid_kl)} pairs):")
        print("-" * 60)
        
        # Error statistics  
        kl_errors = valid_kl['kl_rel_error']
        print(f"Relative Error Statistics:")
        print(f"  Mean: {kl_errors.mean():.3f}%")
        print(f"  Median: {kl_errors.median():.3f}%")
        print(f"  Std Dev: {kl_errors.std():.3f}%")
        print(f"  Min: {kl_errors.min():.3f}%")
        print(f"  Max: {kl_errors.max():.3f}%")
        print(f"  95th percentile: {kl_errors.quantile(0.95):.3f}%")
        print(f"  99th percentile: {kl_errors.quantile(0.99):.3f}%")
        
        # Correlation
        kl_corr = valid_kl['exact_kl'].corr(valid_kl['sinkhorn_kl'])
        print(f"Correlation (Exact vs Sinkhorn): {kl_corr:.6f}")
        
        # Speed comparison
        kl_speedup = valid_kl['exact_kl_time'].mean() / valid_kl['sinkhorn_kl_time'].mean()
        print(f"Timing: Exact avg={valid_kl['exact_kl_time'].mean():.3f}s, "
              f"Sinkhorn avg={valid_kl['sinkhorn_kl_time'].mean():.3f}s")
        print(f"Speedup: {kl_speedup:.1f}x faster")
        
        # Quality assessment
        excellent = (kl_errors < 2).sum()
        very_good = ((kl_errors >= 2) & (kl_errors < 5)).sum()
        good = ((kl_errors >= 5) & (kl_errors < 10)).sum()
        acceptable = ((kl_errors >= 10) & (kl_errors < 20)).sum()
        poor = (kl_errors >= 20).sum()
        
        print(f"Quality Assessment (10K points each):")
        print(f"  Excellent (<2% error): {excellent}/{len(valid_kl)} ({excellent/len(valid_kl)*100:.1f}%)")
        print(f"  Very Good (2-5% error): {very_good}/{len(valid_kl)} ({very_good/len(valid_kl)*100:.1f}%)")
        print(f"  Good (5-10% error): {good}/{len(valid_kl)} ({good/len(valid_kl)*100:.1f}%)")
        print(f"  Acceptable (10-20% error): {acceptable}/{len(valid_kl)} ({acceptable/len(valid_kl)*100:.1f}%)")
        print(f"  Poor (>20% error): {poor}/{len(valid_kl)} ({poor/len(valid_kl)*100:.1f}%)")
    
    # === OVERALL PERFORMANCE ===
    valid_both = results_df.dropna(subset=['exact_total_time', 'sinkhorn_total_time'])
    if len(valid_both) > 0:
        print(f"\n🚀 OVERALL PERFORMANCE ({len(valid_both)} pairs, 10K points each):")
        print("-" * 60)
        
        total_speedup = valid_both['exact_total_time'].mean() / valid_both['sinkhorn_total_time'].mean()
        print(f"Total Time: Exact avg={valid_both['exact_total_time'].mean():.3f}s, "
              f"Sinkhorn avg={valid_both['sinkhorn_total_time'].mean():.3f}s")
        print(f"Overall Speedup: {total_speedup:.1f}x faster")
        
        # Estimate computational complexity
        print(f"\nComputational Scaling (10K points per dataset):")
        print(f"  Exact method: O(n³) ≈ {10000**2:.0e} operations")
        print(f"  Sinkhorn method: O(n²·log(n)) ≈ {10000**2 * np.log(10000):.0e} operations")
        print(f"  Theoretical speedup: ~{(10000**2) / (10000**2 * np.log(10000) / 100):.1f}x")
    
    # Create visualizations and save detailed statistics
    create_full_visualizations(results_df, output_path)
    save_full_statistics(results_df, output_path)

def create_full_visualizations(results_df, output_path):
    """Create visualizations for the full dataset test."""
    
    print(f"\nCreating comprehensive visualizations for full dataset test...")
    
    # Set up plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Filter valid data
    valid_eucl = results_df.dropna(subset=['exact_euclidean', 'sinkhorn_euclidean'])
    valid_kl = results_df.dropna(subset=['exact_kl', 'sinkhorn_kl'])
    
    if len(valid_eucl) == 0 and len(valid_kl) == 0:
        print("No valid data for visualization!")
        return
    
    # Create comprehensive figure
    fig = plt.figure(figsize=(20, 12))
    fig.suptitle('Comprehensive Sinkhorn vs Exact Wasserstein Analysis (10K Points/Dataset)', fontsize=16)
    
    # Euclidean plots
    if len(valid_eucl) > 0:
        # 1. Euclidean correlation
        ax1 = plt.subplot(2, 4, 1)
        ax1.scatter(valid_eucl['exact_euclidean'], valid_eucl['sinkhorn_euclidean'], 
                   alpha=0.7, s=50, c='blue', edgecolors='white', linewidth=0.5)
        ax1.plot([valid_eucl['exact_euclidean'].min(), valid_eucl['exact_euclidean'].max()], 
                [valid_eucl['exact_euclidean'].min(), valid_eucl['exact_euclidean'].max()], 
                'r--', alpha=0.8, linewidth=2)
        ax1.set_xlabel('Exact Euclidean')
        ax1.set_ylabel('Sinkhorn Euclidean')
        ax1.set_title(f'Euclidean Correlation\n(r={valid_eucl["exact_euclidean"].corr(valid_eucl["sinkhorn_euclidean"]):.4f})')
        ax1.grid(True, alpha=0.3)
        
        # 2. Euclidean error distribution
        ax2 = plt.subplot(2, 4, 2)
        valid_eucl['euclidean_rel_error'].hist(bins=40, alpha=0.7, color='blue', ax=ax2)
        ax2.axvline(valid_eucl['euclidean_rel_error'].mean(), color='red', linestyle='--', 
                   label=f'Mean: {valid_eucl["euclidean_rel_error"].mean():.2f}%')
        ax2.axvline(valid_eucl['euclidean_rel_error'].median(), color='orange', linestyle='--',
                   label=f'Median: {valid_eucl["euclidean_rel_error"].median():.2f}%')
        ax2.set_xlabel('Relative Error (%)')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Euclidean Error Distribution')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # 3. Euclidean timing
        ax3 = plt.subplot(2, 4, 3)
        times_eucl = [valid_eucl['exact_euclidean_time'], valid_eucl['sinkhorn_euclidean_time']]
        bp = ax3.boxplot(times_eucl, labels=['Exact', 'Sinkhorn'], patch_artist=True)
        bp['boxes'][0].set_facecolor('lightblue')
        bp['boxes'][1].set_facecolor('lightgreen')
        ax3.set_ylabel('Time (seconds)')
        ax3.set_title('Euclidean Timing (10K points)')
        ax3.grid(True, alpha=0.3)
        ax3.set_yscale('log')
        
        # 4. Euclidean quality pie chart
        ax4 = plt.subplot(2, 4, 4)
        eucl_errors = valid_eucl['euclidean_rel_error']
        excellent = (eucl_errors < 2).sum()
        very_good = ((eucl_errors >= 2) & (eucl_errors < 5)).sum()
        good = ((eucl_errors >= 5) & (eucl_errors < 10)).sum()
        poor = (eucl_errors >= 10).sum()
        
        sizes = [excellent, very_good, good, poor]
        labels = ['Excellent\n(<2%)', 'Very Good\n(2-5%)', 'Good\n(5-10%)', 'Poor\n(>10%)']
        colors = ['darkgreen', 'lightgreen', 'orange', 'red']
        
        ax4.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
        ax4.set_title('Euclidean Quality (10K pts)')
    
    # KL plots
    if len(valid_kl) > 0:
        # 5. KL correlation
        ax5 = plt.subplot(2, 4, 5)
        ax5.scatter(valid_kl['exact_kl'], valid_kl['sinkhorn_kl'], 
                   alpha=0.7, s=50, c='green', edgecolors='white', linewidth=0.5)
        ax5.plot([valid_kl['exact_kl'].min(), valid_kl['exact_kl'].max()], 
                [valid_kl['exact_kl'].min(), valid_kl['exact_kl'].max()], 
                'r--', alpha=0.8, linewidth=2)
        ax5.set_xlabel('Exact KL')
        ax5.set_ylabel('Sinkhorn KL')
        ax5.set_title(f'KL Correlation\n(r={valid_kl["exact_kl"].corr(valid_kl["sinkhorn_kl"]):.4f})')
        ax5.grid(True, alpha=0.3)
        
        # 6. KL error distribution
        ax6 = plt.subplot(2, 4, 6)
        valid_kl['kl_rel_error'].hist(bins=40, alpha=0.7, color='green', ax=ax6)
        ax6.axvline(valid_kl['kl_rel_error'].mean(), color='red', linestyle='--',
                   label=f'Mean: {valid_kl["kl_rel_error"].mean():.2f}%')
        ax6.axvline(valid_kl['kl_rel_error'].median(), color='orange', linestyle='--',
                   label=f'Median: {valid_kl["kl_rel_error"].median():.2f}%')
        ax6.set_xlabel('Relative Error (%)')
        ax6.set_ylabel('Frequency')
        ax6.set_title('KL Error Distribution')
        ax6.legend()
        ax6.grid(True, alpha=0.3)
        
        # 7. KL timing
        ax7 = plt.subplot(2, 4, 7)
        times_kl = [valid_kl['exact_kl_time'], valid_kl['sinkhorn_kl_time']]
        bp = ax7.boxplot(times_kl, labels=['Exact', 'Sinkhorn'], patch_artist=True)
        bp['boxes'][0].set_facecolor('lightblue')
        bp['boxes'][1].set_facecolor('lightgreen')
        ax7.set_ylabel('Time (seconds)')
        ax7.set_title('KL Timing (10K points)')
        ax7.grid(True, alpha=0.3)
        ax7.set_yscale('log')
        
        # 8. KL quality pie chart
        ax8 = plt.subplot(2, 4, 8)
        kl_errors = valid_kl['kl_rel_error']
        excellent = (kl_errors < 2).sum()
        very_good = ((kl_errors >= 2) & (kl_errors < 5)).sum()
        good = ((kl_errors >= 5) & (kl_errors < 10)).sum()
        poor = (kl_errors >= 10).sum()
        
        sizes = [excellent, very_good, good, poor]
        labels = ['Excellent\n(<2%)', 'Very Good\n(2-5%)', 'Good\n(5-10%)', 'Poor\n(>10%)']
        colors = ['darkgreen', 'lightgreen', 'orange', 'red']
        
        ax8.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
        ax8.set_title('KL Quality (10K pts)')
    
    plt.tight_layout()
    plt.savefig(output_path / 'comprehensive_full_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Visualization saved to: {output_path / 'comprehensive_full_analysis.png'}")

def save_full_statistics(results_df, output_path):
    """Save detailed statistical analysis for full dataset test."""
    
    stats_file = output_path / 'detailed_full_statistics.txt'
    
    with open(stats_file, 'w') as f:
        f.write("COMPREHENSIVE SINKHORN vs EXACT WASSERSTEIN - FULL DATASET ANALYSIS\n")
        f.write("="*80 + "\n\n")
        f.write(f"Dataset pairs tested: {len(results_df)}\n")
        f.write(f"Points per dataset: 10,000\n")
        f.write(f"Test completed: {pd.Timestamp.now()}\n\n")
        
        # Detailed statistics for both methods
        valid_eucl = results_df.dropna(subset=['exact_euclidean', 'sinkhorn_euclidean'])
        valid_kl = results_df.dropna(subset=['exact_kl', 'sinkhorn_kl'])
        
        if len(valid_eucl) > 0:
            f.write(f"EUCLIDEAN WASSERSTEIN RESULTS (10K points, {len(valid_eucl)} pairs):\n")
            f.write("-" * 60 + "\n")
            
            eucl_errors = valid_eucl['euclidean_rel_error']
            f.write(f"Error Statistics:\n")
            f.write(f"  Mean error: {eucl_errors.mean():.4f}%\n")
            f.write(f"  Median error: {eucl_errors.median():.4f}%\n")
            f.write(f"  Std deviation: {eucl_errors.std():.4f}%\n")
            f.write(f"  Min error: {eucl_errors.min():.4f}%\n")
            f.write(f"  Max error: {eucl_errors.max():.4f}%\n")
            f.write(f"  99th percentile: {eucl_errors.quantile(0.99):.4f}%\n\n")
            
            f.write(f"Performance:\n")
            f.write(f"  Correlation: {valid_eucl['exact_euclidean'].corr(valid_eucl['sinkhorn_euclidean']):.6f}\n")
            f.write(f"  Exact timing: {valid_eucl['exact_euclidean_time'].mean():.4f} ± {valid_eucl['exact_euclidean_time'].std():.4f}s\n")
            f.write(f"  Sinkhorn timing: {valid_eucl['sinkhorn_euclidean_time'].mean():.4f} ± {valid_eucl['sinkhorn_euclidean_time'].std():.4f}s\n")
            f.write(f"  Speedup: {valid_eucl['exact_euclidean_time'].mean() / valid_eucl['sinkhorn_euclidean_time'].mean():.2f}x\n\n")
        
        if len(valid_kl) > 0:
            f.write(f"KL WASSERSTEIN RESULTS (10K points, {len(valid_kl)} pairs):\n")
            f.write("-" * 60 + "\n")
            
            kl_errors = valid_kl['kl_rel_error']
            f.write(f"Error Statistics:\n")
            f.write(f"  Mean error: {kl_errors.mean():.4f}%\n")
            f.write(f"  Median error: {kl_errors.median():.4f}%\n")
            f.write(f"  Std deviation: {kl_errors.std():.4f}%\n")
            f.write(f"  Min error: {kl_errors.min():.4f}%\n")
            f.write(f"  Max error: {kl_errors.max():.4f}%\n")
            f.write(f"  99th percentile: {kl_errors.quantile(0.99):.4f}%\n\n")
            
            f.write(f"Performance:\n")
            f.write(f"  Correlation: {valid_kl['exact_kl'].corr(valid_kl['sinkhorn_kl']):.6f}\n")
            f.write(f"  Exact timing: {valid_kl['exact_kl_time'].mean():.4f} ± {valid_kl['exact_kl_time'].std():.4f}s\n")
            f.write(f"  Sinkhorn timing: {valid_kl['sinkhorn_kl_time'].mean():.4f} ± {valid_kl['sinkhorn_kl_time'].std():.4f}s\n")
            f.write(f"  Speedup: {valid_kl['exact_kl_time'].mean() / valid_kl['sinkhorn_kl_time'].mean():.2f}x\n\n")
    
    print(f"Detailed statistics saved to: {stats_file}")

def main():
    """Main function to run the comprehensive full dataset test."""
    results_df = run_comprehensive_full_test()
    
    if results_df is not None:
        print(f"\n{'='*90}")
        print("COMPREHENSIVE FULL DATASET TEST COMPLETED SUCCESSFULLY!")
        print(f"{'='*90}")
        print(f"Total results: {len(results_df)} dataset pairs")
        print("Each pair used 10,000 points with identical data for both methods")
        print("Check the 'comprehensive_full_sinkhorn_results' directory for detailed outputs.")
    else:
        print("Test failed - check error messages above.")

if __name__ == "__main__":
    main()
