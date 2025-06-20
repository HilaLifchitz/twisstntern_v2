#!/usr/bin/env python3
"""
Comprehensive Sinkhorn vs Exact Wasserstein Test

Test 100 pairs of datasets to get robust statistics on:
- Sinkhorn 2000 (optimized) vs Exact Wasserstein
- Both Euclidean and KL distances
- Statistical analysis of bias and consistency
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

def create_diverse_test_datasets(data1_path, data2_path, n_pairs=100, base_sample_size=1000):
    """Create 100 diverse pairs of test datasets with different sampling strategies."""
    print(f"Creating {n_pairs} diverse pairs of test datasets...")
    
    # Load the full datasets
    print("Loading base datasets...")
    full_data1 = load_and_sample_data(data1_path, sample_size=base_sample_size*n_pairs//2, random_state=42)
    full_data2 = load_and_sample_data(data2_path, sample_size=base_sample_size*n_pairs//2, random_state=43)
    
    dataset_pairs = []
    
    # Strategy 1: Sequential sampling (50 pairs)
    print("Creating sequential sampling pairs...")
    for i in range(50):
        start_idx = i * (base_sample_size // 2)
        end_idx = (i + 1) * (base_sample_size // 2)
        
        if end_idx <= len(full_data1) and end_idx <= len(full_data2):
            data1_sample = full_data1.iloc[start_idx:end_idx].copy().reset_index(drop=True)
            data2_sample = full_data2.iloc[start_idx:end_idx].copy().reset_index(drop=True)
            
            dataset_pairs.append((data1_sample, data2_sample, f"Sequential_{i+1:02d}"))
    
    # Strategy 2: Random sampling with different seeds (50 pairs)
    print("Creating random sampling pairs...")
    for i in range(50):
        seed1 = 100 + i
        seed2 = 200 + i
        
        data1_sample = full_data1.sample(n=base_sample_size//2, random_state=seed1).reset_index(drop=True)
        data2_sample = full_data2.sample(n=base_sample_size//2, random_state=seed2).reset_index(drop=True)
        
        dataset_pairs.append((data1_sample, data2_sample, f"Random_{i+1:02d}"))
    
    print(f"Created {len(dataset_pairs)} dataset pairs")
    return dataset_pairs

def test_single_pair_comprehensive(data1, data2, pair_name, verbose=False):
    """Test a single pair with both exact and Sinkhorn methods."""
    if verbose:
        print(f"Testing {pair_name} ({len(data1)} vs {len(data2)} points)...")
    
    results = {
        'pair_name': pair_name,
        'n_data1': len(data1),
        'n_data2': len(data2)
    }
    
    # 1. Exact Wasserstein (reference truth)
    try:
        start_time = time.time()
        exact_eucl = compute_wasserstein_euclidean_cpu(data1, data2, max_points=None, verbose=False)
        exact_eucl_time = time.time() - start_time
        
        start_time = time.time()
        exact_kl = compute_wasserstein_kl_cpu(data1, data2, max_points=None, verbose=False)
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
    
    # 2. Optimized Sinkhorn with 2000 points
    try:
        start_time = time.time()
        sinkhorn_eucl = compute_sinkhorn_wasserstein_fast(data1, data2, max_points=2000, reg=0.01, verbose=False)
        sinkhorn_eucl_time = time.time() - start_time
        
        start_time = time.time()
        sinkhorn_kl = compute_sinkhorn_kl_fast(data1, data2, max_points=2000, reg=0.005, verbose=False)
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
        eucl_rel_error = (eucl_error / results['exact_euclidean']) * 100
        results.update({
            'euclidean_abs_error': eucl_error,
            'euclidean_rel_error': eucl_rel_error
        })
    
    if not np.isnan(results.get('exact_kl', np.nan)) and not np.isnan(results.get('sinkhorn_kl', np.nan)):
        kl_error = abs(results['sinkhorn_kl'] - results['exact_kl'])
        kl_rel_error = (kl_error / results['exact_kl']) * 100
        results.update({
            'kl_abs_error': kl_error,
            'kl_rel_error': kl_rel_error
        })
    
    # Add L2 and chi2 for reference
    try:
        L2_dist, chi2_stat, p_val, dof, n_meaningful = compute_l2_chi2_ultrafast(data1, data2, 0.1, verbose=False)
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

def run_comprehensive_test():
    """Run the comprehensive 100-pair test."""
    
    print("="*90)
    print("COMPREHENSIVE SINKHORN vs EXACT WASSERSTEIN TEST (100 PAIRS)")
    print("="*90)
    
    # File paths
    data1_path = Path("migration_topology_weights.csv")
    data2_path = Path("NoMigration_topology_weights.csv")
    
    # Check files exist
    if not data1_path.exists() or not data2_path.exists():
        print("Error: Required CSV files not found!")
        print(f"Looking for: {data1_path} and {data2_path}")
        return None
    
    # Create 100 test datasets
    dataset_pairs = create_diverse_test_datasets(data1_path, data2_path, n_pairs=100, base_sample_size=1000)
    
    # Run tests
    print(f"\nRunning comprehensive tests on {len(dataset_pairs)} dataset pairs...")
    results = []
    
    start_total_time = time.time()
    
    for i, (data1, data2, pair_name) in enumerate(dataset_pairs):
        try:
            # Show progress every 10 tests
            verbose = (i < 3) or (i % 10 == 0)
            
            result = test_single_pair_comprehensive(data1, data2, pair_name, verbose=verbose)
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
    output_path = Path("comprehensive_sinkhorn_results")
    output_path.mkdir(exist_ok=True)
    results_df.to_csv(output_path / "comprehensive_test_results.csv", index=False)
    
    print(f"\nTotal testing time: {total_time:.1f} seconds")
    print(f"Results saved to: {output_path}")
    
    # Analyze results
    analyze_comprehensive_results(results_df, output_path)
    
    return results_df

def analyze_comprehensive_results(results_df, output_path):
    """Comprehensive analysis of the 100-pair test results."""
    
    print(f"\n{'='*90}")
    print("COMPREHENSIVE ANALYSIS - 100 DATASET PAIRS")
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
        print(f"  Mean: {eucl_errors.mean():.2f}%")
        print(f"  Median: {eucl_errors.median():.2f}%")
        print(f"  Std Dev: {eucl_errors.std():.2f}%")
        print(f"  Min: {eucl_errors.min():.2f}%")
        print(f"  Max: {eucl_errors.max():.2f}%")
        print(f"  95th percentile: {eucl_errors.quantile(0.95):.2f}%")
        
        # Correlation
        eucl_corr = valid_eucl['exact_euclidean'].corr(valid_eucl['sinkhorn_euclidean'])
        print(f"Correlation (Exact vs Sinkhorn): {eucl_corr:.6f}")
        
        # Speed comparison
        eucl_speedup = valid_eucl['exact_euclidean_time'].mean() / valid_eucl['sinkhorn_euclidean_time'].mean()
        print(f"Speed: Exact avg={valid_eucl['exact_euclidean_time'].mean():.3f}s, "
              f"Sinkhorn avg={valid_eucl['sinkhorn_euclidean_time'].mean():.3f}s")
        print(f"Speedup: {eucl_speedup:.1f}x faster")
        
        # Quality assessment
        excellent = (eucl_errors < 5).sum()
        good = ((eucl_errors >= 5) & (eucl_errors < 10)).sum()
        acceptable = ((eucl_errors >= 10) & (eucl_errors < 20)).sum()
        poor = (eucl_errors >= 20).sum()
        
        print(f"Quality Assessment:")
        print(f"  Excellent (<5% error): {excellent}/{len(valid_eucl)} ({excellent/len(valid_eucl)*100:.1f}%)")
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
        print(f"  Mean: {kl_errors.mean():.2f}%")
        print(f"  Median: {kl_errors.median():.2f}%")
        print(f"  Std Dev: {kl_errors.std():.2f}%")
        print(f"  Min: {kl_errors.min():.2f}%")
        print(f"  Max: {kl_errors.max():.2f}%")
        print(f"  95th percentile: {kl_errors.quantile(0.95):.2f}%")
        
        # Correlation
        kl_corr = valid_kl['exact_kl'].corr(valid_kl['sinkhorn_kl'])
        print(f"Correlation (Exact vs Sinkhorn): {kl_corr:.6f}")
        
        # Speed comparison
        kl_speedup = valid_kl['exact_kl_time'].mean() / valid_kl['sinkhorn_kl_time'].mean()
        print(f"Speed: Exact avg={valid_kl['exact_kl_time'].mean():.3f}s, "
              f"Sinkhorn avg={valid_kl['sinkhorn_kl_time'].mean():.3f}s")
        print(f"Speedup: {kl_speedup:.1f}x faster")
        
        # Quality assessment
        excellent = (kl_errors < 5).sum()
        good = ((kl_errors >= 5) & (kl_errors < 10)).sum()
        acceptable = ((kl_errors >= 10) & (kl_errors < 20)).sum()
        poor = (kl_errors >= 20).sum()
        
        print(f"Quality Assessment:")
        print(f"  Excellent (<5% error): {excellent}/{len(valid_kl)} ({excellent/len(valid_kl)*100:.1f}%)")
        print(f"  Good (5-10% error): {good}/{len(valid_kl)} ({good/len(valid_kl)*100:.1f}%)")
        print(f"  Acceptable (10-20% error): {acceptable}/{len(valid_kl)} ({acceptable/len(valid_kl)*100:.1f}%)")
        print(f"  Poor (>20% error): {poor}/{len(valid_kl)} ({poor/len(valid_kl)*100:.1f}%)")
    
    # === OVERALL PERFORMANCE ===
    valid_both = results_df.dropna(subset=['exact_total_time', 'sinkhorn_total_time'])
    if len(valid_both) > 0:
        print(f"\n🚀 OVERALL PERFORMANCE ({len(valid_both)} pairs):")
        print("-" * 60)
        
        total_speedup = valid_both['exact_total_time'].mean() / valid_both['sinkhorn_total_time'].mean()
        print(f"Total Time: Exact avg={valid_both['exact_total_time'].mean():.3f}s, "
              f"Sinkhorn avg={valid_both['sinkhorn_total_time'].mean():.3f}s")
        print(f"Overall Speedup: {total_speedup:.1f}x faster")
    
    # Create comprehensive visualizations
    create_comprehensive_visualizations(results_df, output_path)
    
    # Save detailed statistics
    save_detailed_statistics(results_df, output_path)

def create_comprehensive_visualizations(results_df, output_path):
    """Create comprehensive visualizations for the 100-pair test."""
    
    print(f"\nCreating comprehensive visualizations...")
    
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
    fig = plt.figure(figsize=(20, 16))
    
    # === EUCLIDEAN PLOTS ===
    if len(valid_eucl) > 0:
        # 1. Euclidean: Exact vs Sinkhorn scatter
        ax1 = plt.subplot(3, 4, 1)
        ax1.scatter(valid_eucl['exact_euclidean'], valid_eucl['sinkhorn_euclidean'], 
                   alpha=0.6, s=30, c='blue')
        ax1.plot([valid_eucl['exact_euclidean'].min(), valid_eucl['exact_euclidean'].max()], 
                [valid_eucl['exact_euclidean'].min(), valid_eucl['exact_euclidean'].max()], 
                'r--', alpha=0.8, linewidth=2)
        ax1.set_xlabel('Exact Euclidean')
        ax1.set_ylabel('Sinkhorn Euclidean')
        ax1.set_title(f'Euclidean: Exact vs Sinkhorn\n(r={valid_eucl["exact_euclidean"].corr(valid_eucl["sinkhorn_euclidean"]):.4f})')
        ax1.grid(True, alpha=0.3)
        
        # 2. Euclidean error distribution
        ax2 = plt.subplot(3, 4, 2)
        valid_eucl['euclidean_rel_error'].hist(bins=30, alpha=0.7, color='blue', ax=ax2)
        ax2.axvline(valid_eucl['euclidean_rel_error'].mean(), color='red', linestyle='--', 
                   label=f'Mean: {valid_eucl["euclidean_rel_error"].mean():.1f}%')
        ax2.axvline(valid_eucl['euclidean_rel_error'].median(), color='orange', linestyle='--',
                   label=f'Median: {valid_eucl["euclidean_rel_error"].median():.1f}%')
        ax2.set_xlabel('Relative Error (%)')
        ax2.set_ylabel('Frequency')
        ax2.set_title('Euclidean Error Distribution')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # 3. Euclidean timing comparison
        ax3 = plt.subplot(3, 4, 3)
        times_eucl = [valid_eucl['exact_euclidean_time'], valid_eucl['sinkhorn_euclidean_time']]
        ax3.boxplot(times_eucl, labels=['Exact', 'Sinkhorn'])
        ax3.set_ylabel('Time (seconds)')
        ax3.set_title('Euclidean Timing Comparison')
        ax3.grid(True, alpha=0.3)
        ax3.set_yscale('log')
        
        # 4. Euclidean error vs magnitude
        ax4 = plt.subplot(3, 4, 4)
        ax4.scatter(valid_eucl['exact_euclidean'], valid_eucl['euclidean_rel_error'], 
                   alpha=0.6, s=30, c='blue')
        ax4.set_xlabel('Exact Euclidean Value')
        ax4.set_ylabel('Relative Error (%)')
        ax4.set_title('Euclidean Error vs Magnitude')
        ax4.grid(True, alpha=0.3)
    
    # === KL PLOTS ===
    if len(valid_kl) > 0:
        # 5. KL: Exact vs Sinkhorn scatter
        ax5 = plt.subplot(3, 4, 5)
        ax5.scatter(valid_kl['exact_kl'], valid_kl['sinkhorn_kl'], 
                   alpha=0.6, s=30, c='green')
        ax5.plot([valid_kl['exact_kl'].min(), valid_kl['exact_kl'].max()], 
                [valid_kl['exact_kl'].min(), valid_kl['exact_kl'].max()], 
                'r--', alpha=0.8, linewidth=2)
        ax5.set_xlabel('Exact KL')
        ax5.set_ylabel('Sinkhorn KL')
        ax5.set_title(f'KL: Exact vs Sinkhorn\n(r={valid_kl["exact_kl"].corr(valid_kl["sinkhorn_kl"]):.4f})')
        ax5.grid(True, alpha=0.3)
        
        # 6. KL error distribution
        ax6 = plt.subplot(3, 4, 6)
        valid_kl['kl_rel_error'].hist(bins=30, alpha=0.7, color='green', ax=ax6)
        ax6.axvline(valid_kl['kl_rel_error'].mean(), color='red', linestyle='--',
                   label=f'Mean: {valid_kl["kl_rel_error"].mean():.1f}%')
        ax6.axvline(valid_kl['kl_rel_error'].median(), color='orange', linestyle='--',
                   label=f'Median: {valid_kl["kl_rel_error"].median():.1f}%')
        ax6.set_xlabel('Relative Error (%)')
        ax6.set_ylabel('Frequency')
        ax6.set_title('KL Error Distribution')
        ax6.legend()
        ax6.grid(True, alpha=0.3)
        
        # 7. KL timing comparison
        ax7 = plt.subplot(3, 4, 7)
        times_kl = [valid_kl['exact_kl_time'], valid_kl['sinkhorn_kl_time']]
        ax7.boxplot(times_kl, labels=['Exact', 'Sinkhorn'])
        ax7.set_ylabel('Time (seconds)')
        ax7.set_title('KL Timing Comparison')
        ax7.grid(True, alpha=0.3)
        ax7.set_yscale('log')
        
        # 8. KL error vs magnitude
        ax8 = plt.subplot(3, 4, 8)
        ax8.scatter(valid_kl['exact_kl'], valid_kl['kl_rel_error'], 
                   alpha=0.6, s=30, c='green')
        ax8.set_xlabel('Exact KL Value')
        ax8.set_ylabel('Relative Error (%)')
        ax8.set_title('KL Error vs Magnitude')
        ax8.grid(True, alpha=0.3)
    
    # === COMBINED PLOTS ===
    # 9. Error comparison (if both available)
    if len(valid_eucl) > 0 and len(valid_kl) > 0:
        # Find common pairs
        common_pairs = set(valid_eucl['pair_name']) & set(valid_kl['pair_name'])
        if len(common_pairs) > 0:
            ax9 = plt.subplot(3, 4, 9)
            
            eucl_subset = valid_eucl[valid_eucl['pair_name'].isin(common_pairs)]
            kl_subset = valid_kl[valid_kl['pair_name'].isin(common_pairs)]
            
            ax9.scatter(eucl_subset['euclidean_rel_error'], kl_subset['kl_rel_error'], 
                       alpha=0.6, s=30, c='purple')
            ax9.set_xlabel('Euclidean Relative Error (%)')
            ax9.set_ylabel('KL Relative Error (%)')
            ax9.set_title(f'Error Correlation\n({len(common_pairs)} common pairs)')
            ax9.grid(True, alpha=0.3)
    
    # 10. Overall timing summary
    valid_both = results_df.dropna(subset=['exact_total_time', 'sinkhorn_total_time'])
    if len(valid_both) > 0:
        ax10 = plt.subplot(3, 4, 10)
        times_total = [valid_both['exact_total_time'], valid_both['sinkhorn_total_time']]
        ax10.boxplot(times_total, labels=['Exact Total', 'Sinkhorn Total'])
        ax10.set_ylabel('Time (seconds)')
        ax10.set_title('Overall Timing Comparison')
        ax10.grid(True, alpha=0.3)
        ax10.set_yscale('log')
    
    # 11. Quality summary pie chart
    if len(valid_eucl) > 0:
        ax11 = plt.subplot(3, 4, 11)
        eucl_errors = valid_eucl['euclidean_rel_error']
        excellent = (eucl_errors < 5).sum()
        good = ((eucl_errors >= 5) & (eucl_errors < 10)).sum()
        acceptable = ((eucl_errors >= 10) & (eucl_errors < 20)).sum()
        poor = (eucl_errors >= 20).sum()
        
        sizes = [excellent, good, acceptable, poor]
        labels = ['Excellent\n(<5%)', 'Good\n(5-10%)', 'Acceptable\n(10-20%)', 'Poor\n(>20%)']
        colors = ['green', 'lightgreen', 'orange', 'red']
        
        ax11.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
        ax11.set_title('Euclidean Quality Distribution')
    
    # 12. Dataset size vs performance
    ax12 = plt.subplot(3, 4, 12)
    if len(valid_eucl) > 0:
        ax12.scatter(valid_eucl['n_data1'], valid_eucl['euclidean_rel_error'], 
                    alpha=0.6, s=30, c='blue', label='Euclidean')
    if len(valid_kl) > 0:
        ax12.scatter(valid_kl['n_data1'], valid_kl['kl_rel_error'], 
                    alpha=0.6, s=30, c='green', label='KL')
    ax12.set_xlabel('Dataset Size')
    ax12.set_ylabel('Relative Error (%)')
    ax12.set_title('Error vs Dataset Size')
    ax12.legend()
    ax12.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_path / 'comprehensive_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Comprehensive visualization saved to: {output_path / 'comprehensive_analysis.png'}")

def save_detailed_statistics(results_df, output_path):
    """Save detailed statistical analysis to file."""
    
    stats_file = output_path / 'detailed_statistics.txt'
    
    with open(stats_file, 'w') as f:
        f.write("COMPREHENSIVE SINKHORN vs EXACT WASSERSTEIN ANALYSIS\n")
        f.write("="*70 + "\n\n")
        f.write(f"Total dataset pairs tested: {len(results_df)}\n")
        f.write(f"Test completed: {pd.Timestamp.now()}\n\n")
        
        # Euclidean statistics
        valid_eucl = results_df.dropna(subset=['exact_euclidean', 'sinkhorn_euclidean'])
        if len(valid_eucl) > 0:
            f.write(f"EUCLIDEAN WASSERSTEIN RESULTS ({len(valid_eucl)} valid pairs):\n")
            f.write("-" * 50 + "\n")
            
            eucl_errors = valid_eucl['euclidean_rel_error']
            f.write(f"Relative Error Statistics:\n")
            f.write(f"  Count: {len(eucl_errors)}\n")
            f.write(f"  Mean: {eucl_errors.mean():.4f}%\n")
            f.write(f"  Median: {eucl_errors.median():.4f}%\n")
            f.write(f"  Std Dev: {eucl_errors.std():.4f}%\n")
            f.write(f"  Min: {eucl_errors.min():.4f}%\n")
            f.write(f"  Max: {eucl_errors.max():.4f}%\n")
            f.write(f"  25th percentile: {eucl_errors.quantile(0.25):.4f}%\n")
            f.write(f"  75th percentile: {eucl_errors.quantile(0.75):.4f}%\n")
            f.write(f"  95th percentile: {eucl_errors.quantile(0.95):.4f}%\n")
            f.write(f"  99th percentile: {eucl_errors.quantile(0.99):.4f}%\n\n")
            
            f.write(f"Correlation: {valid_eucl['exact_euclidean'].corr(valid_eucl['sinkhorn_euclidean']):.6f}\n")
            f.write(f"Exact timing: {valid_eucl['exact_euclidean_time'].mean():.4f} ± {valid_eucl['exact_euclidean_time'].std():.4f}s\n")
            f.write(f"Sinkhorn timing: {valid_eucl['sinkhorn_euclidean_time'].mean():.4f} ± {valid_eucl['sinkhorn_euclidean_time'].std():.4f}s\n")
            f.write(f"Speedup: {valid_eucl['exact_euclidean_time'].mean() / valid_eucl['sinkhorn_euclidean_time'].mean():.2f}x\n\n")
        
        # KL statistics
        valid_kl = results_df.dropna(subset=['exact_kl', 'sinkhorn_kl'])
        if len(valid_kl) > 0:
            f.write(f"KL WASSERSTEIN RESULTS ({len(valid_kl)} valid pairs):\n")
            f.write("-" * 50 + "\n")
            
            kl_errors = valid_kl['kl_rel_error']
            f.write(f"Relative Error Statistics:\n")
            f.write(f"  Count: {len(kl_errors)}\n")
            f.write(f"  Mean: {kl_errors.mean():.4f}%\n")
            f.write(f"  Median: {kl_errors.median():.4f}%\n")
            f.write(f"  Std Dev: {kl_errors.std():.4f}%\n")
            f.write(f"  Min: {kl_errors.min():.4f}%\n")
            f.write(f"  Max: {kl_errors.max():.4f}%\n")
            f.write(f"  25th percentile: {kl_errors.quantile(0.25):.4f}%\n")
            f.write(f"  75th percentile: {kl_errors.quantile(0.75):.4f}%\n")
            f.write(f"  95th percentile: {kl_errors.quantile(0.95):.4f}%\n")
            f.write(f"  99th percentile: {kl_errors.quantile(0.99):.4f}%\n\n")
            
            f.write(f"Correlation: {valid_kl['exact_kl'].corr(valid_kl['sinkhorn_kl']):.6f}\n")
            f.write(f"Exact timing: {valid_kl['exact_kl_time'].mean():.4f} ± {valid_kl['exact_kl_time'].std():.4f}s\n")
            f.write(f"Sinkhorn timing: {valid_kl['sinkhorn_kl_time'].mean():.4f} ± {valid_kl['sinkhorn_kl_time'].std():.4f}s\n")
            f.write(f"Speedup: {valid_kl['exact_kl_time'].mean() / valid_kl['sinkhorn_kl_time'].mean():.2f}x\n\n")
    
    print(f"Detailed statistics saved to: {stats_file}")

def main():
    """Main function to run the comprehensive test."""
    results_df = run_comprehensive_test()
    
    if results_df is not None:
        print(f"\n{'='*90}")
        print("COMPREHENSIVE TEST COMPLETED SUCCESSFULLY!")
        print(f"{'='*90}")
        print(f"Total results: {len(results_df)} dataset pairs")
        print("Check the 'comprehensive_sinkhorn_results' directory for detailed outputs.")
    else:
        print("Test failed - check error messages above.")

if __name__ == "__main__":
    main()
