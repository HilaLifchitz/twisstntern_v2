#!/usr/bin/env python3
"""
Comprehensive tuning for KL Wasserstein parameters to minimize bias.
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
    HAS_OT = True
except ImportError:
    print("Error: POT library not found. Please install: pip install POT")
    HAS_OT = False
    sys.exit(1)

from fast_topology_metrics import (
    load_topology_data_fast,
    compute_wasserstein_kl_cpu
)

def compute_sinkhorn_kl_custom(data1, data2, max_points=2000, reg=0.01, 
                              numItermax=1000, stopThr=1e-9, method='sinkhorn', 
                              verbose=False):
    """
    Custom KL Wasserstein with tunable Sinkhorn parameters.
    """
    if verbose:
        print(f"Computing KL Sinkhorn (reg={reg}, iter={numItermax}, method={method})...")
    
    start_time = time.time()
    
    try:
        # Sample data if needed
        if len(data1) > max_points:
            data1 = data1.sample(n=max_points, random_state=42)
        if len(data2) > max_points:
            data2 = data2.sample(n=max_points, random_state=42)
        
        # Convert to arrays
        data1_arr = data1[['T1', 'T2', 'T3']].values.astype(np.float64)
        data2_arr = data2[['T1', 'T2', 'T3']].values.astype(np.float64)
        
        # Add epsilon to avoid log(0)
        data1_eps = data1_arr + 1e-8
        data1_eps = data1_eps / data1_eps.sum(1, keepdims=True)
        
        data2_eps = data2_arr + 1e-8
        data2_eps = data2_eps / data2_eps.sum(1, keepdims=True)
        
        # Calculate symmetric KL divergence cost matrix
        kl_pq = (data1_eps[None, :, :] * (np.log(data1_eps[None, :, :]) - np.log(data2_eps[:, None, :]))).sum(-1)
        kl_qp = (data2_eps[None, :, :] * (np.log(data2_eps[None, :, :]) - np.log(data1_eps[:, None, :]))).sum(-1)
        M_kl = (kl_pq + kl_qp.T).astype(np.float64)
        
        # Uniform weights
        a = np.ones(len(data1_eps), dtype=np.float64) / len(data1_eps)
        b = np.ones(len(data2_eps), dtype=np.float64) / len(data2_eps)
        
        # Apply Sinkhorn with custom parameters
        distance = ot.sinkhorn2(a, b, M_kl, reg, 
                              numItermax=numItermax,
                              stopThr=stopThr, 
                              method=method,
                              warn=False)
        
        computation_time = time.time() - start_time
        
        if verbose:
            print(f"  KL Sinkhorn computed in {computation_time:.3f}s")
        
        return float(distance), computation_time
        
    except Exception as e:
        if verbose:
            print(f"  Error computing KL Sinkhorn: {e}")
        return np.nan, np.nan

def test_kl_parameter_combinations():
    """Test comprehensive combinations of KL Wasserstein parameters."""
    
    print("="*80)
    print("KL WASSERSTEIN PARAMETER TUNING")
    print("="*80)
    
    # Load test datasets
    print("Loading test datasets...")
    data1 = load_topology_data_fast("migration_topology_weights.csv", verbose=False)
    data2 = load_topology_data_fast("NoMigration_topology_weights.csv", verbose=False)
    
    # Sample to consistent size for testing
    data1 = data1.sample(n=1000, random_state=42)
    data2 = data2.sample(n=1000, random_state=42)
    
    # Compute exact reference
    print("Computing exact KL reference...")
    exact_kl = compute_wasserstein_kl_cpu(data1, data2, verbose=False)
    print(f"Exact KL: {exact_kl:.6f}")
    
    # Parameter combinations to test
    reg_values = [0.0001, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05]
    numItermax_values = [50, 100, 500, 1000, 2000]
    methods = ['sinkhorn', 'sinkhorn_stabilized', 'sinkhorn_epsilon_scaling']
    stopThr_values = [1e-6, 1e-7, 1e-8, 1e-9]
    sample_sizes = [1000, 2000]
    
    results = []
    total_tests = len(reg_values) * len(methods) * len(sample_sizes)
    test_count = 0
    
    print(f"\nTesting {total_tests} parameter combinations...")
    
    for sample_size in sample_sizes:
        print(f"\n--- Sample Size: {sample_size} ---")
        
        for method in methods:
            print(f"\n  Method: {method}")
            
            for reg in reg_values:
                test_count += 1
                
                try:
                    # Use default values for other parameters, focus on reg and method
                    distance, comp_time = compute_sinkhorn_kl_custom(
                        data1, data2, 
                        max_points=sample_size,
                        reg=reg,
                        numItermax=1000,  # Use higher iterations for accuracy
                        stopThr=1e-9,
                        method=method,
                        verbose=False
                    )
                    
                    if not np.isnan(distance):
                        error = abs(distance - exact_kl)
                        rel_error = (error / exact_kl) * 100
                        
                        results.append({
                            'sample_size': sample_size,
                            'method': method,
                            'reg': reg,
                            'numItermax': 1000,
                            'stopThr': 1e-9,
                            'sinkhorn_kl': distance,
                            'exact_kl': exact_kl,
                            'error': error,
                            'rel_error': rel_error,
                            'time': comp_time
                        })
                        
                        print(f"    reg={reg:6.4f}: KL={distance:.4f} ({rel_error:+5.1f}%), Time={comp_time:.3f}s")
                    else:
                        print(f"    reg={reg:6.4f}: FAILED")
                
                except Exception as e:
                    print(f"    reg={reg:6.4f}: ERROR - {e}")
                    continue
    
    # Now test numItermax and stopThr for best reg/method combinations
    if results:
        results_df = pd.DataFrame(results)
        
        # Find best reg/method combinations
        best_combos = results_df.nsmallest(3, 'rel_error')[['reg', 'method']].drop_duplicates()
        
        print(f"\n--- Testing iterations and thresholds for best combinations ---")
        
        for _, combo in best_combos.iterrows():
            reg = combo['reg']
            method = combo['method']
            
            print(f"\nTesting {method} with reg={reg}:")
            
            for numIter in numItermax_values:
                for stopThr in stopThr_values:
                    try:
                        distance, comp_time = compute_sinkhorn_kl_custom(
                            data1, data2,
                            max_points=1000,
                            reg=reg,
                            numItermax=numIter,
                            stopThr=stopThr,
                            method=method,
                            verbose=False
                        )
                        
                        if not np.isnan(distance):
                            error = abs(distance - exact_kl)
                            rel_error = (error / exact_kl) * 100
                            
                            results.append({
                                'sample_size': 1000,
                                'method': method,
                                'reg': reg,
                                'numItermax': numIter,
                                'stopThr': stopThr,
                                'sinkhorn_kl': distance,
                                'exact_kl': exact_kl,
                                'error': error,
                                'rel_error': rel_error,
                                'time': comp_time
                            })
                            
                            print(f"  iter={numIter:4d}, thr={stopThr:.0e}: "
                                  f"KL={distance:.4f} ({rel_error:+5.1f}%), Time={comp_time:.3f}s")
                    
                    except Exception as e:
                        continue
    
    # Convert to final DataFrame and analyze
    results_df = pd.DataFrame(results)
    
    if len(results_df) > 0:
        # Save results
        results_df.to_csv('kl_wasserstein_tuning_results.csv', index=False)
        
        # Analysis
        analyze_kl_results(results_df, exact_kl)
        create_kl_visualizations(results_df)
    
    return results_df

def analyze_kl_results(results_df, exact_kl):
    """Analyze KL Wasserstein tuning results."""
    
    print(f"\n{'='*80}")
    print("KL WASSERSTEIN TUNING ANALYSIS")
    print(f"{'='*80}")
    
    if len(results_df) == 0:
        print("No valid results to analyze!")
        return
    
    print(f"Analyzed {len(results_df)} successful parameter combinations")
    print(f"Exact KL reference: {exact_kl:.6f}")
    
    # Best overall accuracy
    best_accuracy = results_df.loc[results_df['rel_error'].abs().idxmin()]
    print(f"\nBEST ACCURACY:")
    print(f"  Method: {best_accuracy['method']}")
    print(f"  reg: {best_accuracy['reg']}")
    print(f"  numItermax: {best_accuracy['numItermax']}")
    print(f"  stopThr: {best_accuracy['stopThr']:.0e}")
    print(f"  Error: {best_accuracy['rel_error']:+.1f}%")
    print(f"  Time: {best_accuracy['time']:.3f}s")
    
    # Best speed
    fastest = results_df.loc[results_df['time'].idxmin()]
    print(f"\nFASTEST:")
    print(f"  Method: {fastest['method']}")
    print(f"  reg: {fastest['reg']}")
    print(f"  Error: {fastest['rel_error']:+.1f}%")
    print(f"  Time: {fastest['time']:.3f}s")
    
    # Good compromises (< 20% error, < 1s)
    good_compromises = results_df[
        (results_df['rel_error'].abs() < 20) & 
        (results_df['time'] < 1.0)
    ].nsmallest(5, 'rel_error')
    
    if len(good_compromises) > 0:
        print(f"\nBEST COMPROMISES (<20% error, <1s):")
        for _, row in good_compromises.iterrows():
            print(f"  {row['method']}, reg={row['reg']}, iter={row['numItermax']}: "
                  f"{row['rel_error']:+.1f}% error, {row['time']:.3f}s")
    
    # Analysis by method
    print(f"\nANALYSIS BY METHOD:")
    for method in results_df['method'].unique():
        method_data = results_df[results_df['method'] == method]
        best_method = method_data.loc[method_data['rel_error'].abs().idxmin()]
        print(f"  {method}: Best {best_method['rel_error']:+.1f}% error "
              f"(reg={best_method['reg']}, {best_method['time']:.3f}s)")
    
    # Analysis by regularization
    print(f"\nANALYSIS BY REGULARIZATION:")
    reg_analysis = results_df.groupby('reg').agg({
        'rel_error': ['mean', 'std', 'min'],
        'time': 'mean'
    }).round(2)
    print(reg_analysis)

def create_kl_visualizations(results_df):
    """Create visualizations for KL Wasserstein tuning."""
    
    if len(results_df) == 0:
        print("No data to visualize!")
        return
    
    plt.style.use('default')
    sns.set_palette("viridis")
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('KL Wasserstein Parameter Tuning Results', fontsize=16)
    
    # 1. Error vs Regularization by Method
    ax = axes[0, 0]
    for method in results_df['method'].unique():
        method_data = results_df[results_df['method'] == method]
        ax.semilogx(method_data['reg'], method_data['rel_error'], 'o-', label=method, alpha=0.7)
    ax.set_xlabel('Regularization Parameter')
    ax.set_ylabel('Relative Error (%)')
    ax.set_title('Error vs Regularization by Method')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
    
    # 2. Time vs Regularization by Method
    ax = axes[0, 1]
    for method in results_df['method'].unique():
        method_data = results_df[results_df['method'] == method]
        ax.semilogx(method_data['reg'], method_data['time'], 'o-', label=method, alpha=0.7)
    ax.set_xlabel('Regularization Parameter')
    ax.set_ylabel('Computation Time (s)')
    ax.set_title('Time vs Regularization by Method')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 3. Error vs Time Tradeoff
    ax = axes[0, 2]
    methods = results_df['method'].unique()
    colors = plt.cm.Set2(np.linspace(0, 1, len(methods)))
    
    for method, color in zip(methods, colors):
        method_data = results_df[results_df['method'] == method]
        ax.scatter(method_data['time'], method_data['rel_error'].abs(), 
                  label=method, alpha=0.7, s=60, c=[color])
    
    ax.set_xlabel('Computation Time (s)')
    ax.set_ylabel('|Relative Error| (%)')
    ax.set_title('Error vs Time Tradeoff')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_yscale('log')
    
    # 4. Error by Method (boxplot)
    ax = axes[1, 0]
    methods_for_box = [results_df[results_df['method'] == m]['rel_error'] for m in results_df['method'].unique()]
    ax.boxplot(methods_for_box, labels=results_df['method'].unique())
    ax.set_ylabel('Relative Error (%)')
    ax.set_title('Error Distribution by Method')
    ax.grid(True, alpha=0.3)
    ax.tick_params(axis='x', rotation=45)
    
    # 5. Heatmap: Error by reg and method
    ax = axes[1, 1]
    pivot_data = results_df.pivot_table(values='rel_error', 
                                       index='method', columns='reg', aggfunc='mean')
    
    sns.heatmap(pivot_data, annot=True, fmt='.1f', cmap='RdYlBu_r', ax=ax)
    ax.set_title('Mean Error Heatmap (%) by Method and Regularization')
    ax.set_xlabel('Regularization')
    ax.set_ylabel('Method')
    
    # 6. Best parameters scatter
    ax = axes[1, 2]
    # Show only best results (< 50% error)
    good_results = results_df[results_df['rel_error'].abs() < 50]
    
    if len(good_results) > 0:
        scatter = ax.scatter(good_results['reg'], good_results['rel_error'], 
                           c=good_results['time'], s=80, alpha=0.7, cmap='plasma')
        ax.set_xscale('log')
        ax.set_xlabel('Regularization Parameter')
        ax.set_ylabel('Relative Error (%)')
        ax.set_title('Best Results (colored by time)')
        plt.colorbar(scatter, ax=ax, label='Time (s)')
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig('kl_wasserstein_tuning_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"KL tuning visualization saved to: kl_wasserstein_tuning_analysis.png")

if __name__ == "__main__":
    if not HAS_OT:
        print("Error: POT library required")
        sys.exit(1)
    
    results = test_kl_parameter_combinations() 