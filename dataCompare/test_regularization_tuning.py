#!/usr/bin/env python3
"""
Test different regularization parameters for Sinkhorn to minimize bias.
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

# Import the fast topology metrics functions
from fast_topology_metrics import (
    load_topology_data_fast,
    compute_wasserstein_euclidean_cpu,
    compute_wasserstein_kl_cpu,
    compute_sinkhorn_wasserstein_fast
)

def test_regularization_parameters():
    """Test different regularization parameters to find optimal bias/speed tradeoff."""
    
    print("="*80)
    print("SINKHORN REGULARIZATION PARAMETER TUNING")
    print("="*80)
    
    # Load test datasets
    print("Loading test datasets...")
    data1 = load_topology_data_fast("migration_topology_weights.csv", verbose=False)
    data2 = load_topology_data_fast("NoMigration_topology_weights.csv", verbose=False)
    
    # Sample to 1000 points for consistent testing
    data1 = data1.sample(n=1000, random_state=42)
    data2 = data2.sample(n=1000, random_state=42)
    
    # Compute exact reference values
    print("Computing exact reference values...")
    exact_eucl = compute_wasserstein_euclidean_cpu(data1, data2, verbose=False)
    exact_kl = compute_wasserstein_kl_cpu(data1, data2, verbose=False)
    
    print(f"Exact Euclidean: {exact_eucl:.6f}")
    print(f"Exact KL: {exact_kl:.6f}")
    
    # Test different regularization parameters
    reg_values = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
    sample_sizes = [1000, 2000, 5000]
    
    results = []
    
    print(f"\nTesting {len(reg_values)} regularization values × {len(sample_sizes)} sample sizes...")
    
    for sample_size in sample_sizes:
        print(f"\n--- Sample Size: {sample_size} ---")
        
        for reg in reg_values:
            try:
                start_time = time.time()
                
                # Test Euclidean (current function only does Euclidean)
                sinkhorn_eucl = compute_sinkhorn_wasserstein_fast(
                    data1, data2, max_points=sample_size, reg=reg, verbose=False
                )
                
                # For KL, we'll need to implement a separate function or modify the existing one
                # For now, let's skip KL testing
                sinkhorn_kl = np.nan
                
                computation_time = time.time() - start_time
                
                # Calculate errors
                eucl_error = abs(sinkhorn_eucl - exact_eucl)
                eucl_rel_error = (eucl_error / exact_eucl) * 100
                
                kl_error = abs(sinkhorn_kl - exact_kl)
                kl_rel_error = (kl_error / exact_kl) * 100
                
                results.append({
                    'sample_size': sample_size,
                    'reg': reg,
                    'sinkhorn_eucl': sinkhorn_eucl,
                    'sinkhorn_kl': sinkhorn_kl,
                    'eucl_error': eucl_error,
                    'eucl_rel_error': eucl_rel_error,
                    'kl_error': kl_error,
                    'kl_rel_error': kl_rel_error,
                    'time': computation_time,
                    'exact_eucl': exact_eucl,
                    'exact_kl': exact_kl
                })
                
                print(f"  reg={reg:5.3f}: Eucl={sinkhorn_eucl:.4f} ({eucl_rel_error:+5.1f}%), "
                      f"KL={sinkhorn_kl:.4f} ({kl_rel_error:+5.1f}%), Time={computation_time:.3f}s")
                
            except Exception as e:
                print(f"  reg={reg:5.3f}: ERROR - {e}")
                continue
    
    results_df = pd.DataFrame(results)
    
    # Save results
    results_df.to_csv('regularization_tuning_results.csv', index=False)
    
    # Analysis
    print(f"\n{'='*80}")
    print("REGULARIZATION TUNING ANALYSIS")
    print(f"{'='*80}")
    
    if len(results_df) > 0:
        # Find best parameters (skip NaN values)
        valid_eucl = results_df.dropna(subset=['eucl_rel_error'])
        valid_kl = results_df.dropna(subset=['kl_rel_error'])
        
        if len(valid_eucl) > 0:
            best_eucl = valid_eucl.loc[valid_eucl['eucl_rel_error'].abs().idxmin()]
        else:
            best_eucl = None
            
        if len(valid_kl) > 0:
            best_kl = valid_kl.loc[valid_kl['kl_rel_error'].abs().idxmin()]
        else:
            best_kl = None
        
        print(f"\nBEST PARAMETERS:")
        if best_eucl is not None:
            print(f"Euclidean - reg={best_eucl['reg']}, size={best_eucl['sample_size']}, "
                  f"error={best_eucl['eucl_rel_error']:+.1f}%, time={best_eucl['time']:.3f}s")
        else:
            print("Euclidean - No valid results")
            
        if best_kl is not None:
            print(f"KL        - reg={best_kl['reg']}, size={best_kl['sample_size']}, "
                  f"error={best_kl['kl_rel_error']:+.1f}%, time={best_kl['time']:.3f}s")
        else:
            print("KL        - No valid results (skipped for this test)")
        
        # Show good compromises (< 10% error, < 1s time)
        good_eucl = results_df[(results_df['eucl_rel_error'].abs() < 10) & (results_df['time'] < 1.0)]
        good_kl = results_df[(results_df['kl_rel_error'].abs() < 20) & (results_df['time'] < 1.0)]
        
        if len(good_eucl) > 0:
            print(f"\nGOOD EUCLIDEAN COMPROMISES (<10% error, <1s):")
            for _, row in good_eucl.iterrows():
                print(f"  reg={row['reg']}, size={row['sample_size']}, "
                      f"error={row['eucl_rel_error']:+.1f}%, time={row['time']:.3f}s")
        
        if len(good_kl) > 0:
            print(f"\nGOOD KL COMPROMISES (<20% error, <1s):")
            for _, row in good_kl.iterrows():
                print(f"  reg={row['reg']}, size={row['sample_size']}, "
                      f"error={row['kl_rel_error']:+.1f}%, time={row['time']:.3f}s")
    
    # Create visualizations
    create_tuning_visualizations(results_df)
    
    return results_df

def create_tuning_visualizations(results_df):
    """Create visualizations of regularization tuning results."""
    
    if len(results_df) == 0:
        print("No data to visualize!")
        return
    
    # Set up plotting
    plt.style.use('default')
    sns.set_palette("husl")
    
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('Sinkhorn Regularization Parameter Tuning', fontsize=16)
    
    # 1. Error vs Regularization (Euclidean)
    ax = axes[0, 0]
    for size in results_df['sample_size'].unique():
        subset = results_df[results_df['sample_size'] == size]
        ax.semilogx(subset['reg'], subset['eucl_rel_error'], 'o-', label=f'{size} points')
    ax.set_xlabel('Regularization Parameter')
    ax.set_ylabel('Euclidean Relative Error (%)')
    ax.set_title('Euclidean Error vs Regularization')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
    
    # 2. Error vs Regularization (KL)
    ax = axes[0, 1]
    for size in results_df['sample_size'].unique():
        subset = results_df[results_df['sample_size'] == size]
        ax.semilogx(subset['reg'], subset['kl_rel_error'], 'o-', label=f'{size} points')
    ax.set_xlabel('Regularization Parameter')
    ax.set_ylabel('KL Relative Error (%)')
    ax.set_title('KL Error vs Regularization')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='red', linestyle='--', alpha=0.5)
    
    # 3. Time vs Regularization
    ax = axes[0, 2]
    for size in results_df['sample_size'].unique():
        subset = results_df[results_df['sample_size'] == size]
        ax.semilogx(subset['reg'], subset['time'], 'o-', label=f'{size} points')
    ax.set_xlabel('Regularization Parameter')
    ax.set_ylabel('Computation Time (s)')
    ax.set_title('Time vs Regularization')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 4. Error vs Time tradeoff (Euclidean)
    ax = axes[1, 0]
    scatter = ax.scatter(results_df['time'], results_df['eucl_rel_error'].abs(), 
                        c=results_df['reg'], s=60, alpha=0.7, cmap='viridis')
    ax.set_xlabel('Computation Time (s)')
    ax.set_ylabel('|Euclidean Relative Error| (%)')
    ax.set_title('Euclidean: Error vs Time Tradeoff')
    plt.colorbar(scatter, ax=ax, label='Regularization')
    ax.grid(True, alpha=0.3)
    
    # 5. Error vs Time tradeoff (KL)
    ax = axes[1, 1]
    scatter = ax.scatter(results_df['time'], results_df['kl_rel_error'].abs(), 
                        c=results_df['reg'], s=60, alpha=0.7, cmap='viridis')
    ax.set_xlabel('Computation Time (s)')
    ax.set_ylabel('|KL Relative Error| (%)')
    ax.set_title('KL: Error vs Time Tradeoff')
    plt.colorbar(scatter, ax=ax, label='Regularization')
    ax.grid(True, alpha=0.3)
    
    # 6. Heatmap of errors by reg and sample size
    ax = axes[1, 2]
    pivot_eucl = results_df.pivot_table(values='eucl_rel_error', 
                                       index='reg', columns='sample_size', aggfunc='mean')
    sns.heatmap(pivot_eucl, annot=True, fmt='.1f', cmap='RdBu_r', center=0, ax=ax)
    ax.set_title('Euclidean Error Heatmap (%)')
    ax.set_xlabel('Sample Size')
    ax.set_ylabel('Regularization')
    
    plt.tight_layout()
    plt.savefig('regularization_tuning_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"Visualization saved to: regularization_tuning_analysis.png")

if __name__ == "__main__":
    results = test_regularization_parameters() 