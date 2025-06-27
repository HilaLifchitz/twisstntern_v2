#!/usr/bin/env python3
"""
Compare two topology weights CSV files using twisstntern grid analysis and print L2, chi2, Wasserstein (Euclidean), and Wasserstein (KL) metrics.

Usage:
    python metrics.py file1.csv file2.csv [--granularity 0.1]
"""
import sys
import argparse
import numpy as np
import pandas as pd
from scipy.stats import chi2
try:
    import ot  # POT: Python Optimal Transport
except ImportError:
    print("Installing POT (Python Optimal Transport)...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "POT"])
    import ot
import time

from twisstntern.utils import dump_data
from twisstntern.analysis import triangles_analysis

def wasserstein_distance_euclidean(data1, data2):
    if isinstance(data1, pd.DataFrame):
        data1 = data1[['T1', 'T2', 'T3']].values
    if isinstance(data2, pd.DataFrame):
        data2 = data2[['T1', 'T2', 'T3']].values
    a = np.ones((data1.shape[0],)) / data1.shape[0]
    b = np.ones((data2.shape[0],)) / data2.shape[0]
    M = ot.dist(data1, data2)
    # Use Sinkhorn for speed (reg=0.01). Switch back to ot.emd2 for exact.
    dist = ot.sinkhorn2(a, b, M, reg=0.01)
    return dist

def simplex_wasserstein_distance_kl(data1, data2):
    if isinstance(data1, pd.DataFrame):
        data1 = data1[['T1', 'T2', 'T3']].values
    if isinstance(data2, pd.DataFrame):
        data2 = data2[['T1', 'T2', 'T3']].values
    data1_eps = data1 + 1e-8
    data1_eps = data1_eps / data1_eps.sum(1, keepdims=True)
    data2_eps = data2 + 1e-8
    data2_eps = data2_eps / data2_eps.sum(1, keepdims=True)
    kl_pq = (data1_eps[None, :, :] * (np.log(data1_eps[None, :, :]) - np.log(data2_eps[:, None, :]))).sum(-1)
    kl_qp = (data2_eps[None, :, :] * (np.log(data2_eps[None, :, :]) - np.log(data1_eps[:, None, :]))).sum(-1)
    kl_symm = (kl_pq + kl_qp.T)
    # Use Sinkhorn for speed (reg=0.01). Switch back to ot.emd2 for exact.
    dist = ot.sinkhorn2(np.ones(data2_eps.shape[0])/data2_eps.shape[0], 
                       np.ones(data1_eps.shape[0])/data1_eps.shape[0],
                       kl_symm, reg=0.01)
    return dist

def compute_metrics(data1, data2, granularity):
    """
    Compute L2, chi2, Wasserstein (Euclidean), and Wasserstein (KL) metrics between two datasets.
    Args:
        data1: DataFrame or array for first dataset
        data2: DataFrame or array for second dataset
        granularity: grid granularity (float or str)
    Returns:
        dict: metrics with keys 'l2', 'chi2_stat', 'p_value', 'wasserstein_euclidean', 'wasserstein_kl'
    """
    # Run grid analysis on full datasets
    grid1 = triangles_analysis(data1, granularity) # grid1 is the grid of the first dataset (data1) "data"
    grid2 = triangles_analysis(data2, granularity) # grid2 is the grid of the second dataset (data2) "model"

    # Match triangles by index
    merged = pd.merge(grid1, grid2, on="index", suffixes=("_1", "_2"))
    # Use n_right + n_left as count per triangle
    counts1 = merged["n_right_1"] + merged["n_left_1"]
    counts2 = merged["n_right_2"] + merged["n_left_2"]
    # Normalize to proportions
    prop1 = counts1 / counts1.sum()
    prop2 = counts2 / counts2.sum()

    # L2 distance
        # NOTICE: L2 distance under uniform measure: sqrt( (1 / len(I)) * sum over i of (f[i] - g[i])**2 ), I= index of partition of the simplex
    # This is the TRUE L2 distance between the two distributions
    # It treats your discrete space as a probability space with uniform measure.
    # It makes our norm stable under changes to partition granularity-- i.e. different alpha
   
    l2 = np.linalg.norm(prop1 - prop2)
    l2 = l2* 1/len(grid1) 


    # Chi-square
    # Standard Chi-square test statistic:
    #     chi2 = sum( (O_i - E_i)**2 / E_i )
    # where:
    #     O_i = observed count (e.g. from data)
    #     E_i = expected count (e.g. from model)
    with np.errstate(divide='ignore', invalid='ignore'):
        chi2_stat = np.nansum((counts1 - counts2) ** 2 / (counts1 + 1e-8))
    dof = (prop1 != 0).sum() - 1
    p_value = 1 - chi2.cdf(chi2_stat, dof)

    # Wasserstein distances - SAMPLE for speed!
    print("Sampling data for Wasserstein calculations (1000 points each)...")
    sample_size = 1000
    
    if isinstance(data1, pd.DataFrame):
        data1_sample = data1.sample(n=min(sample_size, len(data1)), random_state=42)
    else:
        indices = np.random.choice(len(data1), size=min(sample_size, len(data1)), replace=False)
        data1_sample = data1[indices]
    
    if isinstance(data2, pd.DataFrame):
        data2_sample = data2.sample(n=min(sample_size, len(data2)), random_state=42)
    else:
        indices = np.random.choice(len(data2), size=min(sample_size, len(data2)), replace=False)
        data2_sample = data2[indices]
    
    print(f"Computing Wasserstein on {len(data1_sample)} × {len(data2_sample)} points...")
    w_euclidean = wasserstein_distance_euclidean(data1_sample, data2_sample)
    w_kl = simplex_wasserstein_distance_kl(data1_sample, data2_sample)

    return {
        'l2': l2,
        'chi2_stat': chi2_stat,
        'p_value': p_value,
        'wasserstein_euclidean': w_euclidean,
        'wasserstein_kl': w_kl
    }

def main():
    parser = argparse.ArgumentParser(description='Compare two topology weights CSV files using twisstntern grid analysis and print L2, chi2, Wasserstein (Euclidean), and Wasserstein (KL) metrics.')
    parser.add_argument('file1', type=str, help='First CSV file')
    parser.add_argument('file2', type=str, help='Second CSV file')
    parser.add_argument('--granularity', type=float, default=0.1, help='Grid granularity (default: 0.1)')
    args = parser.parse_args()

    # Timing: Data loading
    t0 = time.time()
    data1 = dump_data(args.file1)
    data2 = dump_data(args.file2)
    t1 = time.time()
    print(f"Data loading took {t1-t0:.2f} seconds.")

    # Timing: Metric computation
    t2 = time.time()
    data1 = data1.astype(np.float32)
    data2 = data2.astype(np.float32)
    metrics = compute_metrics(data1, data2, args.granularity)
    t3 = time.time()
    print(f"Metric computation took {t3-t2:.2f} seconds.")

    # Print results
    print(f"L2_distance: {metrics['l2']:.8f}")
    print(f"chi2_statistic: {metrics['chi2_stat']:.8f}")
    print(f"p_value: {metrics['p_value']:.8e}")
    print(f"wasserstein_euclidean: {metrics['wasserstein_euclidean']:.8f}")
    print(f"wasserstein_kl: {metrics['wasserstein_kl']:.8f}")

if __name__ == "__main__":
    main() 