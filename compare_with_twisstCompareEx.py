#!/usr/bin/env python3
"""
Compare results from metrics.py with twisstCompareEx.py functions
"""
import sys
import numpy as np
import pandas as pd
from pathlib import Path

# Add the current directory to path to import twisstCompareEx functions
sys.path.append(str(Path(__file__).parent))

# Import functions from twisstCompareEx
from twisstCompareEx import (
    wasserstein_distance,
    simplex_wasserstein_distance,
    compare_triangle_counts,
    SubtrianglesDataCount
)

from twisstntern.utils import dump_data

def main():
    # Same files as used in metrics.py
    file1 = "/home/hlifchit/projects/twissting_baby/dataCompare/NoMigration_topology_weights_trimmed.csv"
    file2 = "/home/hlifchit/projects/twissting_baby/dataCompare/NoMigration_topology_weights.csv"
    alpha = 0.1
    
    print("Loading data using twisstntern.utils.dump_data...")
    
    # Load data using same method as metrics.py
    data1 = dump_data(file1)
    data2 = dump_data(file2)
    
    # Convert to DataFrames for twisstCompareEx functions
    df1 = pd.DataFrame(data1, columns=['T1', 'T2', 'T3'])
    df2 = pd.DataFrame(data2, columns=['T1', 'T2', 'T3'])
    
    print(f"Data1 shape: {df1.shape}")
    print(f"Data2 shape: {df2.shape}")
    
    print("\n" + "="*60)
    print("TWISSTCOMPAREEX.PY RESULTS")
    print("="*60)
    
    # Calculate Wasserstein distances using twisstCompareEx functions
    print("\nCalculating Wasserstein distances...")
    
    # Euclidean Wasserstein
    w_euclidean_tc = wasserstein_distance(df1.values, df2.values)
    print(f"Wasserstein (Euclidean): {w_euclidean_tc:.8f}")
    
    # KL Wasserstein
    w_kl_tc = simplex_wasserstein_distance(df1.values, df2.values)
    print(f"Wasserstein (KL): {w_kl_tc:.8f}")
    
    # Calculate triangle-based metrics
    print("\nCalculating triangle-based metrics...")
    
    # Get triangle counts for both datasets
    triangles1 = SubtrianglesDataCount(df1, alpha)
    triangles2 = SubtrianglesDataCount(df2, alpha)
    
    # Compare triangle counts
    comparison = compare_triangle_counts(df1, df2, alpha)
    
    # Calculate L2 distance from triangle proportions
    total_points1 = triangles1["nTruth"].sum()
    total_points2 = triangles2["nTruth"].sum()
    
    # Normalize counts to proportions
    prop1 = triangles1["nTruth"] / total_points1
    prop2 = triangles2["nTruth"] / total_points2
    
    # L2 distance
    l2_tc = np.linalg.norm(prop1 - prop2)
    print(f"L2 distance: {l2_tc:.8f}")
    
    # Chi-square test
    # Use the comparison results which has both counts
    counts1 = comparison["nTruth"]
    counts2 = comparison["nTruth.1"]  # This is the second dataset count
    
    # Chi-square statistic
    with np.errstate(divide='ignore', invalid='ignore'):
        chi2_stat_tc = np.nansum((counts1 - counts2) ** 2 / (counts1 + counts2 + 1e-8))
    
    # Degrees of freedom
    dof = (prop1 != 0).sum() - 1
    from scipy.stats import chi2
    p_value_tc = 1 - chi2.cdf(chi2_stat_tc, dof)
    
    print(f"Chi-square statistic: {chi2_stat_tc:.8f}")
    print(f"P-value: {p_value_tc:.8e}")
    
    print("\n" + "="*60)
    print("COMPARISON SUMMARY")
    print("="*60)
    print("METRICS.PY RESULTS:")
    print("  L2_distance: 0.02407842")
    print("  chi2_statistic: 6299.47021840")
    print("  p_value: 0.00000000e+00")
    print("  wasserstein_euclidean: 0.00021630")
    print("  wasserstein_kl: 0.00095829")
    
    print("\nTWISSTCOMPAREEX.PY RESULTS:")
    print(f"  L2_distance: {l2_tc:.8f}")
    print(f"  chi2_statistic: {chi2_stat_tc:.8f}")
    print(f"  p_value: {p_value_tc:.8e}")
    print(f"  wasserstein_euclidean: {w_euclidean_tc:.8f}")
    print(f"  wasserstein_kl: {w_kl_tc:.8f}")
    
    print("\nDIFFERENCES:")
    print(f"  L2 difference: {abs(0.02407842 - l2_tc):.8f}")
    print(f"  Chi2 difference: {abs(6299.47021840 - chi2_stat_tc):.8f}")
    print(f"  Wasserstein Euclidean difference: {abs(0.00021630 - w_euclidean_tc):.8f}")
    print(f"  Wasserstein KL difference: {abs(0.00095829 - w_kl_tc):.8f}")

if __name__ == "__main__":
    main() 