#!/usr/bin/env python3
"""
Metrics Computation

A streamlined script for efficiently computing comparison metrics between
two topology weight CSV files. Focuses on speed and accuracy.

Metrics computed:
- L² distance
- χ² statistic with correct degrees of freedom
- Wasserstein distance (Euclidean)
- Wasserstein distance (KL-divergence)



 🔧 How to Use This Script

 Option 1: Run from the command line (bash)
 -------------------------------------------
 Example:
   python metrics.py data1.csv data2.csv --alpha 0.1 --verbose

 Required arguments:
   - file1: path to first CSV file
   - file2: path to second CSV file

 Optional flags:
   --alpha                Grid granularity (default: 0.1)
   --max-wasserstein-points  Limit number of points for Wasserstein (default: 20000)
   --skip-wasserstein     Skip Wasserstein calculations for faster execution
   --verbose              Print detailed output during processing


 Option 2: Import and use in another Python script
 -------------------------------------------------
 You can also call the main comparison function directly:

 from metrics import compute_all_metrics, load_data

 data1 = load_data("data1.csv")
 data2 = load_data("data2.csv")

  FULL COMMAND:
 results = compute_all_metrics(data1, data2, alpha=0.05, verbose=True,max_wasserstein_points=20000,skip_wasserstein=False)
 # OR JUST:
results = compute_all_metrics(data1, data2, alpha=0.05, verbose=True)
L2, chi2, p_value, wasserstein, wasserstein_kl = results['L2_distance'], results['chi2_statistic'], results['p_value'], results['wasserstein_euclidean'], results['wasserstein_kl']
more compact:
L2, chi2, p_value, wasserstein, wasserstein_kl = (lambda r: (r['L2_distance'], r['chi2_statistic'], r['p_value'], r['wasserstein_euclidean'], r['wasserstein_kl']))(compute_all_metrics(data1, data2, alpha=alpha))






 Inputs:
   - file1, file2: CSV file paths with columns ['T1', 'T2', 'T3'] for ternary data
   - alpha: granularity of the triangle grid for L2/χ² comparison
   - max_wasserstein_points: cap on number of sampled points in Wasserstein computations
   - skip_wasserstein: if True, skip Wasserstein distances to save time
   - verbose: print detailed progress and results


 =============================================================================
"""

import sys
import argparse
import time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import chi2
import warnings

# Optimal transport - handle import gracefully with GPU support
try:
    import ot
    HAS_OT = True
    # Try to import cupy for GPU acceleration
    try:
        import cupy as cp
        HAS_GPU = True
        print("GPU acceleration available via CuPy")
    except ImportError:
        HAS_GPU = False
        print("GPU acceleration not available (CuPy not found)")
except ImportError:
    print("Warning: 'ot' library not found. Installing POT...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "POT"])
    import ot
    HAS_OT = True
    HAS_GPU = False

# Try importing torch for additional optimization
try:
    import torch
    HAS_TORCH = torch.cuda.is_available()
    if HAS_TORCH:
        print(f"PyTorch GPU acceleration available on {torch.cuda.get_device_name()}")
except ImportError:
    HAS_TORCH = False

# Try importing numba for JIT compilation speedup
try:
    from numba import jit, prange
    HAS_NUMBA = True
    print("Numba JIT compilation available for speed boost")
except ImportError:
    HAS_NUMBA = False
    # Define dummy jit decorator
    def jit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator
    prange = range

# Add parent directory for twisstntern import
sys.path.append(str(Path(__file__).parent.parent))
from twisstntern.utils import dump_data

# Suppress warnings for clean output
warnings.filterwarnings("ignore")


def load_data(filepath, verbose=False):
    """
    Fast data loading using twisstntern's dump_data.
    
    Args:
        filepath (str): Path to CSV file
        verbose (bool): Print loading info
        
    Returns:
        pd.DataFrame: Processed data with T1, T2, T3 columns
    """
    if verbose:
        print(f"Loading: {filepath}")
    
    start_time = time.time()
    
    try:
        # Use twisstntern's optimized dump_data
        data_array = dump_data(str(filepath))
        data = pd.DataFrame(data_array, columns=['T1', 'T2', 'T3'])
        
        if verbose:
            load_time = time.time() - start_time
            print(f"  Loaded {len(data)} points in {load_time:.3f}s")
        
        return data
        
    except Exception as e:
        if verbose:
            print(f"Error with dump_data: {e}")
            print("Falling back to manual loading...")
        
        # Fallback manual loading
        data = pd.read_csv(filepath)
        
        # Handle topology definitions row
        if len(data) > 1:
            first_row = data.iloc[0]
            if any(isinstance(val, str) and '(' in str(val) for val in first_row):
                data = data.iloc[1:]
        
        # Ensure T1, T2, T3 columns
        if 'T1' in data.columns and 'T2' in data.columns and 'T3' in data.columns:
            data = data[['T1', 'T2', 'T3']]
        else:
            data = data.iloc[:, :3]
            data.columns = ['T1', 'T2', 'T3']
        
        # Convert to numeric
        for col in ['T1', 'T2', 'T3']:
            data[col] = pd.to_numeric(data[col], errors='coerce')
        
        # Clean and normalize
        data = data.dropna()
        row_sums = data.sum(axis=1)
        data = data.div(row_sums, axis=0)
        data = data[data['T2'] != data['T3']]  # Remove y-axis points
        
        if verbose:
            load_time = time.time() - start_time
            print(f"  Loaded {len(data)} points in {load_time:.3f}s (manual)")
        
        return data


def count_points_in_triangle_fast(a1, b1, a2, b2, a3, b3, data):
    """
    Optimized point counting in triangle bounds.
    
    Args:
        a1, b1, a2, b2, a3, b3: Triangle bounds for T1, T2, T3
        data: DataFrame with T1, T2, T3 columns
        
    Returns:
        int: Number of points in triangle
    """
    # Vectorized conditions for speed
    if a1 == 0:
        cond_t1 = (data['T1'] >= a1) & (data['T1'] <= b1)
    else:
        cond_t1 = (data['T1'] > a1) & (data['T1'] <= b1)
    
    if a2 == 0:
        cond_t2 = (data['T2'] >= a2) & (data['T2'] <= b2)
    else:
        cond_t2 = (data['T2'] > a2) & (data['T2'] <= b2)
    
    if a3 == 0:
        cond_t3 = (data['T3'] >= a3) & (data['T3'] <= b3)
    else:
        cond_t3 = (data['T3'] > a3) & (data['T3'] <= b3)
    
    return (cond_t1 & cond_t2 & cond_t3).sum()


def generate_triangles_fast(alpha):
    """
    Fast triangle grid generation.
    
    Args:
        alpha (float): Grid granularity
        
    Returns:
        list: Triangle coordinate dictionaries
    """
    triangles = []
    steps = int(1 / alpha)
    
    for k in range(steps):
        a1 = round(k * alpha, 10)
        b1 = round((k + 1) * alpha, 10)
        
        T2_upper_limit = round(1 - k * alpha, 10)
        T2_steps = round(T2_upper_limit / alpha)
        
        a3_1 = round(1 - (k + 1) * alpha, 10)
        b3_1 = round(1 - k * alpha, 10)
        
        for T2_step in range(T2_steps):
            a2 = round(T2_step * alpha, 10)
            b2 = round((T2_step + 1) * alpha, 10)
            
            if a3_1 >= 0:
                triangles.append({
                    'T1': (a1, b1),
                    'T2': (a2, b2),
                    'T3': (a3_1, b3_1)
                })
            
            a3_2 = round(a3_1 - alpha, 10)
            b3_2 = round(b3_1 - alpha, 10)
            
            if a3_2 >= 0:
                triangles.append({
                    'T1': (a1, b1),
                    'T2': (a2, b2),
                    'T3': (a3_2, b3_2)
                })
            
            a3_1 = a3_2
            b3_1 = b3_2
    
    return triangles



def compute_l2_chi2(data1, data2, alpha, verbose=False):
    """
    Fast computation of L² distance and χ² statistic.
    
    Args:
        data1, data2: DataFrames with topology weights
        alpha: Grid granularity
        verbose: Print timing info
        
    Returns:
        tuple: (L2_distance, chi2_statistic, p_value, degrees_freedom, num_meaningful_triangles)
    """
    # if verbose:
    #     print(f"Computing L² and χ² with α={alpha}...")
    
    start_time = time.time()

    # Generate triangle grid
    triangles = generate_triangles_fast(alpha)

    # Count points in each triangle for both datasets
    counts1 = []
    counts2 = []

    for triangle in triangles:
        count1 = count_points_in_triangle_fast(
            triangle['T1'][0], triangle['T1'][1],
            triangle['T2'][0], triangle['T2'][1],
            triangle['T3'][0], triangle['T3'][1],
            data1
        )
        count2 = count_points_in_triangle_fast(
            triangle['T1'][0], triangle['T1'][1],
            triangle['T2'][0], triangle['T2'][1],
            triangle['T3'][0], triangle['T3'][1],
            data2
        )
        counts1.append(count1)
        counts2.append(count2)

    # Convert to numpy arrays for speed
    counts1 = np.array(counts1)
    counts2 = np.array(counts2)

    # Calculate proportions
    # when calculating  (L2(f,g))^2= \sum_{trinagles} |f-g|^2, f= distr of data1, g=dist of data2
    n1, n2 = len(data1), len(data2)
    props1 = counts1 / n1
    props2 = counts2 / n2

    # Filter meaningful triangles (at least one dataset has points)

    meaningful_mask = (counts1 > 0) | (counts2 > 0)
    # normalized counts- i.e. the values as distrbution per bin -- for the L2 distance
    meaningful_props1 = props1[meaningful_mask]
    meaningful_props2 = props2[meaningful_mask]

    # not normalized counts- i.e. the values as counts per bin -- for the chi2 statistic
    meaningful_counts1 = counts1[meaningful_mask]
    meaningful_counts2 = counts2[meaningful_mask]

    num_meaningful = meaningful_mask.sum()

    if num_meaningful == 0:
        return 0.0, np.nan, np.nan, 0, 0



    # L² distance (using meaningful triangles only)
    residuals_squared = (meaningful_props1 - meaningful_props2) ** 2
    L2_distance = np.sqrt(residuals_squared.sum() / num_meaningful)
    # L2 distance = sqrt[sum over i of (f[i] - g[i])**2 )* dmu(i)]  dmu(i)=1/#triangles - uniform distrbution
    # NOTICE: L2 distance under uniform measure: sqrt( (1 / len(I)) * sum over i of (f[i] - g[i])**2 ), I= index of partition of the simplex
    # This is the TRUE L2 distance between the two distributions
    # It treats your discrete space as a probability space with uniform measure.
    # It makes our norm stable under changes to partition granularity-- i.e. different alpha
    L2_distance = L2_distance * np.sqrt(1/len(triangles))


    # χ² statistic (only for triangles with non-zero expected counts)

    # Standard χ² test statistic:
    #     chi2 = sum( (O_i - E_i)**2 / E_i )
    # where:
    #     O_i = observed count (e.g. from data <-> data1) -- not normalized! not proportions!
    #     E_i = expected count (e.g. from model <-> data2) -- not normalized! not proportions!

    # Check for problematic bins (E_i = 0, O_i ≠ 0)
    problematic_bins = (meaningful_counts1 == 0) & (meaningful_counts2 != 0)
    n_problematic = np.sum(problematic_bins)
    if n_problematic > 0:
        print("⚠️ Warning: Chi-square test assumptions are violated.")
        print(f"   → {n_problematic} bin(s) have expected count E_i = 0 but observed count O_i ≠ 0.")
        print("   → Treat the result as diagnostic rather than a formal chi-square test.")


    if meaningful_mask.any():
        with np.errstate(divide='ignore', invalid='ignore'):
            chi2_stat = np.nansum((meaningful_counts1 - meaningful_counts2) ** 2 / (meaningful_counts1 + 1e-8))
            degrees_freedom = num_meaningful - 1
            from scipy.stats import chi2
            p_value = 1 - chi2.cdf(chi2_stat, degrees_freedom)
    else:
        chi2_stat = np.nan
        p_value = np.nan
        degrees_freedom = 0

    if verbose:
        elapsed = time.time() - start_time

        print(f"  L² and χ² computed in {elapsed:.3f}s")
        print(f"  Meaningful triangles: {num_meaningful}/{len(triangles)}")


    return (
        float(L2_distance),
        float(chi2_stat),
        float(p_value),
        int(degrees_freedom)
    )




def compute_wasserstein_euclidean_gpu(data1, data2, verbose=False):
    """
    GPU-accelerated Wasserstein distance computation with Euclidean cost.
    
    Args:
        data1, data2: DataFrames with topology weights
        verbose: Print timing info
        
    Returns:
        float: Wasserstein distance
    """
    if verbose:
        print("Computing Wasserstein (Euclidean) with GPU acceleration...")
    
    start_time = time.time()
    
    try:
        if HAS_TORCH:
            # Use PyTorch GPU acceleration
            device = torch.device("cuda")
            
            # Convert to torch tensors on GPU
            arr1 = torch.tensor(data1[['T1', 'T2', 'T3']].values, dtype=torch.float32, device=device)
            arr2 = torch.tensor(data2[['T1', 'T2', 'T3']].values, dtype=torch.float32, device=device)
            
            # Compute pairwise distances on GPU
            # Use batched computation to handle memory efficiently
            batch_size = 2000
            n1, n2 = len(arr1), len(arr2)
            
            # Initialize distance matrix
            M = torch.zeros((n1, n2), dtype=torch.float32, device=device)
            
            for i in range(0, n1, batch_size):
                end_i = min(i + batch_size, n1)
                batch1 = arr1[i:end_i]
                
                for j in range(0, n2, batch_size):
                    end_j = min(j + batch_size, n2)
                    batch2 = arr2[j:end_j]
                    
                    # Compute Euclidean distances for this batch
                    diff = batch1.unsqueeze(1) - batch2.unsqueeze(0)  # Broadcasting
                    dist_batch = torch.sqrt(torch.sum(diff ** 2, dim=2))
                    M[i:end_i, j:end_j] = dist_batch
            
            # Move back to CPU for OT computation
            M_cpu = M.cpu().numpy()
            
        elif HAS_GPU:
            # Use CuPy acceleration
            arr1_gpu = cp.array(data1[['T1', 'T2', 'T3']].values, dtype=cp.float32)
            arr2_gpu = cp.array(data2[['T1', 'T2', 'T3']].values, dtype=cp.float32)
            
            # Compute distance matrix on GPU
            M_gpu = cp.sqrt(cp.sum((arr1_gpu[:, cp.newaxis, :] - arr2_gpu[cp.newaxis, :, :]) ** 2, axis=2))
            M_cpu = cp.asnumpy(M_gpu)
            
        else:
            # Fallback to CPU with optimized numpy
            arr1 = data1[['T1', 'T2', 'T3']].values.astype(np.float32)
            arr2 = data2[['T1', 'T2', 'T3']].values.astype(np.float32)
            M_cpu = ot.dist(arr1, arr2, metric='euclidean')
        
        # Uniform weights
        a = np.ones(len(data1), dtype=np.float32) / len(data1)
        b = np.ones(len(data2), dtype=np.float32) / len(data2)
        
        # Solve optimal transport on CPU
        distance = ot.emd2(a, b, M_cpu, numItermax=1000000)
        
        if verbose:
            elapsed = time.time() - start_time
            print(f"  Wasserstein (Euclidean) computed in {elapsed:.3f}s with GPU")
        
        return float(distance)
        
    except Exception as e:
        if verbose:
            print(f"  GPU computation failed: {e}")
            print("  Falling back to CPU computation...")
        # Fallback to CPU version
        return compute_wasserstein_euclidean_cpu(data1, data2, verbose)


def compute_wasserstein_euclidean_cpu(data1, data2, max_points=None, verbose=False):
    """
    Fast Wasserstein distance computation with Euclidean cost.
    
    Args:
        data1, data2: DataFrames with topology weights
        max_points: Maximum points to use (for speed), None for all
        verbose: Print timing info
        
    Returns:
        float: Wasserstein distance
    """
    if verbose:
        print("Computing Wasserstein (Euclidean)...")
    
    start_time = time.time()
    
    try:
        # Sample data if requested for speed
        if max_points is not None:
            if len(data1) > max_points:
                data1 = data1.sample(n=max_points, random_state=42)
                if verbose:
                    print(f"  Sampled data1 to {max_points} points")
            if len(data2) > max_points:
                data2 = data2.sample(n=max_points, random_state=42)
                if verbose:
                    print(f"  Sampled data2 to {max_points} points")
        
        # Convert to arrays
        arr1 = data1[['T1', 'T2', 'T3']].values.astype(np.float32)
        arr2 = data2[['T1', 'T2', 'T3']].values.astype(np.float32)
        
        # Uniform weights
        a = np.ones(len(arr1), dtype=np.float32) / len(arr1)
        b = np.ones(len(arr2), dtype=np.float32) / len(arr2)
        
        # Euclidean cost matrix
        M = ot.dist(arr1, arr2, metric='euclidean')
        
        # Solve optimal transport
        distance = ot.emd2(a, b, M, numItermax=1000000)
        
        if verbose:
            elapsed = time.time() - start_time
            print(f"  Wasserstein (Euclidean) computed in {elapsed:.3f}s")
        
        return float(distance)
        
    except Exception as e:
        if verbose:
            print(f"  Error computing Wasserstein (Euclidean): {e}")
        return np.nan


"""
wasserstein_kl:

    Compute the Earth Mover's Distance (EMD) between two sets of points on the ternary simplex,
    using symmetric Kullback-Leibler (KL) divergence as the ground cost.

    Parameters
    ----------
    data1 : np.ndarray or pd.DataFrame of shape (n1, 3)
        First distribution on the simplex, where each row is a point (T1, T2, T3)
        summing approximately to 1. Represents discrete probability weights.
    
    data2 : np.ndarray or pd.DataFrame of shape (n2, 3)
        Second distribution on the simplex, same format as `data1`.

    Returns
    -------
    dist : float
        The optimal transport cost (EMD) between the two distributions,
        where the cost of transporting mass from a point in `data1` to a point in `data2`
        is given by the symmetric KL divergence between them.

    Notes
    -----
    - This is not a classical Wasserstein distance using Euclidean geometry.
    - Instead, it uses the symmetric KL divergence:
          KL(p || q) + KL(q || p)
      as a "distance" (cost) between points p and q on the simplex.
    - This function applies the EMD algorithm (`ot.emd2`) from the POT (Python Optimal Transport) package.
    - Each input distribution is assumed to have uniform weights over its support (1/n).

    Use case: Comparing two empirical distributions of topology weights (e.g. from TWISST)
    over the simplex in population genetics or phylogenetics.
    """

def compute_wasserstein_kl_gpu(data1, data2, verbose=False):
    """
    GPU-accelerated Wasserstein distance computation with KL-divergence cost.
    
    Args:
        data1, data2: DataFrames with topology weights
        verbose: Print timing info
        
    Returns:
        float: Wasserstein distance
    """
    if verbose:
        print("Computing Wasserstein (KL-divergence) with GPU acceleration...")
    
    start_time = time.time()
    
    try:
        if HAS_TORCH:
            # Use PyTorch GPU acceleration
            device = torch.device("cuda")
            
            # Convert to torch tensors and normalize
            arr1 = torch.tensor(data1[['T1', 'T2', 'T3']].values, dtype=torch.float32, device=device) + 1e-8
            arr2 = torch.tensor(data2[['T1', 'T2', 'T3']].values, dtype=torch.float32, device=device) + 1e-8
            
            # Renormalize
            arr1 = arr1 / torch.sum(arr1, dim=1, keepdim=True)
            arr2 = arr2 / torch.sum(arr2, dim=1, keepdim=True)
            
            # Compute symmetric KL divergence in batches to manage memory
            batch_size = 1000  # Smaller batches for KL computation
            n1, n2 = len(arr1), len(arr2)
            
            kl_symm = torch.zeros((n1, n2), dtype=torch.float32, device=device)
            
            for i in range(0, n1, batch_size):
                end_i = min(i + batch_size, n1)
                batch1 = arr1[i:end_i]
                
                for j in range(0, n2, batch_size):
                    end_j = min(j + batch_size, n2)
                    batch2 = arr2[j:end_j]
                    
                    # Compute KL divergences for this batch
                    # KL(P||Q) = sum(P * log(P/Q))
                    log_ratio_pq = torch.log(batch1.unsqueeze(1) / batch2.unsqueeze(0))
                    log_ratio_qp = torch.log(batch2.unsqueeze(0) / batch1.unsqueeze(1))
                    
                    kl_pq = torch.sum(batch1.unsqueeze(1) * log_ratio_pq, dim=2)
                    kl_qp = torch.sum(batch2.unsqueeze(0) * log_ratio_qp, dim=2).transpose(0, 1)
                    
                    kl_symm[i:end_i, j:end_j] = kl_pq + kl_qp
            
            # Move to CPU for OT computation
            kl_symm_cpu = kl_symm.cpu().numpy()
            
        elif HAS_GPU:
            # Use CuPy acceleration
            arr1_gpu = cp.array(data1[['T1', 'T2', 'T3']].values, dtype=cp.float32) + 1e-8
            arr2_gpu = cp.array(data2[['T1', 'T2', 'T3']].values, dtype=cp.float32) + 1e-8
            
            # Renormalize
            arr1_gpu = arr1_gpu / cp.sum(arr1_gpu, axis=1, keepdims=True)
            arr2_gpu = arr2_gpu / cp.sum(arr2_gpu, axis=1, keepdims=True)
            
            # Compute symmetric KL divergence
            kl_pq = cp.sum(arr1_gpu[:, cp.newaxis, :] * 
                          (cp.log(arr1_gpu[:, cp.newaxis, :]) - cp.log(arr2_gpu[cp.newaxis, :, :])), axis=2)
            kl_qp = cp.sum(arr2_gpu[cp.newaxis, :, :] * 
                          (cp.log(arr2_gpu[cp.newaxis, :, :]) - cp.log(arr1_gpu[:, cp.newaxis, :])), axis=2)
            kl_symm_gpu = kl_pq + kl_qp.T
            kl_symm_cpu = cp.asnumpy(kl_symm_gpu)
            
        else:
            # Fallback to CPU computation
            return compute_wasserstein_kl_cpu(data1, data2, None, verbose)
        
        # Uniform weights - fix dimension ordering
        a = np.ones(len(data2), dtype=np.float32) / len(data2)  # rows of kl_symm
        b = np.ones(len(data1), dtype=np.float32) / len(data1)  # cols of kl_symm
        
        # Solve optimal transport on CPU
        distance = ot.emd2(a, b, kl_symm_cpu, numItermax=int(3e6))
        
        if verbose:
            elapsed = time.time() - start_time
            print(f"  Wasserstein (KL) computed in {elapsed:.3f}s with GPU")
        
        return float(distance)
        
    except Exception as e:
        if verbose:
            print(f"  GPU computation failed: {e}")
            print("  Falling back to CPU computation...")
        # Fallback to CPU version
        return compute_wasserstein_kl_cpu(data1, data2, None, verbose)


def compute_wasserstein_kl_cpu(data1, data2, max_points=None, verbose=False):
    """
    Fast Wasserstein distance computation with KL-divergence cost.
    
    Args:
        data1, data2: DataFrames with topology weights
        max_points: Maximum points to use (for speed), None for all
        verbose: Print timing info
        
    Returns:
        float: Wasserstein distance
    """
    if verbose:
        print("Computing Wasserstein (KL-divergence) on CPU...")
    
    start_time = time.time()
    
    try:
        # Sample data if requested for speed
        if max_points is not None:
            if len(data1) > max_points:
                data1 = data1.sample(n=max_points, random_state=42)
                if verbose:
                    print(f"  Sampled data1 to {max_points} points")
            if len(data2) > max_points:
                data2 = data2.sample(n=max_points, random_state=42)
                if verbose:
                    print(f"  Sampled data2 to {max_points} points")
        
        # Convert to arrays and add epsilon
        arr1 = data1[['T1', 'T2', 'T3']].values.astype(np.float32) + 1e-8
        arr2 = data2[['T1', 'T2', 'T3']].values.astype(np.float32) + 1e-8
        
        # Renormalize
        arr1 = arr1 / arr1.sum(axis=1, keepdims=True)
        arr2 = arr2 / arr2.sum(axis=1, keepdims=True)
        
        # Symmetric KL divergence cost matrix
        # This is the computational bottleneck - optimize it
        kl_pq = (arr1[None, :, :] * (np.log(arr1[None, :, :]) - np.log(arr2[:, None, :]))).sum(-1)
        kl_qp = (arr2[None, :, :] * (np.log(arr2[None, :, :]) - np.log(arr1[:, None, :]))).sum(-1)
        kl_symm = kl_pq + kl_qp.T
        
        # Uniform weights - fix dimension ordering to match kl_symm
        a = np.ones(len(arr2), dtype=np.float32) / len(arr2)  # rows of kl_symm
        b = np.ones(len(arr1), dtype=np.float32) / len(arr1)  # cols of kl_symm
        
        # Solve optimal transport with corrected dimensions
        distance = ot.emd2(a, b, kl_symm, numItermax=int(3e6))
        
        if verbose:
            elapsed = time.time() - start_time
            print(f"  Wasserstein (KL) computed in {elapsed:.3f}s")
        
        return float(distance)
        
    except Exception as e:
        if verbose:
            print(f"  Error computing Wasserstein (KL): {e}")
        return np.nan


def compute_sinkhorn_wasserstein_fast(data1, data2, max_points=2000, reg=0.01, verbose=False):
    """
    Ultra-fast Sinkhorn Wasserstein approximation (Euclidean distance).
    
    Args:
        data1, data2: DataFrames with topology weights
        max_points: Maximum points to use (default: 2000 for speed)
        reg: Regularization parameter (default: 0.01 - optimal for Euclidean)
        verbose: Print timing info
        
    Returns:
        float: Approximate Wasserstein distance
    """
    if verbose:
        print(f"Computing Sinkhorn Wasserstein Euclidean (reg={reg}, max_points={max_points})...")
    
    start_time = time.time()
    
    try:
        # Aggressive sampling for speed
        if len(data1) > max_points:
            data1 = data1.sample(n=max_points, random_state=42)
        if len(data2) > max_points:
            data2 = data2.sample(n=max_points, random_state=42)
        
        # Convert to arrays
        arr1 = data1[['T1', 'T2', 'T3']].values.astype(np.float32)
        arr2 = data2[['T1', 'T2', 'T3']].values.astype(np.float32)
        
        # Uniform weights
        a = np.ones(len(arr1), dtype=np.float32) / len(arr1)
        b = np.ones(len(arr2), dtype=np.float32) / len(arr2)
        
        # Euclidean cost matrix
        M = ot.dist(arr1, arr2, metric='euclidean')
        
        # Fast Sinkhorn approximation with optimal parameters
        distance = ot.sinkhorn2(a, b, M, reg, numItermax=1000)
        
        if verbose:
            elapsed = time.time() - start_time
            print(f"  Sinkhorn Wasserstein Euclidean computed in {elapsed:.3f}s")
        
        return float(distance)
        
    except Exception as e:
        if verbose:
            print(f"  Error computing Sinkhorn Euclidean: {e}")
        return np.nan


def compute_sinkhorn_kl_fast(data1, data2, max_points=2000, reg=0.005, verbose=False):
    """
    Ultra-fast Sinkhorn KL Wasserstein approximation with optimal parameters.
    
    Args:
        data1, data2: DataFrames with topology weights
        max_points: Maximum points to use (default: 2000 for speed)
        reg: Regularization parameter (default: 0.005 - optimal for KL)
        verbose: Print timing info
        
    Returns:
        float: Approximate KL Wasserstein distance
    """
    if verbose:
        print(f"Computing Sinkhorn KL Wasserstein (reg={reg}, max_points={max_points})...")
    
    start_time = time.time()
    
    try:
        # Aggressive sampling for speed
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
        
        # Fast Sinkhorn approximation with optimal KL parameters
        distance = ot.sinkhorn2(a, b, M_kl, reg, numItermax=50, method='sinkhorn')
        
        if verbose:
            elapsed = time.time() - start_time
            print(f"  Sinkhorn KL Wasserstein computed in {elapsed:.3f}s")
        
        return float(distance)
        
    except Exception as e:
        if verbose:
            print(f"  Error computing Sinkhorn KL: {e}")
        return np.nan


@jit(nopython=True, parallel=True, cache=True)
def count_points_vectorized(T1_vals, T2_vals, T3_vals, bounds_array):
    """
    Ultra-fast vectorized point counting using Numba JIT.
    
    Args:
        T1_vals, T2_vals, T3_vals: Point coordinates as numpy arrays
        bounds_array: Array of triangle bounds [N, 6] where columns are [a1,b1,a2,b2,a3,b3]
        
    Returns:
        counts: Array of point counts for each triangle
    """
    n_triangles = bounds_array.shape[0]
    n_points = len(T1_vals)
    counts = np.zeros(n_triangles, dtype=np.int32)
    
    for i in prange(n_triangles):
        a1, b1, a2, b2, a3, b3 = bounds_array[i]
        count = 0
        
        for j in range(n_points):
            t1, t2, t3 = T1_vals[j], T2_vals[j], T3_vals[j]
            
            # Check bounds with same logic as original
            if a1 == 0:
                cond_t1 = (t1 >= a1) and (t1 <= b1)
            else:
                cond_t1 = (t1 > a1) and (t1 <= b1)
            
            if a2 == 0:
                cond_t2 = (t2 >= a2) and (t2 <= b2)
            else:
                cond_t2 = (t2 > a2) and (t2 <= b2)
            
            if a3 == 0:
                cond_t3 = (t3 >= a3) and (t3 <= b3)
            else:
                cond_t3 = (t3 > a3) and (t3 <= b3)
            
            if cond_t1 and cond_t2 and cond_t3:
                count += 1
        
        counts[i] = count
    
    return counts


def compute_l2_chi2_ultrafast(data1, data2, alpha, verbose=False):
    """
    Ultra-fast L² and χ² computation using vectorized operations.
    """
    if verbose:
        print(f"Computing L² and χ² with vectorized operations (α={alpha})...")
    
    start_time = time.time()
    
    # Generate triangle grid
    triangles = generate_triangles_fast(alpha)
    
    # Convert triangles to bounds array for vectorized counting
    bounds_array = np.array([[t['T1'][0], t['T1'][1], t['T2'][0], t['T2'][1], 
                             t['T3'][0], t['T3'][1]] for t in triangles], dtype=np.float64)
    
    # Extract coordinate arrays
    T1_vals1 = data1['T1'].values.astype(np.float64)
    T2_vals1 = data1['T2'].values.astype(np.float64)
    T3_vals1 = data1['T3'].values.astype(np.float64)
    
    T1_vals2 = data2['T1'].values.astype(np.float64)
    T2_vals2 = data2['T2'].values.astype(np.float64)
    T3_vals2 = data2['T3'].values.astype(np.float64)
    
    # Count points using vectorized function
    if HAS_NUMBA:
        counts1 = count_points_vectorized(T1_vals1, T2_vals1, T3_vals1, bounds_array)
        counts2 = count_points_vectorized(T1_vals2, T2_vals2, T3_vals2, bounds_array)
    else:
        # Fallback to original method if numba not available
        counts1 = []
        counts2 = []
        for triangle in triangles:
            count1 = count_points_in_triangle_fast(
                triangle['T1'][0], triangle['T1'][1],
                triangle['T2'][0], triangle['T2'][1],
                triangle['T3'][0], triangle['T3'][1],
                data1
            )
            count2 = count_points_in_triangle_fast(
                triangle['T1'][0], triangle['T1'][1],
                triangle['T2'][0], triangle['T2'][1],
                triangle['T3'][0], triangle['T3'][1],
                data2
            )
            counts1.append(count1)
            counts2.append(count2)
        counts1 = np.array(counts1)
        counts2 = np.array(counts2)
    


   
    # Calculate proportions
    # when calculating  (L2(f,g))^2= \sum_{trinagles} |f-g|^2, f= distr of data1, g=dist of data2
    n1, n2 = len(data1), len(data2)
    props1 = counts1 / n1
    props2 = counts2 / n2

    # Filter meaningful triangles (at least one dataset has points)

    meaningful_mask = (counts1 > 0) | (counts2 > 0)
    # normalized counts- i.e. the values as distrbution per bin -- for the L2 distance
    meaningful_props1 = props1[meaningful_mask]
    meaningful_props2 = props2[meaningful_mask]

    # not normalized counts- i.e. the values as counts per bin -- for the chi2 statistic
    meaningful_counts1 = counts1[meaningful_mask]
    meaningful_counts2 = counts2[meaningful_mask]

    num_meaningful = meaningful_mask.sum()

    if num_meaningful == 0:
        return 0.0, np.nan, np.nan, 0, 0


    # L² distance (using meaningful triangles only)
    residuals_squared = (meaningful_props1 - meaningful_props2) ** 2
    L2_distance = np.sqrt(residuals_squared.sum() / num_meaningful)
    # L2 distance = sqrt[sum over i of (f[i] - g[i])**2 )* dmu(i)]  dmu(i)=1/#triangles - uniform distrbution
    # NOTICE: L2 distance under uniform measure: sqrt( (1 / len(I)) * sum over i of (f[i] - g[i])**2 ), I= index of partition of the simplex
    # This is the TRUE L2 distance between the two distributions
    # It treats your discrete space as a probability space with uniform measure.
    # It makes our norm stable under changes to partition granularity-- i.e. different alpha
    L2_distance = L2_distance * np.sqrt(1/len(triangles))


    # χ² statistic (only for triangles with non-zero expected counts)

    # Standard χ² test statistic:
    #     chi2 = sum( (O_i - E_i)**2 / E_i )
    # where:
    #     O_i = observed count (e.g. from data <-> data1) -- not normalized! not proportions!
    #     E_i = expected count (e.g. from model <-> data2) -- not normalized! not proportions!

    # Check for problematic bins (E_i = 0, O_i ≠ 0)
    problematic_bins = (meaningful_counts1 == 0) & (meaningful_counts2 != 0)
    n_problematic = np.sum(problematic_bins)
    if n_problematic > 0:
        print("⚠️ Warning: Chi-square test assumptions are violated.")
        print(f"   → {n_problematic} bin(s) have expected count E_i = 0 but observed count O_i ≠ 0.")
        print("   → Treat the result as diagnostic rather than a formal chi-square test.")


    if meaningful_mask.any():
        with np.errstate(divide='ignore', invalid='ignore'):
            chi2_stat = np.nansum((meaningful_counts1 - meaningful_counts2) ** 2 / (meaningful_counts1 + 1e-8))
            degrees_freedom = num_meaningful - 1
            from scipy.stats import chi2
            p_value = 1 - chi2.cdf(chi2_stat, degrees_freedom)
    else:
        chi2_stat = np.nan
        p_value = np.nan
        degrees_freedom = 0

    if verbose:
        elapsed = time.time() - start_time

        print(f"  L² and χ² computed in {elapsed:.3f}s")
        print(f"  Meaningful triangles: {num_meaningful}/{len(triangles)}")


    return (
        float(L2_distance),
        float(chi2_stat),
        float(p_value),
        int(degrees_freedom)
    )


def compute_all_metrics(data1, data2, alpha=0.1, max_wasserstein_points=30000, skip_wasserstein=False, verbose=False):
    """
    Compute all comparison metrics between two topology weight files.
    
    Args:
        data1, data2: Dataframes with columns ['T1', 'T2', 'T3'] for ternary data
        alpha: Grid granularity
        verbose: Print progress and timing
        
    Returns:
        dict: All computed metrics
    """
    total_start = time.time()
    
    if verbose:
        print("="*60)
        print("METRICS COMPUTATION")
        print("="*60)
        print(f"Data 1: {data1}")
        print(f"Data 2: {data2}")
        print(f"Granularity: {alpha}")
        print()
    
    
    if verbose:
        print()
    
    # Compute L² and χ² metrics with ultra-fast vectorized operations
    L2_dist, chi2_stat, p_val, dof = compute_l2_chi2(
        data1, data2, alpha, verbose)
    
    # Compute Wasserstein distances if not skipped
    if skip_wasserstein:
        wasserstein_eucl = np.nan
        wasserstein_kl = np.nan
        if verbose:
            print("Skipping Wasserstein distances for faster computation")
    else:
        # Smart selection of computation method
        n1, n2 = len(data1), len(data2)
        total_points = n1 + n2
        
        # Use GPU if available and datasets are large enough to benefit
        # GPU overhead is worth it for datasets > 5000 points each
        use_gpu = (HAS_GPU or HAS_TORCH) and n1 > 1000 and n2 > 1000
        
        # Choose computation method based on dataset size and speed requirements
        if total_points > 15000:  # Use ultra-fast Sinkhorn for large datasets
            if verbose:
                print(f"Large dataset ({total_points} points) - using ultra-fast Sinkhorn approximation")
            # Use aggressive sampling for maximum speed
            sinkhorn_points = min(2000, max_wasserstein_points // 5)
            wasserstein_eucl = compute_sinkhorn_wasserstein_fast(data1, data2, sinkhorn_points, 0.1, verbose)
            wasserstein_kl = compute_sinkhorn_kl_fast(data1, data2, sinkhorn_points, 0.05, verbose)
        elif use_gpu and total_points <= 20000:  # GPU for medium datasets
            if verbose:
                print(f"Using GPU acceleration for {n1} x {n2} point comparison")
            wasserstein_eucl = compute_wasserstein_euclidean_gpu(data1, data2, verbose)
            wasserstein_kl = compute_wasserstein_kl_gpu(data1, data2, verbose)
        else:
            # Use CPU with moderate sampling for smaller datasets
            if total_points > 8000:
                max_points = min(max_wasserstein_points // 2, 5000)  # More aggressive sampling
                if verbose:
                    print(f"Medium dataset ({total_points} points) - using CPU with {max_points} point sampling")
            else:
                max_points = None
                if verbose:
                    print(f"Small dataset - using exact CPU computation for {n1} x {n2} points")
            
            wasserstein_eucl = compute_wasserstein_euclidean_cpu(data1, data2, max_points, verbose)
            wasserstein_kl = compute_wasserstein_kl_cpu(data1, data2, max_points, verbose)
    
    # Compile results
    metrics = {
        'L2_distance': L2_dist,
        'chi2_statistic': chi2_stat,
        'p_value': p_val,
        'degrees_freedom': dof,
        'wasserstein_euclidean': wasserstein_eucl,
        'wasserstein_kl': wasserstein_kl,
        'n_data1': len(data1),
        'n_data2': len(data2),
        'granularity': alpha
    }
    
    if verbose:
        total_time = time.time() - total_start
        print()
        print("="*60)
        print("RESULTS")
        print("="*60)
        print(f"L² distance:              {metrics['L2_distance']:.6f}")
        print(f"χ² statistic:             {metrics['chi2_statistic']:.6f}")
        print(f"χ² p-value:               {metrics['p_value']:.2e}")
        print(f"Degrees of freedom:       {metrics['degrees_freedom']}")
        print(f"Wasserstein (Euclidean):  {metrics['wasserstein_euclidean']:.6f}")
        print(f"Wasserstein (KL):         {metrics['wasserstein_kl']:.6f}")
        print()
        print("DATA SUMMARY:")
        print(f"Data 1 points:            {metrics['n_data1']:,}")
        print(f"Data 2 points:            {metrics['n_data2']:,}")
        print(f"Meaningful triangles:     {metrics['degrees_freedom']+1}")
        print(f"Total computation time:   {total_time:.3f}s")
        print("="*60)
    
    return metrics





if __name__ == "__main__":
  
    parser = argparse.ArgumentParser(description="Compute metrics between two datasets.")
    parser.add_argument("file1", help="Path to first CSV file-- assumed as the 'truth'/data/expected counts")
    parser.add_argument("file2", help="Path to second CSV file-- assumed as the 'test'/model/observed counts")
    parser.add_argument("--alpha", type=float, default=0.1, help="Grid granularity")
    parser.add_argument("--max-wasserstein-points", type=int, default=20000, help="Maximum points for Wasserstein computation")
    parser.add_argument("--skip-wasserstein", action="store_true", help="Skip Wasserstein distance calculations for faster computation")
    parser.add_argument("--verbose", action="store_true", help="Print verbose output")
    args = parser.parse_args()

    data1 = load_data(args.file1)
    data2 = load_data(args.file2)

    results = compute_all_metrics(data1, data2,alpha=args.alpha,max_wasserstein_points=args.max_wasserstein_points,skip_wasserstein=args.skip_wasserstein,verbose=args.verbose)

    print("📊 Metric Results:")
    for key, value in results.items():
        print(f"{key}: {value}")


# def main():
#     """Main function with command-line interface."""
#     parser = argparse.ArgumentParser(description='Fast topology metrics computation')
#     parser.add_argument('file1', nargs='?', help='First topology weights CSV file')
#     parser.add_argument('file2', nargs='?', help='Second topology weights CSV file')
#     parser.add_argument('--granularity', '-g', type=float, default=0.1,
#                        help='Grid granularity (default: 0.1)')
#     parser.add_argument('--quiet', '-q', action='store_true',
#                        help='Suppress verbose output')
#     parser.add_argument('--output', '-o', type=str,
#                        help='Save results to JSON file')
#     parser.add_argument('--max-wasserstein-points', type=int, default=10000,
#                        help='Maximum points for Wasserstein computation (default: 10000)')
#     parser.add_argument('--skip-wasserstein', action='store_true',
#                        help='Skip Wasserstein distance calculations for faster computation')
#     parser.add_argument('--force-gpu', action='store_true',
#                        help='Force GPU usage even for smaller datasets')
#     parser.add_argument('--gpu-info', action='store_true',
#                        help='Show GPU availability information and exit')
#     parser.add_argument('--ultra-fast', action='store_true',
#                        help='Use ultra-fast Sinkhorn approximation for all Wasserstein computations')
#     parser.add_argument('--sinkhorn-points', type=int, default=2000,
#                        help='Number of points for Sinkhorn approximation (default: 2000)')
    
#     args = parser.parse_args()
    
#     # Show GPU info if requested
#     if args.gpu_info:
#         print("GPU Acceleration Information:")
#         print("=" * 40)
#         print(f"CuPy available: {HAS_GPU}")
#         print(f"PyTorch GPU available: {HAS_TORCH}")
#         if HAS_TORCH:
#             print(f"GPU device: {torch.cuda.get_device_name()}")
#             print(f"GPU memory: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
#         if HAS_GPU:
#             try:
#                 print(f"CuPy device: {cp.cuda.runtime.getDeviceProperties(0)['name'].decode()}")
#             except:
#                 print("CuPy device info not available")
#         if not (HAS_GPU or HAS_TORCH):
#             print("No GPU acceleration available")
#             print("Install PyTorch with CUDA or CuPy for GPU support")
#         sys.exit(0)
    
#     # Check that files are provided when not using --gpu-info
#     if not args.file1 or not args.file2:
#         print("Error: Both file1 and file2 are required")
#         parser.print_help()
#         sys.exit(1)
    
#     # Check files exist
#     file1_path = Path(args.file1)
#     file2_path = Path(args.file2)
    
#     if not file1_path.exists():
#         print(f"Error: File not found: {file1_path}")
#         sys.exit(1)
#     if not file2_path.exists():
#         print(f"Error: File not found: {file2_path}")
#         sys.exit(1)
    
#     # Compute metrics
#     verbose = not args.quiet
#     metrics = compute_all_metrics_fast(file1_path, file2_path, 
#                                      args.granularity, args.max_wasserstein_points,
#                                      args.skip_wasserstein, verbose)
    
#     # Save results if requested
#     if args.output:
#         import json
#         output_path = Path(args.output)
        
#         # Convert numpy types to native Python types for JSON serialization
#         json_metrics = {}
#         for key, value in metrics.items():
#             if isinstance(value, np.ndarray):
#                 json_metrics[key] = value.tolist()
#             elif isinstance(value, (np.integer, np.floating)):
#                 json_metrics[key] = value.item()
#             else:
#                 json_metrics[key] = value
        
#         with open(output_path, 'w') as f:
#             json.dump(json_metrics, f, indent=2)
        
#         if verbose:
#             print(f"\nResults saved to: {output_path}")
    
#     return metrics





# if __name__ == "__main__":
#     main() 