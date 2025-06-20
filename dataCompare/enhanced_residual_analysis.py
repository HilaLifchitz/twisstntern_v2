#!/usr/bin/env python3
"""
Enhanced Residual Analysis for Topology Weights Comparison

This enhanced version uses more functions from the twisstntern library and 
includes Wasserstein distance calculations from twisstCompareEx.py while 
maintaining the same visualization outputs as the original residual_analysis.py.

Key enhancements:
- Uses twisstntern.utils functions for data loading and processing
- Uses twisstntern plotting functions for grid creation  
- Integrates Wasserstein distance calculations (Euclidean and KL-divergence)
- Uses SubtrianglesDataCount and compare_triangle_counts from twisstCompareEx.py
- Prints comprehensive comparison metrics including Wasserstein distances

Usage:
    python enhanced_residual_analysis.py [--granularity 0.1] [--output results/]
"""

import sys
import os
import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import ListedColormap, Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
from scipy import stats
from scipy.stats import chi2
try:
    import ot  # optimal transport library for Wasserstein
    HAS_OT = True
except ImportError:
    print("Warning: 'ot' (optimal transport) library not found. Installing...")
    import subprocess
    import sys
    subprocess.check_call([sys.executable, "-m", "pip", "install", "POT"])
    import ot
    HAS_OT = True

import math

# Add the parent directory to the path to import twisstntern
sys.path.append(str(Path(__file__).parent.parent))
import twisstntern
import twisstntern.utils
from twisstntern.utils import dump_data

# Suppress warnings for cleaner output
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# Set modern vibrant plotting parameters
sns.set_style("darkgrid")
sns.set_palette("bright")
# For now, disable LaTeX to avoid issues and use Unicode symbols instead
plt.rcParams.update({
    "text.usetex": False,  # Disable LaTeX for now
    "font.family": "DejaVu Sans", 
    "font.size": 11,
    "axes.spines.left": True,
    "axes.spines.bottom": True,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.alpha": 0.4,
    "figure.facecolor": "white"
})
HAS_LATEX = False


def latex_format(text, bold=False):
    """Format text for LaTeX or fallback to regular formatting."""
    if HAS_LATEX:
        if bold:
            return fr'$\mathbf{{{text}}}$'
        else:
            return fr'${text}$'
    else:
        if bold:
            return text
        else:
            return text


def math_format(text):
    """Format mathematical expressions."""
    if HAS_LATEX:
        return fr'${text}$'
    else:
        # Replace common LaTeX symbols with Unicode
        text = text.replace(r'\mu', 'μ')
        text = text.replace(r'\chi^2', 'χ²')
        text = text.replace(r'L^2', 'L²')
        text = text.replace(r'\alpha', 'α')
        return text


def load_topology_data_twiss(filepath, sample_size=None):
    """
    Load and process topology weights CSV file using twisstntern.utils.dump_data.
    
    Args:
        filepath (str): Path to the CSV file
        sample_size (int, optional): Sample this many points if specified
        
    Returns:
        pd.DataFrame: Processed data with T1, T2, T3 columns
    """
    print(f"Loading data using twisstntern.utils.dump_data from: {filepath}")
    
    try:
        # Use twisstntern's dump_data function
        data_array = dump_data(str(filepath))
        
        # Convert to DataFrame
        data = pd.DataFrame(data_array, columns=['T1', 'T2', 'T3'])
        
        # Sample if requested
        if sample_size and len(data) > sample_size:
            data = data.sample(n=sample_size, random_state=42)
            print(f"Sampled {sample_size} points from {len(data_array)} total points")
        
        print(f"Loaded {len(data)} data points after processing")
        return data
        
    except Exception as e:
        print(f"Error loading with twisstntern.utils.dump_data: {e}")
        print("Falling back to manual CSV loading...")
        
        # Fallback to manual loading
        data = pd.read_csv(filepath)
        
        # Skip the topology definitions row (row 1) and use only numeric data
        if len(data) > 1:
            first_row = data.iloc[0]
            if any(isinstance(val, str) and '(' in str(val) for val in first_row):
                data = data.iloc[1:]
        
        # Ensure we have T1, T2, T3 columns
        if 'T1' in data.columns and 'T2' in data.columns and 'T3' in data.columns:
            data = data[['T1', 'T2', 'T3']]
        else:
            data = data.iloc[:, :3]
            data.columns = ['T1', 'T2', 'T3']
        
        # Convert to numeric and handle any issues
        for col in ['T1', 'T2', 'T3']:
            data[col] = pd.to_numeric(data[col], errors='coerce')
        
        # Remove any rows with NaN values
        data = data.dropna()
        
        # Normalize each row to sum to 1
        row_sums = data.sum(axis=1)
        data = data.div(row_sums, axis=0)
        
        # Remove points where T2 == T3 (on y-axis)
        data = data[data['T2'] != data['T3']]
        
        # Sample if requested
        if sample_size and len(data) > sample_size:
            data = data.sample(n=sample_size, random_state=42)
            print(f"Sampled {sample_size} points from total points")
        
        print(f"Loaded {len(data)} data points after processing")
        return data


# ============================================================================
# FUNCTIONS FROM TWISSTCOMPAREEX.PY
# ============================================================================

def wasserstein_distance_euclidean(data1, data2):
    """Wasserstein distance with Euclidean cost matrix from twisstCompareEx.py"""
    # Convert to numpy arrays if needed
    if isinstance(data1, pd.DataFrame):
        data1 = data1[['T1', 'T2', 'T3']].values
    if isinstance(data2, pd.DataFrame):
        data2 = data2[['T1', 'T2', 'T3']].values
    
    # Uniform weights for each point
    a = np.ones((data1.shape[0],)) / data1.shape[0]
    b = np.ones((data2.shape[0],)) / data2.shape[0]
    
    # Cost matrix: Euclidean distance between all points
    M = ot.dist(data1, data2)
    
    # Compute the Wasserstein distance
    dist = ot.emd2(a, b, M, numItermax=1000000)
    
    return dist


def simplex_wasserstein_distance_kl(data1, data2):
    """Simplex Wasserstein distance with KL cost matrix from twisstCompareEx.py"""
    # Convert to numpy arrays if needed
    if isinstance(data1, pd.DataFrame):
        data1 = data1[['T1', 'T2', 'T3']].values
    if isinstance(data2, pd.DataFrame):
        data2 = data2[['T1', 'T2', 'T3']].values
    
    # Adding an epsilon value to avoid division by zero in the KL calculation
    data1_plus_eps = data1 + 1e-8
    data1_plus_eps = data1_plus_eps / data1_plus_eps.sum(1, keepdims=True)

    data2_plus_eps = data2 + 1e-8
    data2_plus_eps = data2_plus_eps / data2_plus_eps.sum(1, keepdims=True)

    # Calculating the symmetric KL distance between each data point
    kl_pq = (data1_plus_eps[None, :, :] * (np.log(data1_plus_eps[None, :, :]) - np.log(data2_plus_eps[:, None, :]))).sum(-1)
    kl_qp = (data2_plus_eps[None, :, :] * (np.log(data2_plus_eps[None, :, :]) - np.log(data1_plus_eps[:, None, :]))).sum(-1)
    kl_symm = (kl_pq + kl_qp.T)

    # Solve optimal transport problem
    dist = ot.emd2(np.ones(data2_plus_eps.shape[0])/data2_plus_eps.shape[0], 
            np.ones(data1_plus_eps.shape[0])/data1_plus_eps.shape[0],
            M=kl_symm, numItermax=int(3e6))
    
    return dist


def n_twisstcompare(a1, b1, a2, b2, a3, b3, data):
    """Count function from twisstCompareEx.py - auxiliary function to count points in subtriangles"""
    # This additional logic is implemented to avoid double-counting of data
    if a1 == 0:
        condition_a1 = a1 <= data.T1
    else:
        condition_a1 = a1 < data.T1
    
    if a2 == 0:
        condition_a2 = a2 <= data.T2
    else:
        condition_a2 = a2 < data.T2
        
    if a3 == 0:
        condition_a3 = a3 <= data.T3
    else:
        condition_a3 = a3 < data.T3
        
    n = len(data[(condition_a1 & (data.T1 <= b1)) &
                 (condition_a2 & (data.T2 <= b2)) &
                 (condition_a3 & (data.T3 <= b3))])
    return n


# ============================================================================
# TRIANGLE GRID FUNCTIONS USING TWISSTNTERN
# ============================================================================

def create_triangular_grid_twiss(alpha):
    """
    Create triangular grid using twisstntern coordinate system.
    
    Args:
        alpha (float): Grid granularity
        
    Returns:
        list: List of triangle coordinate dictionaries
    """
    triangles = []
    steps = int(1 / alpha)
    
    for k in range(steps):
        a1 = round(k * alpha, 10)
        b1 = round((k + 1) * alpha, 10)
        
        # T2 goes from 0 to (1 - (k+1)*alpha) in steps of alpha
        T2_upper_limit = round(1 - k * alpha, 10)
        T2_steps = round(T2_upper_limit / alpha)
        
        # First triangle (T3 from [1 - (k+1)*α, 1 - k*α])
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
            
            # Second triangle - T3 is -alpha from the last coords
            a3_2 = round(a3_1 - alpha, 10)
            b3_2 = round(b3_1 - alpha, 10)
            
            if a3_2 >= 0:
                triangles.append({
                    'T1': (a1, b1),
                    'T2': (a2, b2),
                    'T3': (a3_2, b3_2)
                })
            
            # Update the T3 coordinates for the next T2 step
            a3_1 = a3_2
            b3_1 = b3_2
    
    return triangles


def perform_enhanced_grid_analysis(data1, data2, alpha):
    """
    Perform enhanced grid-based analysis using twisstntern and twisstCompareEx functions.
    
    Args:
        data1 (pd.DataFrame): First dataset (reference/truth)
        data2 (pd.DataFrame): Second dataset (model)
        alpha (float): Grid granularity
        
    Returns:
        tuple: (results_df, statistics_dict)
    """
    print(f"Performing enhanced grid analysis with α = {alpha}")
    
    # Create triangular grid using twisstntern coordinate system
    triangles = create_triangular_grid_twiss(alpha)
    
    # Calculate Wasserstein distances first
    print("Calculating Wasserstein distances...")
    try:
        wasserstein_euclidean = wasserstein_distance_euclidean(data1, data2)
        print(f"✓ Wasserstein (Euclidean): {wasserstein_euclidean:.6f}")
    except Exception as e:
        print(f"Error calculating Wasserstein (Euclidean): {e}")
        wasserstein_euclidean = np.nan
    
    try:
        wasserstein_kl = simplex_wasserstein_distance_kl(data1, data2)
        print(f"✓ Wasserstein (KL-divergence): {wasserstein_kl:.6f}")
    except Exception as e:
        print(f"Error calculating Wasserstein (KL): {e}")
        wasserstein_kl = np.nan
    
    # Count points in each triangle using enhanced counting
    results_df = []
    
    for triangle in triangles:
        # Count points in this triangle for both datasets
        count_data = n_twisstcompare(
            triangle['T1'][0], triangle['T1'][1],
            triangle['T2'][0], triangle['T2'][1],
            triangle['T3'][0], triangle['T3'][1],
            data1
        )
        
        count_model = n_twisstcompare(
            triangle['T1'][0], triangle['T1'][1],
            triangle['T2'][0], triangle['T2'][1],
            triangle['T3'][0], triangle['T3'][1],
            data2
        )
        
        # Calculate proportions
        prop_data = count_data / len(data1)
        prop_model = count_model / len(data2)
        prop_residual = prop_data - prop_model
        residual_squared = prop_residual ** 2
        
        results_df.append({
            'T1_bounds': triangle['T1'],
            'T2_bounds': triangle['T2'],
            'T3_bounds': triangle['T3'],
            'count_data': count_data,
            'count_model': count_model,
            'prop_data': prop_data,
            'prop_model': prop_model,
            'prop_residual': prop_residual,
            'residual_squared': residual_squared
        })
    
    results_df = pd.DataFrame(results_df)
    
    # Filter out triangles where both datasets have 0 points for statistics
    # These are "empty vs empty" comparisons that don't provide meaningful information
    meaningful_triangles = (results_df['count_data'] > 0) | (results_df['count_model'] > 0)
    filtered_results = results_df[meaningful_triangles]
    
    print(f"Total triangles: {len(results_df)}")
    print(f"Meaningful triangles (at least one dataset has points): {len(filtered_results)}")
    print(f"Empty triangles (both datasets have 0 points): {len(results_df) - len(filtered_results)}")
    
    # Calculate L2 distance and chi-square statistics using only meaningful triangles
    if len(filtered_results) > 0:
        L2_distance = np.sqrt(filtered_results['residual_squared'].sum() / len(filtered_results))
        
        # Chi-square test (only for triangles with non-zero expected counts)
        # Corrected: degrees of freedom should be number of meaningful triangles - 1
        non_zero_expected = filtered_results['count_data'] > 0
        if non_zero_expected.sum() > 0:
            chi_squared_components = filtered_results.loc[non_zero_expected, 'residual_squared'] * (len(data1) ** 2) / filtered_results.loc[non_zero_expected, 'count_data']
            chi_statistic = chi_squared_components.sum()
            degrees_freedom = len(filtered_results) - 1  # Corrected: meaningful triangles - 1
            p_value = chi2.sf(chi_statistic, degrees_freedom) if degrees_freedom > 0 else 1.0
        else:
            chi_statistic = np.nan
            p_value = np.nan
        
        # Calculate summary statistics using filtered data
        mean_residual = filtered_results['prop_residual'].mean()
        std_residual = filtered_results['prop_residual'].std()
        max_residual = filtered_results['prop_residual'].abs().max()
        mean_squared_residual = filtered_results['residual_squared'].mean()
    else:
        # Fallback if no meaningful triangles (shouldn't happen in practice)
        L2_distance = 0.0
        chi_statistic = np.nan
        p_value = np.nan
        mean_residual = 0.0
        std_residual = 0.0
        max_residual = 0.0
        mean_squared_residual = 0.0
    
    # Calculate summary statistics
    statistics = {
        'L2_distance': L2_distance,
        'chi2_statistic': chi_statistic,
        'p_value': p_value,
        'wasserstein_euclidean': wasserstein_euclidean,
        'wasserstein_kl': wasserstein_kl,
        'total_data_points': len(data1),
        'total_model_points': len(data2),
        'num_triangles': len(results_df),
        'num_meaningful_triangles': len(filtered_results),
        'num_empty_triangles': len(results_df) - len(filtered_results),
        'degrees_freedom': len(filtered_results) - 1,
        'mean_residual': mean_residual,
        'std_residual': std_residual,
        'max_residual': max_residual,
        'mean_squared_residual': mean_squared_residual
    }
    
    # Print comprehensive comparison metrics
    print("\n" + "="*70)
    print("COMPREHENSIVE COMPARISON METRICS")
    print("="*70)
    print(f"L² Distance:              {statistics['L2_distance']:.6f}")
    print(f"χ² Statistic:             {statistics['chi2_statistic']:.6f}")
    print(f"χ² p-value:               {statistics['p_value']:.6e}")
    print(f"Wasserstein (Euclidean):  {statistics['wasserstein_euclidean']:.6f}")
    print(f"Wasserstein (KL):         {statistics['wasserstein_kl']:.6f}")
    print(f"Mean Residual:            {statistics['mean_residual']:.6f}")
    print(f"Std Residual:             {statistics['std_residual']:.6f}")
    print(f"Max |Residual|:           {statistics['max_residual']:.6f}")
    print(f"Mean Squared Residual:    {statistics['mean_squared_residual']:.6f}")
    print(f"Data Points:              {statistics['total_data_points']}")
    print(f"Model Points:             {statistics['total_model_points']}")
    print(f"Grid Triangles (Total):   {statistics['num_triangles']}")
    print(f"Meaningful Triangles:     {statistics['num_meaningful_triangles']}")
    print(f"Empty Triangles:          {statistics['num_empty_triangles']}")
    print(f"Degrees of Freedom:       {statistics['degrees_freedom']}")
    print("="*70)
    
    return results_df, statistics


# ============================================================================
# ENHANCED PLOTTING FUNCTIONS USING TWISSTNTERN
# ============================================================================

def plot_ternary_base_twiss(ax, alpha):
    """
    Create ternary plot base using twisstntern functions.
    
    Args:
        ax: Matplotlib axis
        alpha: Grid granularity
    """
    h = math.sqrt(3) / 2
    
    # Plot triangle outline
    x_side_T2 = np.linspace(0, 0.5, 100)
    x_side_T3 = np.linspace(-0.5, 0, 100)
    
    ax.plot(x_side_T2, twisstntern.utils.T2(0, x_side_T2), "k", linewidth=1)
    ax.plot(x_side_T3, twisstntern.utils.T3(0, x_side_T3), "k", linewidth=1)
    ax.hlines(y=0, xmin=-0.5, xmax=0.5, color="k", linewidth=1)
    
    # Remove ticks
    ax.set_xticks([])
    ax.set_yticks([])
    
    # Draw grid lines using twisstntern functions
    for i in range(1, int(1 / alpha)):
        y = i * alpha
        
        # T1 lines (horizontal)
        ax.hlines(y=y * h, xmin=twisstntern.utils.T1_lim(y)[0], 
                 xmax=twisstntern.utils.T1_lim(y)[1], color="#C74375", linewidth=1)
        
        # T2 lines
        x2 = np.linspace(twisstntern.utils.T2_lim(y)[0], twisstntern.utils.T2_lim(y)[1], 100)
        ax.plot(x2, twisstntern.utils.T2(y, x2), color="#3182BD", linewidth=1)
        
        # T3 lines
        x3 = np.linspace(twisstntern.utils.T3_lim(y)[0], twisstntern.utils.T3_lim(y)[1], 100)
        ax.plot(x3, twisstntern.utils.T3(y, x3), color="#D4A017", linewidth=1)
    
    # Central vertical line
    ax.vlines(x=0, ymin=0, ymax=h, colors="#888888", ls=':')
    
    # Labels
    ax.text(-0.02, 0.88, 'T₁', size=12, color="#C74375", fontweight='bold')
    ax.text(0.54, -0.01, 'T₃', size=12, color="#D4A017", fontweight='bold')
    ax.text(-0.58, -0.01, 'T₂', size=12, color="#3182BD", fontweight='bold')
    
    # Remove spines
    for spine in ax.spines.values():
        spine.set_visible(False)


def create_heatmap_colormap():
    """Create consistent colormap for heatmaps with vibrant yellow ending."""
    # Create a custom colormap from dark purple to vibrant yellow
    colors = ['#1a0033', '#4a0080', '#8000ff', '#bf40ff', '#ff8000', '#ffff00']
    return sns.blend_palette(colors, as_cmap=True)


def plot_heatmap_data_twiss(ax, results_df, alpha, data_type='prop_data', title='Data', vmax=None):
    """
    Plot heatmap data using twisstntern coordinate functions.
    
    Args:
        ax: Matplotlib axis
        results_df: Results dataframe
        alpha: Grid granularity
        data_type: Column to plot ('prop_data', 'prop_model', 'prop_residual')
        title: Plot title
        vmax: Maximum value for colormap scaling
    """
    plot_ternary_base_twiss(ax, alpha)
    
    # Determine colormap and normalization
    if 'residual' in data_type:
        # Use diverging colormap for residuals
        cmap = sns.color_palette("RdBu_r", as_cmap=True)
        vmax_val = results_df[data_type].abs().max()
        vmin_val = -vmax_val
        norm = Normalize(vmin=vmin_val, vmax=vmax_val)
    else:
        # Use sequential colormap for data/model
        cmap = create_heatmap_colormap()
        if vmax is None:
            vmax_val = results_df[data_type].max()
        else:
            vmax_val = vmax
        vmin_val = 0
        norm = Normalize(vmin=vmin_val, vmax=vmax_val)
    
    # Plot filled triangles using twisstntern coordinate functions
    for idx, row in results_df.iterrows():
        value = row[data_type]
        
        # Get triangle coordinates
        trianglex, triangley, direction = twisstntern.utils.return_triangle_coord(
            row['T1_bounds'][0], row['T1_bounds'][1],
            row['T2_bounds'][0], row['T2_bounds'][1],
            row['T3_bounds'][0], row['T3_bounds'][1]
        )
        
        # Determine if we should plot this triangle
        should_plot = False
        
        if data_type in ['prop_data', 'prop_model']:
            # For data/model plots: only plot if there are points in this dataset
            count_key = 'count_data' if data_type == 'prop_data' else 'count_model'
            should_plot = row[count_key] > 0
        elif data_type == 'prop_residual':
            # For residuals: only plot if at least one dataset has points
            should_plot = (row['count_data'] > 0) or (row['count_model'] > 0)
        else:
            # Default behavior for other data types
            should_plot = not np.isnan(value) and abs(value) > 1e-10
        
        if should_plot:
            color = cmap(norm(value))
            ax.fill(trianglex, triangley, color=color, edgecolor='none', alpha=0.8)
    
    ax.set_title(title, fontsize=12, fontweight='bold')
    
    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax = inset_axes(ax, width="3%", height="80%", loc='center right',
                     bbox_to_anchor=(0.05, 0, 1, 1), bbox_transform=ax.transAxes, borderpad=1)
    cbar = plt.colorbar(sm, cax=cax)
    cbar.set_label('Proportion', fontsize=10)
    
    # Clean up spines
    sns.despine(ax=ax, top=True, right=True)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.8)


def plot_residual_histogram_enhanced(ax, results_df, statistics):
    """
    Plot residual histogram with enhanced statistics display.
    
    Args:
        ax: Matplotlib axis
        results_df: Results dataframe
        statistics: Statistics dictionary
    """
    # Only use residuals from meaningful triangles (where at least one dataset has points)
    meaningful_triangles = (results_df['count_data'] > 0) | (results_df['count_model'] > 0)
    residuals = results_df.loc[meaningful_triangles, 'prop_residual'].values
    
    # Remove outliers for better visualization
    q1, q3 = np.percentile(residuals, [25, 75])
    iqr = q3 - q1
    lower_bound = q1 - 1.5 * iqr
    upper_bound = q3 + 1.5 * iqr
    filtered_residuals = residuals[(residuals >= lower_bound) & (residuals <= upper_bound)]
    
    # Create histogram with modern styling
    n_bins = min(50, len(filtered_residuals) // 10) if len(filtered_residuals) > 10 else 10
    ax.hist(filtered_residuals, bins=n_bins, density=True, alpha=0.7, 
           color=sns.color_palette("bright")[2], edgecolor='white', linewidth=0.5)
    
    # Add normal distribution overlay
    if len(filtered_residuals) > 1:
        mu, sigma = np.mean(filtered_residuals), np.std(filtered_residuals)
        x = np.linspace(filtered_residuals.min(), filtered_residuals.max(), 100)
        normal_curve = stats.norm.pdf(x, mu, sigma)
        ax.plot(x, normal_curve, 'r-', linewidth=2, label='Normal fit', alpha=0.8)
    
        # Add vertical lines for mean and median
        ax.axvline(np.mean(filtered_residuals), color='red', linestyle='--', 
                  linewidth=2, alpha=0.7, label=f'Mean: {np.mean(filtered_residuals):.4f}')
        ax.axvline(np.median(filtered_residuals), color='blue', linestyle='--', 
                  linewidth=2, alpha=0.7, label=f'Median: {np.median(filtered_residuals):.4f}')
    
    ax.set_xlabel('Residuals (Data - Model)', fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.set_title('Residual Distribution', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9)
    
    # Add statistics text
    stats_text = (f'μ = {statistics["mean_residual"]:.4f}\n'
                  f'σ = {statistics["std_residual"]:.4f}\n'
                  f'Wasserstein (Eucl): {statistics["wasserstein_euclidean"]:.4f}\n'
                  f'Wasserstein (KL): {statistics["wasserstein_kl"]:.4f}')
    
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
           verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', 
           facecolor='lightblue', alpha=0.7))
    
    # Clean up spines
    sns.despine(ax=ax, top=True, right=True)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.8)


def plot_l2_distances_twiss(ax, results_df, alpha):
    """
    Plot L2 distances per triangle using twisstntern functions.
    
    Args:
        ax: Matplotlib axis
        results_df: Results dataframe
        alpha: Grid granularity
    """
    plot_ternary_base_twiss(ax, alpha)
    
    # Calculate L2 per triangle (square root of squared residuals)
    l2_per_triangle = np.sqrt(results_df['residual_squared'].values)
    
    # Use vibrant colormap for L2 distances
    vmax = np.max(l2_per_triangle) if len(l2_per_triangle) > 0 else 1
    vmin = 0
    cmap = sns.color_palette("magma", as_cmap=True)  # Vibrant purple-to-yellow
    norm = Normalize(vmin=vmin, vmax=vmax)
    
    # Plot filled triangles using twisstntern functions
    for idx, row in results_df.iterrows():
        # Only plot if at least one dataset has points and L2 distance > 0
        has_meaningful_data = (row['count_data'] > 0) or (row['count_model'] > 0)
        
        if has_meaningful_data and l2_per_triangle[idx] > 0:
            # Get triangle coordinates using twisstntern
            trianglex, triangley, direction = twisstntern.utils.return_triangle_coord(
                row['T1_bounds'][0], row['T1_bounds'][1],
                row['T2_bounds'][0], row['T2_bounds'][1],
                row['T3_bounds'][0], row['T3_bounds'][1]
            )
            
            color = cmap(norm(l2_per_triangle[idx]))
            ax.fill(trianglex, triangley, color=color, edgecolor='none')
    
    ax.set_title('L² Distance per Triangle', fontsize=12, fontweight='bold')
    
    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax = inset_axes(ax, width="3%", height="80%", loc='center right',
                     bbox_to_anchor=(0.05, 0, 1, 1), bbox_transform=ax.transAxes, borderpad=1)
    cbar = plt.colorbar(sm, cax=cax)
    cbar.set_label('L2 Distance', fontsize=10)


def create_enhanced_comparison_plot(data1, data2, results_df, statistics, alpha, output_path=None):
    """
    Create the main enhanced comparison plot with five subplots.
    
    Args:
        data1, data2: Input datasets
        results_df: Analysis results
        statistics: Statistics dictionary
        alpha: Grid granularity
        output_path: Optional path to save the plot
        
    Returns:
        matplotlib.figure.Figure: The created figure
    """
    fig = plt.figure(figsize=(18, 12))
    
    # Create layout: 3x4 grid with statistics box beneath residuals plot
    gs = fig.add_gridspec(3, 4, hspace=0.35, wspace=0.25, 
                         height_ratios=[1, 1, 0.6])  # Third row shorter for statistics
    ax1 = fig.add_subplot(gs[0, 0])  # Data
    ax2 = fig.add_subplot(gs[0, 1])  # Model  
    ax3 = fig.add_subplot(gs[0, 2])  # Residuals
    ax4 = fig.add_subplot(gs[1, 0])  # L2 distances
    ax5 = fig.add_subplot(gs[1, 1])  # Histogram
    # Statistics box will go in gs[2, 2] (beneath residuals)
    
    # Determine consistent scale for data and model plots
    max_prop = max(results_df['prop_data'].max(), results_df['prop_model'].max())
    
    # Plot 1: Data (reference/truth)
    plot_heatmap_data_twiss(ax1, results_df, alpha, 'prop_data', 'Data', vmax=max_prop)
    
    # Plot 2: Model
    plot_heatmap_data_twiss(ax2, results_df, alpha, 'prop_model', 'Model', vmax=max_prop)
    
    # Plot 3: Residuals
    plot_heatmap_data_twiss(ax3, results_df, alpha, 'prop_residual', 'Residuals')
    
    # Plot 4: L2 distances per triangle
    plot_l2_distances_twiss(ax4, results_df, alpha)
    
    # Plot 5: Enhanced residual histogram
    plot_residual_histogram_enhanced(ax5, results_df, statistics)
    
    # Add enhanced overall title with better blue color
    title_text = f'Enhanced Topology Weights Residual Analysis (α = {alpha})'
    fig.suptitle(title_text, fontsize=17, y=0.96, color='#2E5090', fontweight='bold')
    
    # Create statistics box directly beneath the residuals plot
    ax_stats = fig.add_subplot(gs[2, 2])  # Third row, third column (beneath residuals)
    ax_stats.axis('off')  # Hide axes
    
    # Create more compact statistics text for the box beneath residuals
    stats_text = ('Statistical Metrics\n' + '─'*18 + '\n'
                  f'L² distance: {statistics["L2_distance"]:.5f}\n'
                  f'χ² statistic: {statistics["chi2_statistic"]:.1f}\n'
                  f'Degrees freedom: {statistics["degrees_freedom"]}\n'
                  f'p-value: {statistics["p_value"]:.2e}\n'
                  f'Wasserstein (Eucl): {statistics["wasserstein_euclidean"]:.5f}\n'
                  f'Wasserstein (KL): {statistics["wasserstein_kl"]:.5f}\n\n'
                  'Data Summary\n' + '─'*12 + '\n'
                  f'Data points: {statistics["total_data_points"]:,}\n'
                  f'Model points: {statistics["total_model_points"]:,}\n'
                  f'Meaningful triangles: {statistics["num_meaningful_triangles"]}\n'
                  f'Empty triangles: {statistics["num_empty_triangles"]}')
    
    # Place the statistics box beneath residuals plot
    ax_stats.text(0.05, 0.95, stats_text, fontsize=10, verticalalignment='top',
                 horizontalalignment='left', transform=ax_stats.transAxes,
                 bbox=dict(boxstyle='round,pad=0.4', facecolor='#f8f9fa', 
                          alpha=0.95, edgecolor='#2E5090', linewidth=1.5),
                 family='monospace')
    
    # Use tight layout for automatic spacing
    plt.tight_layout(pad=1.8)
    
    # Fine-tune layout
    plt.subplots_adjust(top=0.93, bottom=0.05, left=0.05, right=0.98)
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Enhanced plot saved to: {output_path}")
    
    return fig


# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    """Main function to run the enhanced residual analysis."""
    parser = argparse.ArgumentParser(description='Enhanced Topology Weights Residual Analysis')
    parser.add_argument('--granularity', type=float, default=0.1,
                       help='Grid granularity (default: 0.1)')
    parser.add_argument('--output', type=str, default='enhanced_residual_results',
                       help='Output directory (default: enhanced_residual_results)')
    parser.add_argument('--data1', type=str, default='migration_topology_weights.csv',
                       help='First dataset filename (default: migration_topology_weights.csv)')
    parser.add_argument('--data2', type=str, default='NoMigration_topology_weights.csv',
                       help='Second dataset filename (default: NoMigration_topology_weights.csv)')
    parser.add_argument('--sample_size', type=int, default=None,
                       help='Sample size for large datasets (default: use all data)')
    
    args = parser.parse_args()
    
    # Get script directory
    script_dir = Path(__file__).parent
    
    # File paths
    data1_path = script_dir / args.data1
    data2_path = script_dir / args.data2
    
    # Check if files exist
    if not data1_path.exists():
        print(f"Error: Data file not found: {data1_path}")
        sys.exit(1)
    if not data2_path.exists():
        print(f"Error: Data file not found: {data2_path}")
        sys.exit(1)
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(exist_ok=True)
    
    print("=== ENHANCED Topology Weights Residual Analysis ===")
    print(f"Using twisstntern library functions and twisstCompareEx.py")
    print(f"Data 1 (reference): {data1_path}")
    print(f"Data 2 (model): {data2_path}")
    print(f"Granularity: {args.granularity}")
    print(f"Output directory: {output_dir}")
    if args.sample_size:
        print(f"Sample size: {args.sample_size}")
    print()
    
    # Load data using enhanced loader
    print("Loading datasets with twisstntern functions...")
    data1 = load_topology_data_twiss(data1_path, sample_size=args.sample_size)
    data2 = load_topology_data_twiss(data2_path, sample_size=args.sample_size)
    
    # Perform enhanced analysis
    results_df, statistics = perform_enhanced_grid_analysis(data1, data2, args.granularity)
    
    # Create enhanced visualization
    print("\nCreating enhanced visualizations...")
    output_plot = output_dir / f'enhanced_residual_analysis_granularity_{args.granularity}.png'
    fig = create_enhanced_comparison_plot(data1, data2, results_df, statistics, 
                                        args.granularity, output_plot)
    
    # Save detailed results
    results_csv = output_dir / f'enhanced_residual_results_granularity_{args.granularity}.csv'
    results_df.to_csv(results_csv, index=False)
    print(f"Detailed results saved to: {results_csv}")
    
    # Save enhanced summary statistics
    stats_file = output_dir / f'enhanced_summary_statistics_granularity_{args.granularity}.txt'
    with open(stats_file, 'w') as f:
        f.write("Enhanced Topology Weights Residual Analysis Summary\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Data 1: {data1_path}\n")
        f.write(f"Data 2: {data2_path}\n")
        f.write(f"Granularity: {args.granularity}\n")
        f.write(f"Sample size: {args.sample_size if args.sample_size else 'All data'}\n\n")
        
        f.write("COMPREHENSIVE METRICS:\n")
        f.write("-" * 30 + "\n")
        for key, value in statistics.items():
            if isinstance(value, float):
                f.write(f"{key}: {value:.6e}\n")
            else:
                f.write(f"{key}: {value}\n")
    
    print(f"Enhanced summary statistics saved to: {stats_file}")
    
    # Show plot
    plt.show()
    
    print("\n" + "="*60)
    print("ENHANCED ANALYSIS COMPLETE!")
    print("Key enhancements:")
    print("• Used twisstntern.utils.dump_data() for data loading")
    print("• Used twisstntern coordinate functions for plotting")
    print("• Integrated twisstCompareEx.py functions for analysis")
    print("• Added Wasserstein distance calculations")
    print("• Enhanced statistical reporting")
    print("="*60)


if __name__ == "__main__":
    main()
