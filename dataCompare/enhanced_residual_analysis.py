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
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.colors import ListedColormap, Normalize, LinearSegmentedColormap
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
from matplotlib.patches import Polygon
import matplotlib.patches as mpatches

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

# Add import for compute_metrics
from dataCompare.metrics import compute_metrics

# ============================================================================
# COLORMAP CONFIGURATION - EASY TO CHANGE!
# ============================================================================
# Colormap for Data and Model plots (sequential colormaps work best)
DATA_MODEL_COLORMAP = "viridis"  # Options: "viridis", "plasma", "inferno", "magma", "cividis", "mako", "rocket", "flare"

# Colormap for L2 Distance plot (sequential colormaps work best)  
L2_COLORMAP = "mako_r"  # Options: "mako_r", "viridis", "plasma", "inferno", "magma", "cividis", "rocket", "flare"

# Colormap for Residuals plot (diverging colormaps work best)
RESIDUALS_COLORMAP = "RdBu_r"  # Options: "RdBu_r", "seismic", "coolwarm", "bwr", "PiYG", "PRGn"

# Histogram color
HISTOGRAM_COLOR = "#ffb347"  # Orange color for residuals histogram

# KDE curve color  
KDE_COLOR = "#22223b"  # Dark blue for KDE overlay


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
    """Wasserstein distance with Euclidean cost matrix (using Sinkhorn for speed)."""
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
    """Simplex Wasserstein distance with KL cost matrix (using Sinkhorn for speed)."""
    if isinstance(data1, pd.DataFrame):
        data1 = data1[['T1', 'T2', 'T3']].values
    if isinstance(data2, pd.DataFrame):
        data2 = data2[['T1', 'T2', 'T3']].values
    data1_plus_eps = data1 + 1e-8
    data1_plus_eps = data1_plus_eps / data1_plus_eps.sum(1, keepdims=True)
    data2_plus_eps = data2 + 1e-8
    data2_plus_eps = data2_plus_eps / data2_plus_eps.sum(1, keepdims=True)
    kl_pq = (data1_plus_eps[None, :, :] * (np.log(data1_plus_eps[None, :, :]) - np.log(data2_plus_eps[:, None, :]))).sum(-1)
    kl_qp = (data2_plus_eps[None, :, :] * (np.log(data2_plus_eps[None, :, :]) - np.log(data1_plus_eps[:, None, :]))).sum(-1)
    kl_symm = (kl_pq + kl_qp.T)
    # Use Sinkhorn for speed (reg=0.01). Switch back to ot.emd2 for exact.
    dist = ot.sinkhorn2(np.ones(data2_plus_eps.shape[0])/data2_plus_eps.shape[0], 
                       np.ones(data1_plus_eps.shape[0])/data1_plus_eps.shape[0],
                       kl_symm, reg=0.01)
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
    print(f"Performing enhanced grid analysis with α = {alpha}")
    triangles = create_triangular_grid_twiss(alpha)
    # Count points in each triangle using enhanced counting (run only once)
    results_df = []
    for triangle in triangles:
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
        prop_data = count_data / len(data1)
        prop_model = count_model / len(data2)
        prop_residual = prop_data - prop_model
        residual_squared = prop_residual ** 2
        count_residual = count_data - count_model
        results_df.append({
            'T1_bounds': triangle['T1'],
            'T2_bounds': triangle['T2'],
            'T3_bounds': triangle['T3'],
            'count_data': count_data,
            'count_model': count_model,
            'prop_data': prop_data,
            'prop_model': prop_model,
            'prop_residual': prop_residual,
            'residual_squared': residual_squared,
            'count_residual': count_residual
        })
    results_df = pd.DataFrame(results_df)
    meaningful_triangles = (results_df['count_data'] > 0) | (results_df['count_model'] > 0)
    filtered_results = results_df[meaningful_triangles]
    # Compute metrics from results_df (no double computation)
    # L2 distance
    prop1 = results_df['prop_data']
    prop2 = results_df['prop_model']
    l2 = np.linalg.norm(prop1 - prop2)
    # Chi-square
    counts1 = results_df['count_data']
    counts2 = results_df['count_model']
    with np.errstate(divide='ignore', invalid='ignore'):
        chi2_stat = np.nansum((counts1 - counts2) ** 2 / (counts1 + counts2 + 1e-8))
    dof = (prop1 != 0).sum() - 1
    from scipy.stats import chi2
    p_value = 1 - chi2.cdf(chi2_stat, dof)
    if len(filtered_results) > 0:
        mean_residual = filtered_results['prop_residual'].mean()
        std_residual = filtered_results['prop_residual'].std()
        max_residual = filtered_results['prop_residual'].abs().max()
        mean_squared_residual = filtered_results['residual_squared'].mean()
    else:
        mean_residual = 0.0
        std_residual = 0.0
        max_residual = 0.0
        mean_squared_residual = 0.0
    statistics = {
        'L2_distance': l2,
        'chi2_statistic': chi2_stat,
        'p_value': p_value,
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
    print("\n" + "="*70)
    print("COMPREHENSIVE GRID-BASED METRICS (no Wasserstein)")
    print("="*70)
    print(f"L² Distance:              {statistics['L2_distance']:.6f}")
    print(f"χ² Statistic:             {statistics['chi2_statistic']:.6f}")
    print(f"χ² p-value:               {statistics['p_value']:.6e}")
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

def draw_empty_triangle(ax, trianglex, triangley):
    triangle_coords = list(zip(trianglex, triangley))
    empty_triangle = Polygon(
        triangle_coords,
        closed=True,
        facecolor='white',
        edgecolor='grey',
        hatch='///',
        linewidth=0.5
    )
    ax.add_patch(empty_triangle)


def plot_ternary_base_twiss(ax, alpha):
    """
    Create ternary plot base using twisstntern functions, with all gridlines in grey.
    Args:
        ax: Matplotlib axis
        alpha: Grid granularity
    """
    h = math.sqrt(3) / 2
    # Plot triangle outline
    x_side_T2 = np.linspace(0, 0.5, 100)
    x_side_T3 = np.linspace(-0.5, 0, 100)
    ax.plot(x_side_T2, twisstntern.utils.T2(0, x_side_T2), color="grey", linewidth=1)
    ax.plot(x_side_T3, twisstntern.utils.T3(0, x_side_T3), color="grey", linewidth=1)
    ax.hlines(y=0, xmin=-0.5, xmax=0.5, color="grey", linewidth=1)
    # Remove ticks
    ax.set_xticks([])
    ax.set_yticks([])
    # Draw grid lines using twisstntern functions (all grey)
    for i in range(1, int(1 / alpha)):
        y = i * alpha
        # T1 lines (horizontal)
        ax.hlines(y=y * h, xmin=twisstntern.utils.T1_lim(y)[0], 
                 xmax=twisstntern.utils.T1_lim(y)[1], color="grey", linewidth=1)
        # T2 lines
        x2 = np.linspace(twisstntern.utils.T2_lim(y)[0], twisstntern.utils.T2_lim(y)[1], 100)
        ax.plot(x2, twisstntern.utils.T2(y, x2), color="grey", linewidth=1)
        # T3 lines
        x3 = np.linspace(twisstntern.utils.T3_lim(y)[0], twisstntern.utils.T3_lim(y)[1], 100)
        ax.plot(x3, twisstntern.utils.T3(y, x3), color="grey", linewidth=1)
    # Central vertical line
    ax.vlines(x=0, ymin=0, ymax=h, colors="grey", ls=':')
    # Labels (closer to triangle)
    ax.text(0.0, -0.07, 'T₁', size=13, color="black", fontweight='bold', ha='center', va='top')
    ax.text(0.52, 0.03, 'T₂', size=13, color="black", fontweight='bold', ha='left', va='center')
    ax.text(-0.52, 0.03, 'T₃', size=13, color="black", fontweight='bold', ha='right', va='center')
    # Remove spines
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_ternary_with_colorbar(ax, results_df, alpha, data_type, title, vmin, vmax, cmap, cbar_label, cbar_title=None):
    plot_heatmap_data_twiss(ax, results_df, alpha, data_type, title, vmax=vmax, vmin=vmin, override_cmap=cmap)
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    import matplotlib.pyplot as plt
    norm = Normalize(vmin=vmin, vmax=vmax)
    cax = inset_axes(ax, width="3%", height="60%", loc='right', borderpad=2)
    cb = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
    cb.ax.set_title(cbar_label, fontsize=12, fontweight='normal', pad=10)
    return cb


def plot_heatmap_data_twiss(ax, results_df, alpha, data_type='prop_data', title='Data', vmax=None, vmin=None, override_cmap=None):
    plot_ternary_base_twiss(ax, alpha)
    import matplotlib.pyplot as plt
    if data_type == 'count_residual':
        cmap = sns.color_palette(RESIDUALS_COLORMAP, as_cmap=True) if override_cmap is None else override_cmap
        vmax_val = results_df[data_type].abs().max()
        vmin_val = -vmax_val
        norm = Normalize(vmin=vmin_val, vmax=vmax_val)
        values = results_df[data_type]
        empty_logic = lambda row: row['count_data'] == 0 and row['count_model'] == 0
    elif data_type == 'count_data':
        cmap = plt.get_cmap(DATA_MODEL_COLORMAP) if override_cmap is None else override_cmap
        vmax_val = vmax if vmax is not None else results_df[data_type].max()
        vmin_val = vmin if vmin is not None else results_df[data_type].min()
        norm = Normalize(vmin=vmin_val, vmax=vmax_val)
        values = results_df[data_type]
        empty_logic = lambda row: row['count_data'] == 0
    elif data_type == 'count_model':
        cmap = plt.get_cmap(DATA_MODEL_COLORMAP) if override_cmap is None else override_cmap
        vmax_val = vmax if vmax is not None else results_df[data_type].max()
        vmin_val = vmin if vmin is not None else results_df[data_type].min()
        norm = Normalize(vmin=vmin_val, vmax=vmax_val)
        values = results_df[data_type]
        empty_logic = lambda row: row['count_model'] == 0
    elif data_type == 'l2':
        cmap = override_cmap if override_cmap is not None else plt.get_cmap(L2_COLORMAP)
        l2_per_triangle = np.sqrt(results_df['residual_squared'].values)
        vmax_val = vmax if vmax is not None else np.max(l2_per_triangle)
        vmin_val = vmin if vmin is not None else 0
        norm = Normalize(vmin=vmin_val, vmax=vmax_val)
        values = l2_per_triangle
        empty_logic = lambda row: (row['count_data'] == 0 and row['count_model'] == 0)
    elif 'residual' in data_type:
        cmap = sns.color_palette(RESIDUALS_COLORMAP, as_cmap=True) if override_cmap is None else override_cmap
        vmax_val = results_df[data_type].abs().max()
        vmin_val = -vmax_val
        norm = Normalize(vmin=vmin_val, vmax=vmax_val)
        values = results_df[data_type]
        empty_logic = lambda row: row['count_data'] == 0 and row['count_model'] == 0
    else:
        if data_type == 'prop_data':
            data_type = 'count_data'
            title = title.replace('Proportion', 'Count')
        elif data_type == 'prop_model':
            data_type = 'count_model'
            title = title.replace('Proportion', 'Count')
        cmap = plt.get_cmap(DATA_MODEL_COLORMAP) if override_cmap is None else override_cmap
        vmax_val = vmax if vmax is not None else results_df[data_type].max()
        vmin_val = vmin if vmin is not None else results_df[data_type].min()
        norm = Normalize(vmin=vmin_val, vmax=vmax_val)
        values = results_df[data_type]
        empty_logic = lambda row: row['count_data'] == 0
    for idx, row in results_df.iterrows():
        value = values[idx]
        trianglex, triangley, direction = twisstntern.utils.return_triangle_coord(
            row['T1_bounds'][0], row['T1_bounds'][1],
            row['T2_bounds'][0], row['T2_bounds'][1],
            row['T3_bounds'][0], row['T3_bounds'][1]
        )
        if empty_logic(row):
            draw_empty_triangle(ax, trianglex, triangley)
        else:
            color = cmap(norm(value))
            ax.fill(trianglex, triangley, color=color, edgecolor='none', alpha=0.8)
    ax.set_title(title, fontsize=12, fontweight='bold')
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
    stats_text = (
        f'μ = {statistics["mean_residual"]:.4f}\n'
        f'σ = {statistics["std_residual"]:.4f}\n'
        f'Wasserstein (Eucl): {statistics["wasserstein_euclidean"]:.4f}\n'
        f'Wasserstein (KL): {statistics["wasserstein_kl"]:.4f}'
    )
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=9,
           verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', 
           facecolor='lightblue', alpha=0.7))
    
    # Clean up spines
    sns.despine(ax=ax, top=True, right=True)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.8)


def plot_l2_distances_twiss(ax, results_df, alpha):
    import matplotlib.pyplot as plt
    plot_ternary_base_twiss(ax, alpha)
    l2_per_triangle = np.sqrt(results_df['residual_squared'].values)
    vmax = np.max(l2_per_triangle) if len(l2_per_triangle) > 0 else 1
    vmin = 0
    cmap = plt.get_cmap("viridis")
    norm = Normalize(vmin=vmin, vmax=vmax)
    for idx, row in results_df.iterrows():
        has_meaningful_data = (row['count_data'] > 0) or (row['count_model'] > 0)
        trianglex, triangley, direction = twisstntern.utils.return_triangle_coord(
            row['T1_bounds'][0], row['T1_bounds'][1],
            row['T2_bounds'][0], row['T2_bounds'][1],
            row['T3_bounds'][0], row['T3_bounds'][1]
        )
        if not has_meaningful_data:
            draw_empty_triangle(ax, trianglex, triangley)
        elif l2_per_triangle[idx] > 0:
            color = cmap(norm(l2_per_triangle[idx]))
            ax.fill(trianglex, triangley, color=color, edgecolor='none')
    ax.set_title('L² Distance per Triangle', fontsize=12, fontweight='bold')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    # No colorbar here


def create_enhanced_comparison_plot(data1, data2, results_df, statistics, alpha, output_path=None):
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    import seaborn as sns
    import numpy as np
    fig = plt.figure(figsize=(20, 10))
    gs = gridspec.GridSpec(2, 3, height_ratios=[1, 1.1], hspace=0.25, wspace=0.25)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = fig.add_subplot(gs[1, 0])
    ax5 = fig.add_subplot(gs[1, 1:])
    min_count = min(results_df['count_data'].min(), results_df['count_model'].min())
    max_count = max(results_df['count_data'].max(), results_df['count_model'].max())
    cmap_seq = plt.get_cmap(DATA_MODEL_COLORMAP)
    # Data plot
    plot_ternary_with_colorbar(ax1, results_df, alpha, 'count_data', 'Data', min_count, max_count, cmap_seq, "Count")
    # Model plot
    plot_ternary_with_colorbar(ax2, results_df, alpha, 'count_model', 'Model', min_count, max_count, cmap_seq, "Count")
    # Residuals plot
    plot_heatmap_data_twiss(ax3, results_df, alpha, 'count_residual', 'Residuals')
    vmax_resid = results_df['count_residual'].abs().max()
    norm_resid = Normalize(vmin=-vmax_resid, vmax=vmax_resid)
    cmap_resid = sns.color_palette(RESIDUALS_COLORMAP, as_cmap=True)
    cax3 = inset_axes(ax3, width="3%", height="60%", loc='right', borderpad=2)
    cb3 = plt.colorbar(cm.ScalarMappable(norm=norm_resid, cmap=cmap_resid), cax=cax3)
    cb3.ax.set_title("Count", fontsize=12, fontweight='normal', pad=10)
    # L2 plot
    cmap_l2 = plt.get_cmap(L2_COLORMAP)
    l2_per_triangle = np.sqrt(results_df['residual_squared'].values)
    vmax_l2 = np.max(l2_per_triangle) if len(l2_per_triangle) > 0 else 1
    plot_ternary_with_colorbar(ax4, results_df, alpha, 'l2', 'L2 Distance', 0, vmax_l2, cmap_l2, "L2 Distance")
    # Residuals histogram
    meaningful_triangles = (results_df['count_data'] > 0) | (results_df['count_model'] > 0)
    residuals = results_df.loc[meaningful_triangles, 'count_residual'].values
    # Histogram in configured color
    ax5.hist(residuals, bins=30, color=HISTOGRAM_COLOR, alpha=0.7, edgecolor='white', density=True)
    # Overlay KDE
    from scipy.stats import gaussian_kde
    x_grid = np.linspace(residuals.min(), residuals.max(), 200)
    kde = gaussian_kde(residuals)
    ax5.plot(x_grid, kde(x_grid), color=KDE_COLOR, linewidth=2, label="KDE fit")
    ax5.set_title('Residuals Distribution (Count)', fontsize=13, fontweight='bold')
    ax5.set_xlabel('Count Residual (Data - Model)', fontsize=12)
    ax5.set_ylabel('Density', fontsize=12)
    stats_text = (
        f"L2: {statistics['L2_distance']:.4f}\n"
        f"Chi2: {statistics['chi2_statistic']:.1f}\n"
        f"p-value: {statistics['p_value']:.2e}\n"
        f"Mean Residual: {statistics['mean_residual']:.4f}\n"
        f"Std Residual: {statistics['std_residual']:.4f}\n"
        f"Max |Residual|: {statistics['max_residual']:.4f}\n"
        f"Data Points: {statistics['total_data_points']}\n"
        f"Model Points: {statistics['total_model_points']}\n"
        f"Triangles: {statistics['num_triangles']}\n"
        f"Meaningful: {statistics['num_meaningful_triangles']}\n"
        f"Empty: {statistics['num_empty_triangles']}\n"
        f"DoF: {statistics['degrees_freedom']}"
    )
    ax5.text(0.02, 0.98, stats_text, transform=ax5.transAxes, fontsize=11, va='top', ha='left',
             bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', alpha=0.8, edgecolor='gray'))
    ax5.legend(loc='upper right', fontsize=11)
    sns.despine(ax=ax5)
    fig.suptitle(f'Residual Analysis (α = {alpha})', fontsize=18, y=0.98, color='#2E5090', fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.96])
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
    parser.add_argument('--skip_wasserstein', action='store_true',
                       help='Skip Wasserstein calculations for faster execution')
    
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
    
    # Timing: Data loading
    t0 = time.time()
    print("Loading datasets with twisstntern functions...")
    data1 = load_topology_data_twiss(data1_path, sample_size=args.sample_size)
    data2 = load_topology_data_twiss(data2_path, sample_size=args.sample_size)
    t1 = time.time()
    print(f"Data loading took {t1-t0:.2f} seconds.")
    
    # Timing: Grid analysis
    t2 = time.time()
    results_df, statistics = perform_enhanced_grid_analysis(data1, data2, args.granularity)
    t3 = time.time()
    print(f"Grid analysis took {t3-t2:.2f} seconds.")

    # Import and print metrics from metrics.py
    if not args.skip_wasserstein:
        print("\nImporting and computing full metrics (including Wasserstein) from metrics.py...")
        t4 = time.time()
        metrics = compute_metrics(data1, data2, args.granularity)
        t5 = time.time()
        print(f"Metrics computation took {t5-t4:.2f} seconds.")
        print(f"L2_distance: {metrics['l2']:.8f}")
        print(f"chi2_statistic: {metrics['chi2_stat']:.8f}")
        print(f"p_value: {metrics['p_value']:.8e}")
        print(f"wasserstein_euclidean: {metrics['wasserstein_euclidean']:.8f}")
        print(f"wasserstein_kl: {metrics['wasserstein_kl']:.8f}")
        
        # Add Wasserstein metrics to statistics for plotting
        statistics['wasserstein_euclidean'] = metrics['wasserstein_euclidean']
        statistics['wasserstein_kl'] = metrics['wasserstein_kl']
    else:
        print("\nSkipping Wasserstein calculations for faster execution.")
        print("Grid-based metrics only:")
        print(f"L2_distance: {statistics['L2_distance']:.8f}")
        print(f"chi2_statistic: {statistics['chi2_statistic']:.8f}")
        print(f"p_value: {statistics['p_value']:.8e}")
    
    # Timing: Plotting
    t4 = time.time()
    print("\nCreating enhanced visualizations...")
    output_plot = output_dir / f'enhanced_residual_analysis_granularity_{args.granularity}.png'
    fig = create_enhanced_comparison_plot(data1, data2, results_df, statistics, 
                                        args.granularity, output_plot)
    t5 = time.time()
    print(f"Plotting took {t5-t4:.2f} seconds.")
    
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
