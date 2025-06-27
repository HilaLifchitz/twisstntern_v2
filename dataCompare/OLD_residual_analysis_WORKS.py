#!/usr/bin/env python3
"""
Residual Analysis for Topology Weights Comparison

This script performs a residual analysis comparing two topology weights datasets
using triangular gridding and statistical metrics. It creates visualizations
similar to the provided example with data, model, residuals, and residual distribution.

Usage:
    python residual_analysis.py [--granularity 0.1] [--output results/]
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
plt.rcParams.update({
    "text.usetex": False,  # Set to True if you have LaTeX installed
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


def load_topology_data(filepath):
    """
    Load and process topology weights CSV file.
    
    Args:
        filepath (str): Path to the CSV file
        
    Returns:
        pd.DataFrame: Processed data with T1, T2, T3 columns
    """
    print(f"Loading data from: {filepath}")
    
    # Read the CSV file
    data = pd.read_csv(filepath)
    
    # Skip the topology definitions row (row 1) and use only numeric data
    if len(data) > 1:
        # Check if first row contains topology strings
        first_row = data.iloc[0]
        if any(isinstance(val, str) and '(' in str(val) for val in first_row):
            data = data.iloc[1:]  # Skip topology definition row
    
    # Ensure we have T1, T2, T3 columns
    if 'T1' in data.columns and 'T2' in data.columns and 'T3' in data.columns:
        data = data[['T1', 'T2', 'T3']]
    else:
        # Use first three columns
        data = data.iloc[:, :3]
        data.columns = ['T1', 'T2', 'T3']
    
    # Convert to numeric and handle any issues
    for col in ['T1', 'T2', 'T3']:
        data[col] = pd.to_numeric(data[col], errors='coerce')
    
    # Remove any rows with NaN values
    data = data.dropna()
    
    # Normalize each row to sum to 1 (just in case)
    row_sums = data.sum(axis=1)
    data = data.div(row_sums, axis=0)
    
    # Remove points where T2 == T3 (on y-axis)
    data = data[data['T2'] != data['T3']]
    
    print(f"Loaded {len(data)} data points after processing")
    return data


def create_triangular_grid(alpha):
    """
    Create a triangular grid for binning data points.
    
    Args:
        alpha (float): Grid granularity
        
    Returns:
        list: List of dictionaries containing grid cell boundaries
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


def count_points_in_triangle(data, triangle):
    """
    Count data points within a triangular grid cell.
    
    Args:
        data (pd.DataFrame): Data with T1, T2, T3 columns
        triangle (dict): Triangle boundaries
        
    Returns:
        int: Number of points in the triangle
    """
    T1_bounds = triangle['T1']
    T2_bounds = triangle['T2']
    T3_bounds = triangle['T3']
    
    # Handle boundary conditions to avoid double-counting
    if T1_bounds[0] == 0:
        condition_T1 = (data['T1'] >= T1_bounds[0]) & (data['T1'] <= T1_bounds[1])
    else:
        condition_T1 = (data['T1'] > T1_bounds[0]) & (data['T1'] <= T1_bounds[1])
    
    if T2_bounds[0] == 0:
        condition_T2 = (data['T2'] >= T2_bounds[0]) & (data['T2'] <= T2_bounds[1])
    else:
        condition_T2 = (data['T2'] > T2_bounds[0]) & (data['T2'] <= T2_bounds[1])
    
    if T3_bounds[0] == 0:
        condition_T3 = (data['T3'] >= T3_bounds[0]) & (data['T3'] <= T3_bounds[1])
    else:
        condition_T3 = (data['T3'] > T3_bounds[0]) & (data['T3'] <= T3_bounds[1])
    
    return len(data[condition_T1 & condition_T2 & condition_T3])


def perform_grid_analysis(data1, data2, alpha):
    """
    Perform grid-based analysis comparing two datasets.
    
    Args:
        data1 (pd.DataFrame): First dataset (reference/truth)
        data2 (pd.DataFrame): Second dataset (model)
        alpha (float): Grid granularity
        
    Returns:
        tuple: (results_df, statistics)
    """
    print(f"Performing grid analysis with granularity α = {alpha}")
    
    # Create triangular grid
    triangles = create_triangular_grid(alpha)
    print(f"Created {len(triangles)} triangular grid cells")
    
    # Count points in each triangle for both datasets
    results = []
    for i, triangle in enumerate(triangles):
        count1 = count_points_in_triangle(data1, triangle)
        count2 = count_points_in_triangle(data2, triangle)
        
        results.append({
            'triangle_id': i,
            'T1_bounds': triangle['T1'],
            'T2_bounds': triangle['T2'],
            'T3_bounds': triangle['T3'],
            'count_data': count1,
            'count_model': count2,
            'residual': count1 - count2,
            'residual_squared': (count1 - count2) ** 2
        })
    
    results_df = pd.DataFrame(results)
    
    # Calculate statistics
    total_data = len(data1)
    total_model = len(data2)
    
    # Normalize counts to proportions
    results_df['prop_data'] = results_df['count_data'] / total_data
    results_df['prop_model'] = results_df['count_model'] / total_model
    results_df['prop_residual'] = results_df['prop_data'] - results_df['prop_model']
    
    # Calculate L2 distance
    L2_distance = np.sqrt(results_df['residual_squared'].sum() / len(results_df))
    
    # Calculate chi-square statistic (where expected > 0)
    valid_triangles = results_df[results_df['count_data'] > 0]
    if len(valid_triangles) > 0:
        chi2_statistic = (valid_triangles['residual_squared'] / valid_triangles['count_data']).sum()
        degrees_freedom = len(valid_triangles) - 1
        p_value = chi2.sf(chi2_statistic, degrees_freedom) if degrees_freedom > 0 else 1.0
    else:
        chi2_statistic = 0
        degrees_freedom = 0
        p_value = 1.0
    
    # Calculate residual statistics
    residuals = results_df['prop_residual'].values
    residual_mean = np.mean(residuals)
    residual_std = np.std(residuals)
    residual_min = np.min(residuals)
    residual_max = np.max(residuals)
    
    statistics = {
        'L2_distance': L2_distance,
        'chi2_statistic': chi2_statistic,
        'degrees_freedom': degrees_freedom,
        'p_value': p_value,
        'residual_mean': residual_mean,
        'residual_std': residual_std,
        'residual_min': residual_min,
        'residual_max': residual_max,
        'total_data_points': total_data,
        'total_model_points': total_model,
        'num_triangles': len(triangles)
    }
    
    print(f"Analysis complete:")
    print(f"  L2 distance: {L2_distance:.6f}")
    print(f"  Chi-square statistic: {chi2_statistic:.6f}")
    print(f"  Degrees of freedom: {degrees_freedom}")
    print(f"  P-value: {p_value:.6e}")
    print(f"  Residual mean: {residual_mean:.6f}")
    print(f"  Residual std: {residual_std:.6f}")
    
    return results_df, statistics


def plot_ternary_base(ax, alpha):
    """
    Plot the basic ternary triangle with grid lines.
    
    Args:
        ax: Matplotlib axis
        alpha (float): Grid granularity
    """
    h = math.sqrt(3) / 2
    
    # Draw main triangle edges
    x_side_T2 = np.linspace(0, 0.5, 100)
    x_side_T3 = np.linspace(-0.5, 0, 100)
    
    ax.plot(x_side_T2, twisstntern.utils.T2(0, x_side_T2), "k", linewidth=1)
    ax.plot(x_side_T3, twisstntern.utils.T3(0, x_side_T3), "k", linewidth=1)
    ax.hlines(y=0, xmin=-0.5, xmax=0.5, color="k", linewidth=1)
    
    # Draw grid lines with modern vibrant colors
    grid_color = sns.color_palette("muted")[7]  # Sophisticated gray from muted palette
    for i in range(1, int(1 / alpha)):
        y = i * alpha
        ax.hlines(y=y * h, xmin=twisstntern.utils.T1_lim(y)[0], xmax=twisstntern.utils.T1_lim(y)[1], 
                 color=grid_color, linewidth=0.8, alpha=0.6)
        
        x2 = np.linspace(twisstntern.utils.T2_lim(y)[0], twisstntern.utils.T2_lim(y)[1], 100)
        ax.plot(x2, twisstntern.utils.T2(y, x2), color=grid_color, linewidth=0.8, alpha=0.6)
        
        x3 = np.linspace(twisstntern.utils.T3_lim(y)[0], twisstntern.utils.T3_lim(y)[1], 100)
        ax.plot(x3, twisstntern.utils.T3(y, x3), color=grid_color, linewidth=0.8, alpha=0.6)
    
    # Vertical line through center
    ax.vlines(x=0, ymin=0, ymax=h, colors=grid_color, ls=':', linewidth=1.0, alpha=0.7)
    
    # Labels with vibrant colors
    label_colors = sns.color_palette("bright", 3)  # Use bright palette for vibrant labels
    ax.text(-0.02, 0.88, 'T1', size=13, fontweight='bold', color=label_colors[0])
    ax.text(0.54, -0.02, 'T3', size=13, fontweight='bold', color=label_colors[1]) 
    ax.text(-0.58, -0.02, 'T2', size=13, fontweight='bold', color=label_colors[2])
    
    # Remove ticks and spines
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def create_heatmap_colormap():
    """Create a seaborn-like colormap for heatmap visualization."""
    # Use seaborn's RdBu_r colormap which is similar to the original but more refined
    return plt.cm.RdBu_r


def plot_heatmap_data(ax, results_df, alpha, data_type='prop_data', title='Data', vmax=None):
    """
    Plot triangular heatmap of data.
    
    Args:
        ax: Matplotlib axis
        results_df: Results dataframe
        alpha: Grid granularity
        data_type: Column to plot ('prop_data', 'prop_model', 'prop_residual')
        title: Plot title
        vmax: Maximum value for colormap
    """
    plot_ternary_base(ax, alpha)
    
    # Get data values and set seaborn-like color schemes
    if data_type == 'prop_residual':
        values = results_df[data_type].values
        # Use vibrant diverging colormap for residuals
        vmin, vmax = -np.max(np.abs(values)), np.max(np.abs(values))
        cmap = sns.diverging_palette(145, 300, s=60, as_cmap=True)  # Vibrant teal to magenta
        norm = Normalize(vmin=vmin, vmax=vmax)
    else:
        values = results_df[data_type].values
        if vmax is None:
            vmax = np.max(values)
        vmin = 0
        # Use vibrant sequential colormap
        cmap = sns.color_palette("viridis", as_cmap=True)  # Vibrant green-to-yellow
        norm = Normalize(vmin=vmin, vmax=vmax)
    
    # Plot filled triangles
    for idx, row in results_df.iterrows():
        if values[idx] > 0 or data_type == 'prop_residual':
            # Get triangle coordinates
            trianglex, triangley, direction = twisstntern.utils.return_triangle_coord(
                row['T1_bounds'][0], row['T1_bounds'][1],
                row['T2_bounds'][0], row['T2_bounds'][1],
                row['T3_bounds'][0], row['T3_bounds'][1]
            )
            
            color = cmap(norm(values[idx]))
            ax.fill(trianglex, triangley, color=color, edgecolor='none')
    
    ax.set_title(title, fontsize=12, fontweight='bold')
    
    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax = inset_axes(ax, width="3%", height="80%", loc='center right',
                     bbox_to_anchor=(0.05, 0, 1, 1), bbox_transform=ax.transAxes, borderpad=1)
    cbar = plt.colorbar(sm, cax=cax)
    
    if data_type == 'prop_residual':
        cbar.set_label('Residual', fontsize=10)
    else:
        cbar.set_label('Proportion', fontsize=10)


def plot_residual_histogram(ax, results_df, statistics):
    """
    Plot seaborn-style histogram of residuals.
    
    Args:
        ax: Matplotlib axis
        results_df: Results dataframe
        statistics: Statistics dictionary
    """
    residuals = results_df['prop_residual'].values
    
    # Create vibrant histogram
    vibrant_colors = sns.color_palette("bright", 8)
    n_bins = min(50, int(np.sqrt(len(residuals))))
    
    # Main histogram with vibrant styling
    counts, bins, patches = ax.hist(residuals, bins=n_bins, density=True, 
                                   alpha=0.8, color=vibrant_colors[4], 
                                   edgecolor='white', linewidth=1.2)
    
    # Overlay normal distribution with vibrant color
    x = np.linspace(residuals.min(), residuals.max(), 100)
    normal_dist = stats.norm.pdf(x, statistics['residual_mean'], statistics['residual_std'])
    ax.plot(x, normal_dist, color=vibrant_colors[3], linewidth=3, 
           label='Normal fit', alpha=0.9)
    
    # Add vertical line at mean with vibrant color
    ax.axvline(statistics['residual_mean'], color=vibrant_colors[1], 
              linestyle='--', linewidth=2, alpha=0.8,
              label=f'Mean: {statistics["residual_mean"]:.4f}')
    
    # Seaborn-style labels and formatting
    ax.set_xlabel('Residual', fontsize=12, fontweight='medium')
    ax.set_ylabel('Density', fontsize=12, fontweight='medium')
    ax.set_title('Residual Distribution', fontsize=13, fontweight='bold', pad=15)
    
    # Seaborn-style legend
    legend = ax.legend(frameon=True, fancybox=True, shadow=True, 
                      framealpha=0.9, facecolor='white', edgecolor='gray')
    legend.get_frame().set_linewidth(0.8)
    
    # Clean up spines for seaborn look
    sns.despine(ax=ax, top=True, right=True)
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.8)


def plot_l2_distances(ax, results_df, alpha):
    """
    Plot L2 distances per triangle.
    
    Args:
        ax: Matplotlib axis
        results_df: Results dataframe
        alpha: Grid granularity
    """
    plot_ternary_base(ax, alpha)
    
    # Calculate L2 per triangle (square root of squared residuals)
    l2_per_triangle = np.sqrt(results_df['residual_squared'].values)
    
    # Use vibrant colormap for L2 distances
    vmax = np.max(l2_per_triangle)
    vmin = 0
    cmap = sns.color_palette("magma", as_cmap=True)  # Vibrant purple-to-yellow
    norm = Normalize(vmin=vmin, vmax=vmax)
    
    # Plot filled triangles
    for idx, row in results_df.iterrows():
        if l2_per_triangle[idx] > 0:
            # Get triangle coordinates
            trianglex, triangley, direction = twisstntern.utils.return_triangle_coord(
                row['T1_bounds'][0], row['T1_bounds'][1],
                row['T2_bounds'][0], row['T2_bounds'][1],
                row['T3_bounds'][0], row['T3_bounds'][1]
            )
            
            color = cmap(norm(l2_per_triangle[idx]))
            ax.fill(trianglex, triangley, color=color, edgecolor='none')
    
    ax.set_title('L2 Distance per Triangle', fontsize=12, fontweight='bold')
    
    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax = inset_axes(ax, width="3%", height="80%", loc='center right',
                     bbox_to_anchor=(0.05, 0, 1, 1), bbox_transform=ax.transAxes, borderpad=1)
    cbar = plt.colorbar(sm, cax=cax)
    cbar.set_label('L2 Distance', fontsize=10)


def create_comparison_plot(data1, data2, results_df, statistics, alpha, output_path=None):
    """
    Create the main comparison plot with five subplots.
    
    Args:
        data1, data2: Input datasets
        results_df: Analysis results
        statistics: Statistics dictionary
        alpha: Grid granularity
        output_path: Optional path to save the plot
    """
    fig = plt.figure(figsize=(18, 12))
    
    # Create a more complex subplot layout: 2x3 grid with one subplot spanning two columns
    gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1]) 
    ax3 = fig.add_subplot(gs[0, 2])
    ax4 = fig.add_subplot(gs[1, 0])
    ax5 = fig.add_subplot(gs[1, 1:])  # Histogram spans two columns
    
    # Determine consistent scale for data and model plots
    max_prop = max(results_df['prop_data'].max(), results_df['prop_model'].max())
    
    # Plot 1: Data (reference/truth)
    plot_heatmap_data(ax1, results_df, alpha, 'prop_data', 'Data', vmax=max_prop)
    
    # Plot 2: Model
    plot_heatmap_data(ax2, results_df, alpha, 'prop_model', 'Model', vmax=max_prop)
    
    # Plot 3: Residuals
    plot_heatmap_data(ax3, results_df, alpha, 'prop_residual', 'Residuals')
    
    # Plot 4: L2 distances per triangle
    plot_l2_distances(ax4, results_df, alpha)
    
    # Plot 5: Residual histogram
    plot_residual_histogram(ax5, results_df, statistics)
    
    # Add vibrant overall title and statistics
    fig.suptitle(f'Topology Weights Residual Analysis (α = {alpha})', 
                fontsize=17, fontweight='bold', y=0.96, 
                color=sns.color_palette("bright")[0])
    
    # Add statistics text with seaborn styling
    stats_text = (f'L2 distance: {statistics["L2_distance"]:.6f}\n'
                  f'χ² statistic: {statistics["chi2_statistic"]:.6f}\n'
                  f'p-value: {statistics["p_value"]:.6e}\n'
                  f'Data points: {statistics["total_data_points"]}\n'
                  f'Model points: {statistics["total_model_points"]}\n'
                  f'Grid cells: {statistics["num_triangles"]}')
    
    # Vibrant text box
    stats_color = sns.color_palette("pastel")[5]  # Light yellow from pastel palette
    fig.text(0.02, 0.02, stats_text, fontsize=11, verticalalignment='bottom',
             bbox=dict(boxstyle='round,pad=0.5', facecolor=stats_color, 
                      alpha=0.8, edgecolor='gray', linewidth=1))
    
    plt.tight_layout(pad=1.5)  # More padding for seaborn look
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Plot saved to: {output_path}")
    
    return fig


def main():
    """Main function to run the residual analysis."""
    parser = argparse.ArgumentParser(description='Topology Weights Residual Analysis')
    parser.add_argument('--granularity', type=float, default=0.1,
                       help='Grid granularity (default: 0.1)')
    parser.add_argument('--output', type=str, default='residual_analysis_results',
                       help='Output directory (default: residual_analysis_results)')
    parser.add_argument('--data1', type=str, default='migration_topology_weights.csv',
                       help='First dataset filename (default: migration_topology_weights.csv)')
    parser.add_argument('--data2', type=str, default='NoMigration_topology_weights.csv',
                       help='Second dataset filename (default: NoMigration_topology_weights.csv)')
    
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
    
    print("=== Topology Weights Residual Analysis ===")
    print(f"Data 1 (reference): {data1_path}")
    print(f"Data 2 (model): {data2_path}")
    print(f"Granularity: {args.granularity}")
    print(f"Output directory: {output_dir}")
    print()
    
    # Load data
    print("Loading datasets...")
    data1 = load_topology_data(data1_path)
    data2 = load_topology_data(data2_path)
    
    # Perform analysis
    results_df, statistics = perform_grid_analysis(data1, data2, args.granularity)
    
    # Create visualization
    print("Creating visualizations...")
    output_plot = output_dir / f'residual_analysis_granularity_{args.granularity}.png'
    fig = create_comparison_plot(data1, data2, results_df, statistics, 
                               args.granularity, output_plot)
    
    # Save detailed results
    results_csv = output_dir / f'residual_analysis_results_granularity_{args.granularity}.csv'
    results_df.to_csv(results_csv, index=False)
    print(f"Detailed results saved to: {results_csv}")
    
    # Save summary statistics
    stats_file = output_dir / f'summary_statistics_granularity_{args.granularity}.txt'
    with open(stats_file, 'w') as f:
        f.write("Topology Weights Residual Analysis Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Data 1: {data1_path}\n")
        f.write(f"Data 2: {data2_path}\n")
        f.write(f"Granularity: {args.granularity}\n\n")
        
        for key, value in statistics.items():
            if isinstance(value, float):
                f.write(f"{key}: {value:.6e}\n")
            else:
                f.write(f"{key}: {value}\n")
    
    print(f"Summary statistics saved to: {stats_file}")
    
    # Show plot
    plt.show()
    
    print("\nAnalysis complete!")


if __name__ == "__main__":
    main() 