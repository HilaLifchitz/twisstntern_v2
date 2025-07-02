#!/usr/bin/env python3
"""
Simple KDE Density Plot Script for TWISSTNTERN

This is a simplified version for direct experimentation with density plots.
Just edit the parameters below and run it.

Usage:
    python simple_density_plot.py
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import gaussian_kde
from sklearn.neighbors import NearestNeighbors

# Import necessary functions from twisstntern
from twisstntern.utils import (
    cartizian,
    dump_data,
    T1,
    T2,
    T3,
    T1_lim,
    T2_lim,
    T3_lim,
    h,
)

# Import global style settings from twisstntern
from twisstntern.visualization import (
    T1_color_data,
    T2_color_data, 
    T3_color_data,
    save_figure
)

# ============================================================================
# 🎨 EXPERIMENT WITH THESE PARAMETERS!
# ============================================================================
DATA_FILE = "csv files/gene_flow_sims/ne0.45.csv"  # 🔄 Change this to your data file
GRANULARITY = "fine"                                # 🔄 "superfine", "fine", "coarse", or float (e.g., 0.1)                               # 🔄 Change this to experiment!
DENSITY_METHOD = "neighbors"                    # 🔄 "kde" or "neighbors" 
'''
--- "neighbors": Estimate local point density using a fixed-radius k-nearest-neighbors approach:
 For each point, count how many neighbors lie within a specified radius (`bandwidth`).
 This provides a non-kernel, non-parametric density estimate sensitive to local clustering.

--- "kde": Estimate local point density using a kernel density estimate (KDE) approach:
 For each point, count how many neighbors lie within a specified radius (`bandwidth`).
 This provides a non-kernel, non-parametric density estimate sensitive to local clusteri'''
BANDWIDTH = 0.02       
# Small bandwidth (e.g., 0.01):
# - Very sensitive to local variations
# - Shows fine-grained density patterns  
# - Can be "noisy" with sparse data

# Large bandwidth (e.g., 0.1):  
# - Smoother, more general patterns
# - Averages over larger areas
# - Can miss important local clusters                              # 🔄 Sensitivity (smaller = more sensitive)
OUTPUT_PREFIX = "experiment"                         # 🔄 Output filename prefix

COLORMAP = "viridis_r" 
# Popular colormaps to try:
# "viridis", "plasma", "inferno", "magma", "coolwarm", "RdBu_r", 
# "Blues", "Reds", "Greens", "Oranges", "rocket", "mako", "crest"

GRID = True  # 🔄 True = grey grid lines (granularity 0.1), False = no grid lines at all
POINT_ALPHA = 0.8                                    # 🔄 Point transparency (0.0 = invisible, 1.0 = fully opaque)
# ============================================================================


def draw_grey_grid_lines(ax, alpha=0.1):
    """
    Draw grey grid lines on ternary plot, similar to twisstntern visualization.
    Uses fixed granularity of 0.1 as requested.
    Note: Triangle outline is already drawn in black by the main function.
    Grid lines are drawn with zorder=1 to ensure they appear under data points.
    """
    # Draw grid lines using twisstntern functions (all grey, under data points)
    for i in range(1, int(1 / alpha)):
        y = i * alpha
        # T1 lines (horizontal)
        ax.hlines(y=y * h, xmin=T1_lim(y)[0], 
                 xmax=T1_lim(y)[1], color="grey", linewidth=1, zorder=1)
        # T2 lines
        x2 = np.linspace(T2_lim(y)[0], T2_lim(y)[1], 100)
        ax.plot(x2, T2(y, x2), color="grey", linewidth=1, zorder=1)
        # T3 lines
        x3 = np.linspace(T3_lim(y)[0], T3_lim(y)[1], 100)
        ax.plot(x3, T3(y, x3), color="grey", linewidth=1, zorder=1)
    
    # Central vertical line
    ax.vlines(x=0, ymin=0, ymax=h, colors="grey", ls=':', zorder=1)


def plot_density_colored(data, granularity, file_name, density_method="kde", colormap="viridis", bandwidth=0.05, grid=False, point_alpha=0.8):
    """
    Plot ternary coordinate grid with data points colored by local density.
    """
    
    # Map granularity names to alpha values, or parse user-provided float
    if granularity == "superfine":
        alpha = 0.05
    elif granularity == "fine":
        alpha = 0.1
    elif granularity == "coarse":
        alpha = 0.25
    else:
        alpha = float(granularity)  # if the user provides a float value for granularity

    fig = plt.figure(figsize=(8, 6))
    ax = plt.axes()

    # === Use EXACT same triangle drawing code as twisstntern plot() function ===
    triangle_x = [0, -0.5, 0.5, 0]
    triangle_y = [h, 0, 0, h]
    ax.plot(triangle_x, triangle_y, color="k", linewidth=1, zorder=3)
    ax.set_xticks([])
    ax.set_yticks([])

    # Draw grid lines based on grid parameter
    if grid:
        # Use grey grid lines with fixed granularity 0.1
        draw_grey_grid_lines(ax, alpha=0.1)
    # When grid=False, draw no grid lines at all (just the triangle outline already drawn above)

    # Convert data points from ternary to cartesian coordinates
    x_data, y_data = cartizian(data["T1"], data["T2"], data["T3"])
    
    # Calculate density for each point
    if density_method == "kde":
        points = np.vstack([x_data, y_data])
        kde = gaussian_kde(points, bw_method=bandwidth)
        density = kde(points)
    elif density_method == "neighbors":
        points = np.column_stack([x_data, y_data])
        nn = NearestNeighbors(radius=bandwidth)
        nn.fit(points)
        density = nn.radius_neighbors(points, return_distance=False)
        density = np.array([len(neighbors) for neighbors in density])
    else:
        raise ValueError("density_method must be 'kde' or 'neighbors'")
    
    # Create scatter plot colored by density (drawn on top of grid lines)
    scatter = plt.scatter(
        x_data,
        y_data,
        c=density,
        cmap=colormap,
        alpha=point_alpha,
        s=15,
        edgecolors="none",
        zorder=2,
    )
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8, aspect=20)
    
    # Set appropriate label above the colorbar based on density method
    if density_method == "neighbors":
        cbar.ax.set_title("Count", fontsize=10, pad=10)
    else:  # kde
        cbar.ax.set_title("Density", fontsize=10, pad=10)
    
    # === Labeling ===
    label_color = "black"
    label_size = 12
    plt.text(-0.01, 0.88, r"$\mathbf{T}_1$", size=label_size, color=label_color)
    plt.text(0.51, -0.005, r"$\mathbf{T}_3$", size=label_size, color=label_color)
    plt.text(-0.55, -0.005, r"$\mathbf{T}_2$", size=label_size, color=label_color)

    # Remove plot borders
    for spine in ax.spines.values():
        spine.set_color("none")

    # Save the plot
    title = f"{file_name}_density_{density_method}_granularity_{alpha}.png"
    save_figure(fig, title)
    return fig


def main():
    """Run with parameters defined at top of script."""
    
    data_file = Path(DATA_FILE)
    if not data_file.exists():
        print(f"❌ Error: File '{DATA_FILE}' not found!")
        print("   Please edit the DATA_FILE parameter at the top of this script.")
        return 1
    
    try:
        print(f"📊 Loading data from: {data_file}")
        
        try:
            data = dump_data(str(data_file))
        except:
            data = pd.read_csv(data_file)
            
        print(f"✅ Loaded {len(data)} data points")
        
        # Validate required columns
        required_cols = ["T1", "T2", "T3"]
        missing_cols = [col for col in required_cols if col not in data.columns]
        if missing_cols:
            print(f"❌ Error: Missing required columns: {missing_cols}")
            print(f"   Available columns: {list(data.columns)}")
            return 1
        
        print(f"🎨 Generating density plot...")
        print(f"   Granularity: {GRANULARITY}")
        print(f"   Colormap: {COLORMAP}")
        print(f"   Method: {DENSITY_METHOD}")
        print(f"   Bandwidth: {BANDWIDTH}")
        print(f"   Grid: {'Grey (0.1)' if GRID else 'None'}")
        print(f"   Point Alpha: {POINT_ALPHA}")
        
        # Generate the density plot
        fig = plot_density_colored(
            data=data,
            granularity=GRANULARITY,
            file_name=OUTPUT_PREFIX,
            density_method=DENSITY_METHOD,
            colormap=COLORMAP,
            bandwidth=BANDWIDTH,
            grid=GRID,
            point_alpha=POINT_ALPHA
        )
        
        # Get granularity value for filename
        if GRANULARITY == "superfine":
            alpha = 0.05
        elif GRANULARITY == "fine":
            alpha = 0.1
        elif GRANULARITY == "coarse":
            alpha = 0.25
        else:
            alpha = float(GRANULARITY)
        
        output_file = f"{OUTPUT_PREFIX}_density_{DENSITY_METHOD}_granularity_{alpha}.png"
        print(f"✅ Density plot saved: {output_file}")
        
        return 0
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit(main()) 