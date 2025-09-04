#!/usr/bin/env python
# coding: utf-8

import os
from ete3 import Tree  # Tree object from ete3 to parse and pretty-print Newick trees

os.environ["QT_QPA_PLATFORM"] = "xcb"  # Force XCB backend instead of Wayland
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Add 3D plotting capability
from pathlib import Path
import pandas as pd
from math import sqrt
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Polygon
import seaborn as sns  # For color palettes
from .utils import (
    cartizian,
    return_triangle_coord,
    dump_data,
    T1,
    T2,
    T3,
    T1_lim,
    T2_lim,
    T3_lim,
    T3_lim_symm,
    h,
    mid_point_triangle,
    right_triangle_coordinates_list,
    number_triangles,
)
from .analysis import fundamental_asymmetry, triangles_analysis
from sklearn.neighbors import NearestNeighbors

# ============================================================================
# GLOBAL STYLE SETTINGS - TWEAK THESE FOR VISUAL CUSTOMIZATION
# ============================================================================

# Colors for the isoclines of T1, T2, T3 in the index
T1_color = "#7B1E1E" # for index plot
T2_color = "#277DA1" # for index plot
T3_color = "#F4A261" # for index plot

T1_color_data = "lightgrey" # for the data plot
T2_color_data = "lightgrey" # for the data plot
T3_color_data = "lightgrey" # for the data plot

style = "RdBu_r" # for D-LR plots (=plot results + plot_fundamental_asymmetry)- the colormap
style_heatmap = "viridis_r" # for heatmap (=plot_ternary_heatmap_data)- the colormap

# Configurable colormaps:
# Main colormaps to use
# Professional diverging colormaps (red-blue) for D_LR values (-1 to +1):
# style = "RdBu_r"              # Red-white-blue (inverted)
# style = "RdBu"                # Blue-white-red
# style = "RdYlBu_r"            # Red-yellow-blue (inverted)  
# style = "coolwarm"            # Cool-warm (blue to red)
# style = "seismic"             # Red-white-blue (seismic)
# style = "bwr"                 # Blue-white-red

# Professional sequential colormaps for heatmaps (count data):
# style_heatmap = "viridis_r" # for heatmap (=plot_ternary_heatmap_data)- the colormap
# style_heatmap = "Blues"       # Sequential colormap for heatmaps (white to blue)

# Popular colormap options:
# Diverging: 'RdBu_r', 'coolwarm', 'seismic', 'RdYlBu_r'  
# Sequential: 'viridis', 'plasma', 'inferno', 'Blues', 'Greys'
# For count data: 'viridis_r', 'plasma_r', 'Blues', 'Greys'

# Control the size and style of plotted points
# Point styling for data plots:
# - s = 2                       # Point size (1=tiny, 5=small, 10=medium, 20=large)
# - alpha = 0.3                 # Transparency (0=invisible, 1=opaque)
# - marker = 'o'                # Point shape ('o'=circle, '.'=pixel, 's'=square)
# - edgecolor = 'none'          # Point border ('none', 'black', 'white')

# Isocline (grid line) styling:
# - alpha = 0.4                 # Transparency (0=invisible, 1=opaque)
# - linewidths = 0.2            # Edge line width

# Triangle fill styling for results plots:
# - alpha = 0.8                 # Fill transparency 
# - edgecolor = 'white'         # Triangle border color
# - linewidth = 0.1             # Border thickness

# Colorbar/legend styling:
# - alpha = 0.4, 0.9, 1.0      # Transparency levels
# - fontsize = 8, 10, 12       # Text size
# - pad = 0.02                  # Spacing from plot
# - shrink = 0.8                # Colorbar size relative to plot

# Triangle boundary styling:
# - linewidth = 0.5             # Border width
# - color = 'black'             # Border color

# ============================================================================


# Store image files in png format with:
# - High quality (bbox_inches='tight' removes whitespace)
# - High resolution (dpi=300 for publication quality)
def save_plot(file_name, format="png", dpi=300):
    """
    Save the current plot with professional formatting.
    
    Args:
        file_name: Output filename (without extension)
        format: Image format (default: png)
        dpi: int — resolution in dots per inch (default: 300)
    """
    plt.savefig(f"{file_name}.{format}", 
                dpi=dpi, 
                bbox_inches='tight', 
                facecolor='white', 
                edgecolor='none')


# ======================================================================================
# CORE PLOTTING FUNCTIONS
# ======================================================================================

def plot(data, granularity, file_name):
    """
    Plot ternary data with customizable granularity grid overlay.
    
    Creates a scatter plot of data points on a ternary diagram with optional grid lines
    showing isoclines for constant values of T1, T2, and T3 coordinates.
    
    Args:
        data (pd.DataFrame): Data with T1, T2, T3 columns
        granularity (float or str): Grid spacing ('fine', 'coarse', or numeric value)
        file_name (str): Output filename prefix
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

    # Extract Cartesian coordinates from ternary data
    x, y = cartizian(data["T1"], data["T2"], data["T3"])

    # Create figure
    fig, ax = plt.subplots(figsize=(8, 7))
    
    # Draw triangle boundary
    triangle_x = [0, -0.5, 0.5, 0]  # A → B → C → A
    triangle_y = [h, 0, 0, h]
    ax.plot(triangle_x, triangle_y, color='black', linewidth=1.5)
    
    # Draw optional isoclines (grid lines) if granularity is fine enough
    if alpha <= 0.25:
        draw_isoclines(ax, alpha)
    
    # Plot data points
    ax.scatter(x, y, 
              s=2, 
              alpha=0.3,
              linewidths=0.4,
              color='#21918c',  # Professional teal color
              edgecolors='none')
    
    # Add axis labels
    label_size = 14
    label_color = 'black'
    ax.text(-0.01, 0.88, r"$\mathbf{T}_1$", size=label_size, color=label_color) #T1
    ax.text(0.51, -0.005, r"$\mathbf{T}_3$", size=label_size, color=label_color) #T3
    ax.text(-0.535, -0.005, r"$\mathbf{T}_2$", size=label_size, color=label_color) #T2
    
    # Set equal aspect ratio and remove axes
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Save plot
    save_plot(f"{file_name}_granuality_{alpha}")
    plt.close()


def plot_fundamental_asymmetry(data, file_name):
    """
    Visualize fundamental asymmetry analysis with colored triangular regions.
    
    This function creates a ternary plot showing the left/right division and colors
    each half according to its D_LR asymmetry score, providing an immediate visual
    assessment of directional bias in the data.
    
    Key visual elements:
    - Left triangle: colored according to negative D_LR values
    - Right triangle: colored according to positive D_LR values  
    - Color intensity represents magnitude of asymmetry
    - Statistical significance indicators with star markers
    
    Customization options:
    - Change colormap for D_LR values: see `get_professional_colormap(style=style, ...)`
    - Modify significance star colors/sizes: see star plotting section
    - Change alpha (transparency) of filled triangles: see `ax.fill(..., alpha=0.8)`
    
    Args:
        data (pd.DataFrame): Ternary data with T1, T2, T3 columns
        file_name (str): Output filename prefix
    """
    # Perform fundamental asymmetry analysis
    main_n_r, main_n_l, main_d_lr, main_g_test, main_p_value = fundamental_asymmetry(data)
    
    fig, ax = plt.subplots(figsize=(8, 7))
    
    # Draw triangle boundary
    triangle_x = [0, -0.5, 0.5, 0]  # A → B → C → A
    triangle_y = [h, 0, 0, h]
    ax.plot(triangle_x, triangle_y, color='black', linewidth=1.5)
    
    # Draw vertical dividing line  
    ax.axvline(x=0, color='black', linewidth=1, alpha=0.7)
    
    # Define left and right triangle regions
    trianglex_R = [0, 0, 0.5, 0]
    triangley_R = [0, h, 0, 0]
    trianglex_L = [-0.5, 0, 0, -0.5]
    triangley_L = [0, 0, h, 0]
    
    # Set up professional colormap for D_LR values
    cmap = get_professional_colormap(style=style, truncate=True)  # Custom blue to red
    norm = Normalize(vmin=-1, vmax=1)
    
    # Map D_LR scores to colors (right=positive, left=negative)
    color_R = cmap(norm(abs(main_d_lr)))  # Right gets positive color
    color_L = cmap(norm(-abs(main_d_lr)))  # Left gets negative color
    
    # Color fill triangles by D_LR score using the colormap
    ax.fill(trianglex_R, triangley_R, color=color_R, alpha=0.8)
    ax.fill(trianglex_L, triangley_L, color=color_L, alpha=0.8)
    
    # Add statistical significance indicators with stars
    x, y = 0.15, 0.4 * h
    if 0.001 <= main_p_value < 0.05:
        ax.scatter(x, y, marker="*", s=18, color="#fde724", alpha=1.0)
    elif 1e-5 <= main_p_value < 0.001:
        ax.scatter(x, y, marker="*", s=22, color="#5ec961", alpha=1.0)
    elif main_p_value < 1e-5:
        ax.scatter(x, y, marker="*", s=25, color="black", alpha=1.0)
    
    # Add colorbar
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, 
                        shrink=0.8, pad=0.02)
    cbar.set_label('D_LR', rotation=270, labelpad=15)
    cbar.set_ticks([-1, -0.5, 0, 0.5, 1])
    cbar.ax.set_yticklabels(['-1', '-0.5', '0', '0.5', '1'])
    
    # Add significance legend with star markers
    for y_pos, color, label, size, alpha in [
        (0.8, "black", "$p < 10^{-5}$", 25, 1.0),
        (0.77, "#5ec961", "$p < 0.001$", 22, 1.0),  # Super bright green from viridis
        (0.74, "#fde724", "$p < 0.05$", 18, 1.0),   # Brightest yellow from viridis
    ]:
        ax.scatter(-0.45, y_pos, color=color, marker="*", alpha=alpha, s=size)
        ax.text(-0.43, y_pos, label, size=8)
    
    # Add sample size annotations
    ax.text(-0.5, -0.1, "n =", size=12)
    ax.text(-0.3, -0.1, str(main_n_l), size=12, color="grey")
    ax.text(0.2, -0.1, str(main_n_r), size=12, color="grey")
    
    # Set equal aspect ratio and remove axes  
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Save plot
    save_plot(f"{file_name}_fundamental_asymmetry")
    plt.close()


def plotting_triangle_index(granularity, file_name):
    """
    Create an index plot showing the triangular grid with isoclines and triangle numbering.
    
    This creates a reference diagram that shows:
    - The triangular subdivision at the specified granularity
    - Isoclines for T1, T2, T3 coordinates 
    - Index numbers for each triangle (for cross-referencing with analysis results)
    
    Args:
        granularity (float or str): Grid resolution
        file_name (str): Output filename prefix  
    """
    # Map granularity names to alpha values
    if granularity == "superfine":
        alpha = 0.05
    elif granularity == "fine":
        alpha = 0.1
    elif granularity == "coarse":
        alpha = 0.25
    else:
        alpha = float(granularity)
    
    fig, ax = plt.subplots(figsize=(8, 7))
    
    # Draw triangle boundary
    triangle_x = [0, -0.5, 0.5, 0]
    triangle_y = [h, 0, 0, h]
    ax.plot(triangle_x, triangle_y, color='black', linewidth=1.5)
    
    # Draw isoclines with different colors for each coordinate
    draw_isoclines(ax, alpha, 
                   colors={'T1': T1_color, 'T2': T2_color, 'T3': T3_color})
    
    # Add triangle index numbers
    coordinate_list = right_triangle_coordinates_list(alpha)
    for i, coord_data in coordinate_list.iterrows():
        coords = coord_data[0:3]
        triangle_index = coord_data["index"]
        
        # Get triangle center for label placement
        (a1, b1), (a2, b2), (a3, b3) = coords
        center_x, center_y = mid_point_triangle(a1, b1, a2, b2, a3, b3)
        
        # Add index number
        ax.text(center_x, center_y, str(triangle_index), 
                ha='center', va='center', fontsize=8, 
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
    
    # Add axis labels
    label_size = 14
    ax.text(-0.01, 0.88, r"$\mathbf{T}_1$", size=label_size, color=T1_color)
    ax.text(0.51, -0.005, r"$\mathbf{T}_3$", size=label_size, color=T3_color)
    ax.text(-0.535, -0.005, r"$\mathbf{T}_2$", size=label_size, color=T2_color)
    
    # Set equal aspect ratio and remove axes
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Save plot
    save_plot(f"{file_name}_index_granularity_{alpha}")
    plt.close()


def plot_results(results, granularity, file_name):
    """
    Visualize triangle analysis results with color-coded asymmetry values.
    
    Creates a ternary plot where each subtriangle is filled with a color representing
    its D_LR asymmetry score. This provides a spatial map of asymmetry patterns
    across the ternary space.
    
    Args:
        results (pd.DataFrame): Triangle analysis results with D_LR values
        granularity (float or str): Grid resolution used in analysis
        file_name (str): Output filename prefix
    """
    # Map granularity to alpha value
    if granularity == "superfine":
        alpha = 0.05
    elif granularity == "fine":
        alpha = 0.1
    elif granularity == "coarse":
        alpha = 0.25
    else:
        alpha = float(granularity)
    
    fig, ax = plt.subplots(figsize=(8, 7))
    
    # Set up colormap
    cmap = get_professional_colormap(style=style, truncate=True)
    norm = Normalize(vmin=-1, vmax=1)
    
    # Plot each triangle with its D_LR color
    for i, row in results.iterrows():
        if pd.isna(row['D-LR']):
            continue
            
        coords = row['coord. (T1, T2, T3)']
        d_lr = row['D-LR']
        
        # Get triangle coordinates
        (a1, b1), (a2, b2), (a3, b3) = coords
        trianglex, triangley, _ = return_triangle_coord(a1, b1, a2, b2, a3, b3)
        
        # Color triangle by D_LR value
        color = cmap(norm(d_lr))
        ax.fill(trianglex, triangley, color=color, alpha=0.8, 
                edgecolor='white', linewidth=0.1)
    
    # Draw triangle boundary
    triangle_x = [0, -0.5, 0.5, 0]
    triangle_y = [h, 0, 0, h]
    ax.plot(triangle_x, triangle_y, color='black', linewidth=1.5)
    
    # Add colorbar
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, 
                        shrink=0.8, pad=0.02)
    cbar.set_label('D_LR', rotation=270, labelpad=15)
    
    # Add axis labels
    label_size = 14
    ax.text(-0.01, 0.88, r"$\mathbf{T}_1$", size=label_size)
    ax.text(0.51, -0.005, r"$\mathbf{T}_3$", size=label_size)
    ax.text(-0.535, -0.005, r"$\mathbf{T}_2$", size=label_size)
    
    # Set equal aspect ratio and remove axes
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Save plot
    save_plot(f"{file_name}_analysis_granularity_{alpha}")
    plt.close()


def plot_ternary_heatmap_data(data, granularity, file_name, heatmap_colormap="viridis_r"):
    """
    Create a ternary heatmap showing data density across the triangular space.
    
    Args:
        data (pd.DataFrame): Ternary data with T1, T2, T3 columns
        granularity (float): Grid resolution for binning
        file_name (str): Output filename prefix
        heatmap_colormap (str): Colormap for the heatmap
    """
    # Convert to Cartesian coordinates
    x, y = cartizian(data["T1"], data["T2"], data["T3"])
    
    fig, ax = plt.subplots(figsize=(8, 7))
    
    # Create hexagonal binning for density
    hb = ax.hexbin(x, y, gridsize=int(1/granularity), cmap=heatmap_colormap, alpha=0.8)
    
    # Draw triangle boundary
    triangle_x = [0, -0.5, 0.5, 0]
    triangle_y = [h, 0, 0, h]
    ax.plot(triangle_x, triangle_y, color='black', linewidth=1.5)
    
    # Add colorbar
    cb = fig.colorbar(hb, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label('Count', rotation=270, labelpad=15)
    
    # Add axis labels
    label_size = 14
    ax.text(-0.01, 0.88, r"$\mathbf{T}_1$", size=label_size)
    ax.text(0.51, -0.005, r"$\mathbf{T}_3$", size=label_size)
    ax.text(-0.535, -0.005, r"$\mathbf{T}_2$", size=label_size)
    
    # Set equal aspect ratio and remove axes
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Save plot
    save_plot(f"{file_name}_heatmap_count_granularity_{granularity}")
    plt.close()


def plot_density_colored_radcount(data, file_name):
    """
    Create a density plot with radial count coloring.
    
    Args:
        data (pd.DataFrame): Ternary data with T1, T2, T3 columns
        file_name (str): Output filename prefix
    """
    # Convert to Cartesian coordinates
    x, y = cartizian(data["T1"], data["T2"], data["T3"])
    
    fig, ax = plt.subplots(figsize=(8, 7))
    
    # Create scatter plot with density coloring
    ax.scatter(x, y, c=range(len(x)), cmap='viridis', alpha=0.6, s=2)
    
    # Draw triangle boundary
    triangle_x = [0, -0.5, 0.5, 0]
    triangle_y = [h, 0, 0, h]
    ax.plot(triangle_x, triangle_y, color='black', linewidth=1.5)
    
    # Add axis labels
    label_size = 14
    ax.text(-0.01, 0.88, r"$\mathbf{T}_1$", size=label_size)
    ax.text(0.51, -0.005, r"$\mathbf{T}_3$", size=label_size)
    ax.text(-0.535, -0.005, r"$\mathbf{T}_2$", size=label_size)
    
    # Set equal aspect ratio and remove axes
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Save plot
    save_plot(f"{file_name}_density_radcount")
    plt.close()


# ======================================================================================
# HELPER FUNCTIONS
# ======================================================================================

def get_professional_colormap(style="RdBu_r", truncate=False):
    """
    Get a professional colormap for D_LR visualization.
    
    Args:
        style (str): Colormap name
        truncate (bool): Whether to truncate extreme values
        
    Returns:
        matplotlib colormap object
    """
    cmap = plt.cm.get_cmap(style)
    if truncate:
        # Truncate to avoid pure white in center for better visibility
        from matplotlib.colors import ListedColormap
        colors = cmap(np.linspace(0.1, 0.9, 256))
        cmap = ListedColormap(colors)
    return cmap


def draw_isoclines(ax, alpha, colors=None):
    """
    Draw isoclines (grid lines) on the ternary plot.
    
    Args:
        ax: Matplotlib axes object
        alpha (float): Grid spacing
        colors (dict): Colors for T1, T2, T3 lines (optional)
    """
    if colors is None:
        colors = {'T1': T1_color_data, 'T2': T2_color_data, 'T3': T3_color_data}
    
    # Draw T1 isoclines (horizontal)
    for t1 in np.arange(alpha, 1, alpha):
        x_l, x_r = T1_lim(t1)
        ax.plot([x_l, x_r], [t1*h, t1*h], color=colors['T1'], alpha=0.4, linewidth=0.2)
    
    # Draw T2 isoclines (left-leaning)
    for t2 in np.arange(alpha, 1, alpha):
        x_l, x_r = T2_lim(t2)
        x_vals = np.linspace(x_l, x_r, 100)
        y_vals = T2(t2, x_vals)
        ax.plot(x_vals, y_vals, color=colors['T2'], alpha=0.4, linewidth=0.2)
    
    # Draw T3 isoclines (right-leaning)
    for t3 in np.arange(alpha, 1, alpha):
        x_l, x_r = T3_lim(t3)
        x_vals = np.linspace(x_l, x_r, 100)
        y_vals = T3(t3, x_vals)
        ax.plot(x_vals, y_vals, color=colors['T3'], alpha=0.4, linewidth=0.2)