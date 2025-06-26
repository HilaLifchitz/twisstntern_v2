#!/usr/bin/env python
# coding: utf-8

import os
from ete3 import Tree  # Tree object from ete3 to parse and pretty-print Newick trees

os.environ["QT_QPA_PLATFORM"] = "xcb"  # Force XCB backend instead of Wayland
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
from math import sqrt
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Polygon
import seaborn as sns  # For color palettes
from twisstntern.utils import (
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
from twisstntern.analysis import fundamental_asymmetry, triangles_analysis


# Colors for the isoclines of T1, T2, T3 in the index
T1_color = "#7B1E1E"
T2_color = "#277DA1"
T3_color = "#F4A261"

T1_color_data = "lightgrey"
T2_color_data = "lightgrey"
T3_color_data = "lightgrey"


style = "RdBu_r"
style_heatmap = "Blues"

def save_figure(fig, filename, dpi=300):
    """
    Utility function to save figures with consistent DPI and layout.

    Parameters:
        fig: matplotlib.figure.Figure — the figure object to save
        filename: str — output filename (including .png or .pdf)
        dpi: int — resolution in dots per inch (default: 300)
    """
    # Check if figure has 3D subplots
    has_3d = any(hasattr(ax, "zaxis") for ax in fig.get_axes())

    if not has_3d:
        fig.tight_layout()  # Only use tight_layout for 2D plots

    fig.savefig(filename, dpi=dpi, bbox_inches="tight")


# Initial plotting of data-points in ternary coordinates
def plot(data, granularity, file_name):
    """
    Plot ternary coordinate grid and data points with coordinate lines.

    Args:
        data: Dictionary with keys "T1", "T2", "T3" containing ternary coordinates.
        alpha: Step size for plotting coordinate isoclines.
        file_name: Name prefix for saving the output plot image.

    Returns:
        matplotlib.figure.Figure object
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

    # Draw triangle by connecting vertices A (top), B (left), and C (right)
    triangle_x = [0, -0.5, 0.5, 0]  # A → B → C → A
    triangle_y = [h, 0, 0, h]

    ax.plot(triangle_x, triangle_y, color="k", linewidth=1)

    # Hide X and Y axes tick marks
    ax.set_xticks([])
    ax.set_yticks([])

    steps = np.arange(alpha, 1, alpha)
    for y in steps:
        x1_l, x1_r = T1_lim(y)
        ax.hlines(y=y * h, xmin=x1_l, xmax=x1_r, color=T1_color_data, linewidth=1)

        x2_l, x2_r = T2_lim(y)
        x2 = np.linspace(x2_l, x2_r, 100)
        ax.plot(x2, T2(y, x2), color=T2_color_data, linewidth=1)

        x3_l, x3_r = T3_lim(y)
        x3 = np.linspace(x3_l, x3_r, 100)
        ax.plot(x3, T3(y, x3), color=T3_color_data, linewidth=1)

    ax.vlines(x=0, ymin=0, ymax=h, colors="#3E3E3E", linestyles=":")

    # convert data points from ternary to cartezian coordinates and plot them
    x_data, y_data = cartizian(
        data["T1"], data["T2"], data["T3"]
    )  # coordinates of the data points to be plotted
    # plt.scatter(x_data, y_data, color="lightsteelblue", alpha=0.5, s=9)
    points_color = "#A2C5F2"
    #points_color= "green"
    plt.scatter(
        x_data,
        y_data,
        facecolors=points_color,
        edgecolors="gray",
        alpha=0.4,
        linewidths=0.2,
        s=25,
    )
    label_color = "black" # 25.6    tweek with locations of labels
    label_size = 12
    # Label triangle corners
    plt.text(-0.01, 0.88, r"$\mathbf{T}_1$", size=label_size, color=label_color)
    plt.text(0.51, -0.005, r"$\mathbf{T}_3$", size=label_size, color=label_color)
    plt.text(-0.55, -0.005, r"$\mathbf{T}_2$", size=label_size, color=label_color)

    # # printing the coordinates
    # coord = np.arange(0, 1 + alpha, alpha)
    # T_1 = np.arange(0, 0.5 + alpha / 2, alpha / 2)
    # T_2 = np.arange(-0.5, 0 + alpha / 2, alpha / 2)
    # T_3 = np.arange(-0.5, 0.5 + alpha, alpha)

    # for i in range(len(coord)):
    #     c = coord[i]
    #     label = str(round(c, 2))
    #     label_complement = str(round(1 - c, 2))

    #     # T1 label (crimson, offset right)
    #     x1 = T_1[i]
    #     y1 = T2(0, x1)
    #     plt.text(x1 + 0.01, y1, label_complement, size=7, color=T1_color)

    #     # T2 label (dodgerblue, offset left)
    #     x2 = T_2[i]
    #     y2 = T3(0, x2)
    #     plt.text(x2 - 0.04, y2, label_complement, size=7, color=T2_color)

    #     # T3 label (darkgoldenrod, on x-axis)
    #     x3 = T_3[i]
    #     plt.text(x3, -0.03, label, size=7, color=T3_color)

    # removing the box lines around the plot
    ax.spines["right"].set_color("none")
    ax.spines["left"].set_color("none")
    ax.spines["bottom"].set_color("none")
    ax.spines["top"].set_color("none")

    # saving the plot
    title = file_name + "_granuality_" + str(alpha) + ".png"
    save_figure(fig, title)
    return fig


# plotting the asymmetry between the twi main left and right subtriangles
def plot_fundamental_asymmetry(data, file_name):
    """
    Visualizes fundamental asymmetry between left and right subtriangles in ternary space.

    This plot colors the left and right triangle based on the D_LR statistic, and annotates
    the statistical significance of the asymmetry using p-value thresholds (with star markers).

    Args:
        data (pd.DataFrame): DataFrame with columns T1, T2, T3 (ternary coordinates).
        file_name (str): Output filename prefix for the saved PNG figure.

    Returns:
        tuple: (D_LR value, G-test statistic, p-value)
    """
    fig, ax = plt.subplots(figsize=(7, 5))  # Increased width to accommodate colorbar

    # Draw triangle by connecting vertices A (top), B (left), and C (right)
    triangle_x = [0, -0.5, 0.5, 0]  # A → B → C → A
    triangle_y = [h, 0, 0, h]
    ax.plot(triangle_x, triangle_y, color="k", linewidth=1)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.vlines(0, 0, h, colors="black")

    # Define triangle polygons (left and right)
    trianglex_R = [0, 0, 0.5, 0]
    triangley_R = [0, h, 0, 0]
    trianglex_L = [-0.5, 0, 0, -0.5]
    triangley_L = [0, h, 0, 0]

    # Compute statistics
    main_n_r, main_n_l, main_d_lr, main_g_test, main_p_value = fundamental_asymmetry(
        data
    )
    
    # Set up professional colormap for D_LR values
    cmap = get_professional_colormap(style=style, truncate=True)  # Custom blue to red
    norm = Normalize(vmin=-1, vmax=1)  # D_LR ranges from -1 to 1
    
    # Calculate separate D_LR values for left and right triangles
    # Right triangle: use the main D_LR value (positive = more points on right)
    d_lr_right = main_d_lr
    # Left triangle: use the opposite perspective (negative = more points on left)
    d_lr_left = -main_d_lr  # This shows the asymmetry from the left perspective
    
    # Color fill triangles by D_LR score using the colormap
    color_R = cmap(norm(d_lr_right))
    color_L = cmap(norm(d_lr_left))
    
    ax.fill(trianglex_R, triangley_R, color=color_R, alpha=0.8)
    ax.fill(trianglex_L, triangley_L, color=color_L, alpha=0.8)

    # Annotate p-value significance stars
    x, y = 0.15, 0.4 * h
    if 0.001 <= main_p_value < 0.05:
        ax.scatter(x, y, color="darkgoldenrod", marker="*", alpha=0.4, s=9)
    elif 1e-5 <= main_p_value < 0.001:
        ax.scatter(x, y, color="darkslateblue", marker="*", alpha=0.9, s=22)
    elif main_p_value < 1e-5:
        ax.scatter(x, y, color="black", marker="*", alpha=1, s=25)

    # Add colorbar for D_LR values
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax = inset_axes(ax, width="3%", height="60%", loc='center right',
                     bbox_to_anchor=(0.05, 0, 1, 1), bbox_transform=ax.transAxes, borderpad=1)
    cbar = plt.colorbar(sm, cax=cax)
    cbar.ax.set_title('D_LR', fontsize=10, pad=6)
    cbar.set_ticks([-1, -0.5, 0, 0.5, 1])
    cbar.ax.set_yticklabels(['-1', '-0.5', '0', '0.5', '1'])

    # Add star marker legend - moved to left side to avoid overlap with colorbar
    legend_items = [
        (0.8, "black", "$p < 10^{-5}$", 25, 1.0),
        (0.77, "darkslateblue", "$p < 0.001$", 22, 0.9),
        (0.74, "darkgoldenrod", "$p < 0.05$", 9, 0.4),
    ]
    for y_pos, color, label, size, alpha in legend_items:
        ax.scatter(-0.45, y_pos, color=color, marker="*", alpha=alpha, s=size)
        ax.text(-0.43, y_pos, label, size=8)

    # Label sample sizes
    ax.text(-0.5, -0.1, "n =", size=12)
    ax.text(-0.3, -0.1, str(main_n_l), size=12, color="grey")
    ax.text(0.2, -0.1, str(main_n_r), size=12, color="grey")

    # Hide plot borders
    for spine in ax.spines.values():
        spine.set_color("none")

    # Save figure
    title = f"{file_name}_fundamental_asymmetry.png"
    save_figure(fig, title)
    return main_d_lr, main_g_test, main_p_value


# Generating an index plot with distinct numbers assigned to each subtriangle for clear visualization.
# These numbers will correspond to the indices in the result file.
def plotting_triangle_index(granularity, file_name):
    """
    Generate an index plot of right-hand-side subtriangles in the ternary diagram.
    Each triangle is labeled by its index for visual reference.

    Parameters:
        granularity (str or float): One of {"superfine", "fine", "coarse"} or a float value indicating alpha.
    Returns:
        matplotlib Figure object
    """
    # Map preset granularity labels to alpha values and adjust figure size/font accordingly
    if granularity == "superfine":
        alpha = 0.05
        fig = plt.figure(figsize=(7, 6))
        font_size = 7
    elif granularity == "fine":
        alpha = 0.1
        fig = plt.figure(figsize=(5, 4))
        font_size = 8
    elif granularity == "coarse":
        alpha = 0.25
        fig = plt.figure(figsize=(4, 3))
        font_size = 9
    else:
        alpha = float(granularity)
        fig = plt.figure(figsize=(7, 6))
        font_size = 7

    ax = plt.axes()

    # Draw triangle boundary (T2 side + base)
    x_side_T2 = np.linspace(0, 0.5, 100)
    ax.plot(x_side_T2, T2(0, x_side_T2), color="k", linewidth=1)
    ax.hlines(y=0, xmin=0, xmax=0.5, color="k", linewidth=1)

    # Dotted vertical midline (y-axis in ternary diagram)
    ax.vlines(x=0, ymin=0, ymax=h, color="gray", linestyle=":", linewidth=1.2)

    # Hide axis ticks
    ax.set_xticks([])
    ax.set_yticks([])

    # === Plot isoclines ===
    # T2 isoclines (curves from base to left edge)
    for i in range(1, int(1 / (2 * alpha))):
        y = i * alpha
        x2 = np.linspace(0, T2_lim(y)[1], 100)
        ax.plot(x2, T2(y, x2), color=T2_color, linewidth=1)

    # T1 and T3 isoclines (horizontal and right-edge curves)
    for i in range(1, int(1 / alpha)):
        y = i * alpha
        # T1 isocline: horizontal line
        ax.hlines(y=y * h, xmin=0, xmax=T1_lim(y)[1], color=T1_color, linewidth=1)
        # T3 isocline: curved line
        x3 = np.linspace(*T3_lim_symm(y), 100)
        ax.plot(x3, T3(y, x3), color=T3_color, linewidth=1)

    # === Triangle index labels ===
    admissible_triangles = right_triangle_coordinates_list(alpha)
    for _, row in admissible_triangles.iterrows():
        (a1, b1), (a2, b2), (a3, b3) = row[0], row[1], row[2]
        x, y = mid_point_triangle(a1, b1, a2, b2, a3, b3)
        index = row["index"]
        plt.text(x - 0.01, y, str(index), size=font_size)

    # === Triangle corner labels ===
    plt.text(-0.02, 0.88, "T1", size=12, color=T1_color)
    plt.text(-0.03, -0.01, "T2", size=12, color=T2_color)
    plt.text(0.54, -0.01, "T3", size=12, color=T3_color)

    # === Axis coordinate labels ===
    coord = np.arange(0, 1 + alpha, alpha)
    T_1 = np.arange(0, 0.5 + alpha / 2, alpha / 2)
    T_2_3 = np.arange(0, 0.5 + alpha, alpha)

    # T1 labels
    for i in range(len(T_1)):
        label = str(round(1 - coord[i], 2))
        x = T_1[i]
        y = T2(0, x)
        plt.text(x + 0.01, y, label, size=7, color=T1_color)

    # T2 labels
    for i in range(len(T_2_3) - 2):
        label = str(round(coord[i] + alpha, 2))
        y = T2(coord[i] + alpha, 0)
        plt.text(-0.033, y, label, size=7, color=T2_color)

    # T3 labels
    for i in range(len(T_2_3)):
        label = str(round(coord[i] + 0.5, 2))
        x = T_2_3[i]
        plt.text(x, -0.03, label, size=7, color=T3_color)

    # === Hide plot frame ===
    for spine in ["right", "left", "bottom", "top"]:
        ax.spines[spine].set_color("none")

    # Save and return
    title = f"{file_name}_index_granularity_{alpha}.png"
    save_figure(fig, title)
    return fig


# Plot results from triangles_analysis by color-coding subtriangles by D-LR values
# and marking significant p-values.
def plot_results(res, granularity, file_name):
    """
    Plot results from triangles_analysis by color-coding subtriangles by D-LR values
    and marking significant p-values.

    Parameters:
        res (DataFrame): Output of triangles_analysis(data, granularity)
        granularity (str or float): One of {"coarse", "fine", "superfine"} or a float
        file_name (str): Base name to save the plot
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

    # Create plot canvas with increased width for colorbar
    fig = plt.figure(figsize=(6, 8))
    ax = plt.axes()

    # Draw triangle boundaries (T2 and base)
    x_side_T2 = np.linspace(0, 0.5, 100)
    ax.plot(x_side_T2, T2(0, x_side_T2), "k", linewidth=1)
    ax.hlines(y=0, xmin=0, xmax=0.5, color="k", linewidth=1)

    # Remove axis ticks
    ax.set_xticks([])
    ax.set_yticks([])

    # Draw isoclines
    for i in range(1, int(1 / (2 * alpha))):  # T2 isoclines
        y = i * alpha
        x_vals = np.linspace(0, T2_lim(y)[1], 100)
        ax.plot(x_vals, T2(y, x_vals), color=T2_color, linewidth=1)

    for i in range(1, int(1 / alpha)):  # T1 and T3 isoclines
        y = i * alpha
        ax.hlines(y * h, xmin=0, xmax=T1_lim(y)[1], color=T1_color, linewidth=1)  # T1
        x_vals = np.linspace(*T3_lim_symm(y), 100)  # T3
        ax.plot(x_vals, T3(y, x_vals), color=T3_color, linewidth=1)

    # Vertical center line
    ax.vlines(x=0, ymin=0, ymax=h, colors="gray", linestyle=":", linewidth=1.2)

    # Set up professional colormap for D-LR values !!!
    cmap = get_professional_colormap(style=style, truncate=True)  # Custom blue to red
    norm = Normalize(vmin=-1, vmax=1)  # D-LR ranges from -1 to 1

    # Plot subtriangles and highlight significance
    for i, row in res.iterrows():
        (a1, b1), (a2, b2), (a3, b3) = row["coord. (T1, T2, T3)"]

        trianglex, triangley, _ = return_triangle_coord(a1, b1, a2, b2, a3, b3)

        if np.isnan(row["D-LR"]):
            # Create striped pattern for empty triangles using Polygon with hatch
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
        else:
            # Use RdBu_r colormap for D-LR values
            d_lr_value = row["D-LR"]
            color = cmap(norm(d_lr_value))
            plt.fill(trianglex, triangley, color=color, alpha=0.8)

            # Mark significant midpoints
            x, y = mid_point_triangle(a1, b1, a2, b2, a3, b3)
            p = row["p-value(g-test)"]
            if 0.05 > p >= 0.001:
                ax.scatter(x, y, color="yellow", marker="*", alpha=0.4, s=9)
            elif 0.001 > p >= 1e-5:
                ax.scatter(x, y, color="darkslateblue", marker="*", alpha=0.9, s=22)
            elif p < 1e-5:
                ax.scatter(x, y, color="black", marker="*", alpha=1, s=25)

    # Legend for empty triangles - positioned in middle between triangle and colorbar
    # Create a small striped rectangle for the legend
    legend_rect = Polygon(
        [(0.3, 0.85), (0.35, 0.85), (0.35, 0.9), (0.3, 0.9)],
        closed=True,
        facecolor='white',
        edgecolor='grey',
        hatch='///',
        linewidth=0.5
    )
    ax.add_patch(legend_rect)
    plt.text(0.36, 0.865, "empty triangle", size=8)

    # Legend for p-values - positioned in middle between triangle and colorbar
    ax.scatter(0.3, 0.8, color="black", marker="*", alpha=1, s=25)
    plt.text(0.32, 0.8, "$p < 10^{-5}$", size=10)
    ax.scatter(0.3, 0.77, color="darkslateblue", marker="*", alpha=0.9, s=22)
    plt.text(0.32, 0.77, "$p < 0.001$", size=10)
    ax.scatter(0.3, 0.74, color="darkgoldenrod", marker="*", alpha=0.4, s=9)
    plt.text(0.32, 0.74, "$p < 0.05$", size=10)

    # Add colorbar for D-LR values
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax = inset_axes(ax, width="3%", height="60%", loc='center right',
                     bbox_to_anchor=(0.05, 0, 1, 1), bbox_transform=ax.transAxes, borderpad=1)
    cbar = plt.colorbar(sm, cax=cax)
    cbar.ax.set_title('D_LR', fontsize=10, pad=6)
    cbar.set_ticks([-1, -0.5, 0, 0.5, 1])
    cbar.ax.set_yticklabels(['-1', '-0.5', '0', '0.5', '1'])

    # === Annotate triangle corners ===
    # Removed coordinate labels to avoid overlay on heatmap
    # plt.text(-0.02, 0.88, "T1", size=12, color=T1_color)
    # plt.text(-0.03, -0.01, "T2", size=12, color=T2_color)
    # plt.text(0.54, -0.01, "T3", size=12, color=T3_color)

    # === Add axis labels ===
    # Removed axis coordinate labels to avoid overlay on heatmap
    # coord = np.arange(0, 1 + alpha, alpha)
    # T_1 = np.arange(0, 0.5 + alpha / 2, alpha / 2)
    # T_2_3 = np.arange(0, 0.5 + alpha, alpha)

    # # T1 labels
    # for i, x in enumerate(T_1):
    #     label = str(round(1 - coord[i], 2))
    #     y = T2(0, x)
    #     plt.text(x + 0.01, y, label, size=7, color=T1_color)

    # # T2 labels
    # for i in range(len(T_2_3) - 2):
    #     label = str(round(coord[i] + alpha, 2))
    #     y = T2(coord[i] + alpha, 0)
    #     plt.text(-0.033, y, label, size=7, color=T2_color)

    # # T3 labels
    # for i, x in enumerate(T_2_3):
    #     label = str(round(coord[i] + 0.5, 2))
    #     plt.text(x, -0.03, label, size=7, color=T3_color)

    # Remove plot box frame
    for spine in ["right", "left", "bottom", "top"]:
        ax.spines[spine].set_color("none")
    #plt.title(style)     

    # Save figure
    title = f"{file_name}_analysis_granularity_{alpha}.png"
    save_figure(fig, title)
    return fig


def plot_genome_position_2d(data, file_name, genome_positions=None, colormap="inferno"):
    """
    Plot ternary diagram with points colored by genome position.

    Args:
        data: Dictionary or DataFrame with keys/columns "T1", "T2", "T3" containing ternary coordinates.
        file_name: Name prefix for saving the output plot image.
        genome_positions: Array-like of genome positions. If None, uses row indices as positions.
        colormap: Colormap name for genome position coloring (default: 'viridis').

    Returns:
        matplotlib.figure.Figure object
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    # Draw triangle by connecting vertices A (top), B (left), and C (right)
    triangle_x = [0, -0.5, 0.5, 0]  # A → B → C → A
    triangle_y = [h, 0, 0, h]
    ax.plot(triangle_x, triangle_y, color="k", linewidth=2)

    # Hide X and Y axes tick marks
    ax.set_xticks([])
    ax.set_yticks([])

    # Convert data points from ternary to cartesian coordinates
    if isinstance(data, dict):
        x_data, y_data = cartizian(data["T1"], data["T2"], data["T3"])
        n_points = len(data["T1"])
    else:  # DataFrame
        x_data, y_data = cartizian(data["T1"], data["T2"], data["T3"])
        n_points = len(data)

    # Set up genome positions
    if genome_positions is None:
        genome_positions = np.arange(n_points)

    # Create scatter plot colored by genome position
    scatter = ax.scatter(
        x_data,
        y_data,
        c=genome_positions,
        cmap=colormap,
        alpha=0.7,
        s=15,
        edgecolors="none",
    )

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8, aspect=20)
    cbar.set_label("Genome position", rotation=270, labelpad=20, fontsize=12)

    # Label triangle corners
    plt.text(-0.02, 0.88, "T1", size=14, color=T1_color, weight="bold")
    plt.text(0.54, -0.01, "T3", size=14, color=T3_color, weight="bold")
    plt.text(-0.58, -0.01, "T2", size=14, color=T2_color, weight="bold")

    # Remove the box lines around the plot
    for spine in ax.spines.values():
        spine.set_color("none")

    # Set title
    plt.title("Ternary Coordinates Along Genome", fontsize=14, pad=20)

    # Save the plot
    title = f"{file_name}_genome_position_2d.png"
    save_figure(fig, title)
    return fig


def plot_toblerone_3D(
    data, file_name, genome_positions=None, colormap="viridis", subsample=500
):
    """
    Create 3D ternary prism with genome position as the X-axis (horizontal).
    T1 at top (Y-axis), ternary triangle in Y-Z plane, genome extends along X-axis.

    Args:
        data: Dictionary or DataFrame with keys/columns "T1", "T2", "T3" containing ternary coordinates.
        file_name: Name prefix for saving the output plot image.
        genome_positions: Array-like of genome positions. If None, uses row indices as positions.
        colormap: Colormap name for genome position coloring (default: 'viridis').
        subsample: Number of points to show (default: 500). Use None for all points.

    Returns:
        matplotlib.figure.Figure object
    """
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.colors as mcolors

    fig = plt.figure(figsize=(20, 8))

    # Extract coordinates
    if isinstance(data, dict):
        t1, t2, t3 = data["T1"], data["T2"], data["T3"]
        n_points = len(data["T1"])
    else:  # DataFrame
        t1, t2, t3 = data["T1"], data["T2"], data["T3"]
        n_points = len(data)

    # Convert to numpy arrays
    t1, t2, t3 = np.array(t1), np.array(t2), np.array(t3)

    # Set up genome positions
    if genome_positions is None:
        genome_positions = np.linspace(2.7, 25.0, n_points)

    # Subsample data for performance
    if subsample is not None and n_points > subsample:
        indices = np.linspace(0, n_points - 1, subsample, dtype=int)
        t1_sub, t2_sub, t3_sub = t1[indices], t2[indices], t3[indices]
        genome_pos_sub = np.array(genome_positions)[indices]
        print(
            f"Subsampling {n_points} points to {subsample} for ternary prism visualization"
        )
    else:
        t1_sub, t2_sub, t3_sub = t1, t2, t3
        genome_pos_sub = genome_positions

    # Normalize genome positions to X-axis (horizontal length)
    genome_normalized = (genome_pos_sub - np.min(genome_pos_sub)) / (
        np.max(genome_pos_sub) - np.min(genome_pos_sub)
    )
    genome_x = genome_normalized * 25  # Length along X-axis

    # Color setup
    normalize = mcolors.Normalize(
        vmin=np.min(genome_pos_sub), vmax=np.max(genome_pos_sub)
    )
    try:
        colormap_obj = plt.colormaps[colormap]
    except AttributeError:
        colormap_obj = plt.cm.get_cmap(colormap)

    # Define ternary triangle vertices in Y-Z plane
    # T1 at top (high Y), T2 and T3 at bottom left/right (in Z direction)
    triangle_vertices = np.array(
        [
            [h, 0],  # T1 vertex (top, Y-Z coordinates)
            [0, -0.5],  # T2 vertex (bottom-left in Z)
            [0, 0.5],  # T3 vertex (bottom-right in Z)
        ]
    )

    # Create three different viewing angles
    views = [
        {"elev": 15, "azim": 30, "title": "Perspective View"},
        {"elev": 0, "azim": 0, "title": "Side View"},
        {"elev": 0, "azim": 90, "title": "Front View"},
    ]

    for view_idx, view in enumerate(views):
        ax = fig.add_subplot(1, 3, view_idx + 1, projection="3d")

        # Draw the triangular prism frame
        # Front triangle (at x=0)
        front_triangle = np.column_stack([np.zeros(3), triangle_vertices])
        # Back triangle (at x=max_genome)
        back_triangle = np.column_stack(
            [np.full(3, np.max(genome_x)), triangle_vertices]
        )

        # Draw front triangle edges
        front_edges = [[0, 1], [1, 2], [2, 0]]
        for edge in front_edges:
            points = front_triangle[edge]
            ax.plot(
                points[:, 0], points[:, 1], points[:, 2], "k-", linewidth=2, alpha=0.8
            )

        # Draw back triangle edges
        for edge in front_edges:
            points = back_triangle[edge]
            ax.plot(
                points[:, 0], points[:, 1], points[:, 2], "k-", linewidth=2, alpha=0.8
            )

        # Draw connecting edges between front and back
        for i in range(3):
            ax.plot(
                [front_triangle[i, 0], back_triangle[i, 0]],
                [front_triangle[i, 1], back_triangle[i, 1]],
                [front_triangle[i, 2], back_triangle[i, 2]],
                "k-",
                linewidth=2,
                alpha=0.8,
            )

        # Convert ternary coordinates to Y-Z plane coordinates
        # For each point: X = genome position, Y = T1, Z = (T3-T2)/2
        x_points = genome_x  # Genome position along X-axis
        y_points = t1_sub * h  # T1 scaled to triangle height
        z_points = (t3_sub - t2_sub) / 2  # T2-T3 difference scaled

        # Plot all points in the 3D ternary prism
        scatter = ax.scatter(
            x_points,
            y_points,
            z_points,
            c=genome_pos_sub,
            cmap=colormap,
            s=25,
            alpha=0.8,
            edgecolors="black",
            linewidths=0.3,
        )

        # Set viewing angle
        ax.view_init(elev=view["elev"], azim=view["azim"])

        # Add vertex labels
        ax.text(0, h + 0.1, 0, "T1", fontsize=12, color=T1_color, weight="bold")
        ax.text(
            np.max(genome_x),
            -0.1,
            -0.6,
            "T2",
            fontsize=12,
            color=T2_color,
            weight="bold",
        )
        ax.text(
            np.max(genome_x),
            -0.1,
            0.6,
            "T3",
            fontsize=12,
            color=T3_color,
            weight="bold",
        )

        # Add genome position labels along X-axis
        for i, frac in enumerate([0, 0.25, 0.5, 0.75, 1]):
            pos_idx = int(frac * (len(genome_pos_sub) - 1))
            x_pos = frac * np.max(genome_x)
            ax.text(
                x_pos,
                -0.3,
                0,
                f"{genome_pos_sub[pos_idx]:.1f}kb",
                fontsize=8,
                alpha=0.8,
                ha="center",
            )

        ax.set_title(view["title"], fontsize=12)

        # Set axis limits
        ax.set_xlim(-2, np.max(genome_x) + 2)
        ax.set_ylim(-0.4, h + 0.3)
        ax.set_zlim(-0.8, 0.8)

        # Clean up axes
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_zlabel("")

    # Add colorbar
    cbar = fig.colorbar(scatter, ax=fig.get_axes(), shrink=0.8, aspect=40, pad=0.1)
    cbar.set_label("Genome Position (kb)", rotation=270, labelpad=20, fontsize=12)

    # Add genome position scale at bottom
    ax_bottom = fig.add_axes([0.1, 0.02, 0.8, 0.03])
    gradient = np.linspace(0, 1, 256).reshape(1, -1)
    ax_bottom.imshow(gradient, aspect="auto", cmap=colormap)
    ax_bottom.set_xlim(0, 256)
    ax_bottom.set_xticks([0, 64, 128, 192, 256])
    ax_bottom.set_xticklabels(
        [
            f"{genome_pos_sub[0]:.1f}kb",
            f"{genome_pos_sub[len(genome_pos_sub)//4]:.1f}kb",
            f"{genome_pos_sub[len(genome_pos_sub)//2]:.1f}kb",
            f"{genome_pos_sub[3*len(genome_pos_sub)//4]:.1f}kb",
            f"{genome_pos_sub[-1]:.1f}kb",
        ]
    )
    ax_bottom.set_yticks([])
    ax_bottom.set_xlabel("Genome Position", fontsize=12)

    plt.suptitle(
        "Ternary Coordinates in 3D Prism Along Genome (X-axis)", fontsize=16, y=0.95
    )

    # Save the plot
    title = f"{file_name}_genome_ternary_prism_x.png"
    save_figure(fig, title)
    return fig


def plot_toblerone_single_3D(
    data,
    file_name,
    genome_positions=None,
    colormap="viridis",
    subsample=1000,
    view_angle="perspective",
):
    """
    Create single-view 3D ternary prism with genome as X-axis for better detail.

    Args:
        data: Dictionary or DataFrame with keys/columns "T1", "T2", "T3" containing ternary coordinates.
        file_name: Name prefix for saving the output plot image.
        genome_positions: Array-like of genome positions. If None, uses row indices as positions.
        colormap: Colormap name for genome position coloring (default: 'viridis').
        subsample: Number of points to show (default: 1000). Use None for all points.
        view_angle: One of 'perspective', 'side', 'front' for viewing angle.

    Returns:
        matplotlib.figure.Figure object
    """
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.colors as mcolors

    fig = plt.figure(figsize=(15, 10))
    ax = fig.add_subplot(111, projection="3d")

    # Extract coordinates
    if isinstance(data, dict):
        t1, t2, t3 = data["T1"], data["T2"], data["T3"]
        n_points = len(data["T1"])
    else:  # DataFrame
        t1, t2, t3 = data["T1"], data["T2"], data["T3"]
        n_points = len(data)

    # Convert to numpy arrays
    t1, t2, t3 = np.array(t1), np.array(t2), np.array(t3)

    # Set up genome positions
    if genome_positions is None:
        genome_positions = np.linspace(2.7, 25.0, n_points)

    # Subsample data for performance
    if subsample is not None and n_points > subsample:
        indices = np.linspace(0, n_points - 1, subsample, dtype=int)
        t1_sub, t2_sub, t3_sub = t1[indices], t2[indices], t3[indices]
        genome_pos_sub = np.array(genome_positions)[indices]
        print(
            f"Subsampling {n_points} points to {subsample} for ternary prism visualization"
        )
    else:
        t1_sub, t2_sub, t3_sub = t1, t2, t3
        genome_pos_sub = genome_positions

    # Normalize genome positions to X-axis
    genome_normalized = (genome_pos_sub - np.min(genome_pos_sub)) / (
        np.max(genome_pos_sub) - np.min(genome_pos_sub)
    )
    genome_x = genome_normalized * 30  # Length along X-axis

    # Define ternary triangle vertices in Y-Z plane
    triangle_vertices = np.array(
        [
            [h, 0],  # T1 vertex (top in Y-Z plane)
            [0, -0.5],  # T2 vertex (bottom-left in Z)
            [0, 0.5],  # T3 vertex (bottom-right in Z)
        ]
    )

    # Draw the triangular prism frame
    front_triangle = np.column_stack([np.zeros(3), triangle_vertices])
    back_triangle = np.column_stack([np.full(3, np.max(genome_x)), triangle_vertices])

    # Draw front and back triangle edges
    edges = [[0, 1], [1, 2], [2, 0]]
    for edge in edges:
        # Front triangle
        points = front_triangle[edge]
        ax.plot(points[:, 0], points[:, 1], points[:, 2], "k-", linewidth=2)
        # Back triangle
        points = back_triangle[edge]
        ax.plot(points[:, 0], points[:, 1], points[:, 2], "k-", linewidth=2)

    # Draw connecting edges
    for i in range(3):
        ax.plot(
            [front_triangle[i, 0], back_triangle[i, 0]],
            [front_triangle[i, 1], back_triangle[i, 1]],
            [front_triangle[i, 2], back_triangle[i, 2]],
            "k-",
            linewidth=2,
        )

    # Convert ternary to 3D coordinates with genome as X-axis
    x_points = genome_x  # Genome position along X-axis
    y_points = t1_sub * h  # T1 scaled to triangle height (Y-axis)
    z_points = (t3_sub - t2_sub) / 2  # T2-T3 difference (Z-axis)

    scatter = ax.scatter(
        x_points,
        y_points,
        z_points,
        c=genome_pos_sub,
        cmap=colormap,
        s=30,
        alpha=0.8,
        edgecolors="black",
        linewidths=0.3,
    )

    # Set viewing angle
    if view_angle == "perspective":
        ax.view_init(elev=15, azim=30)
    elif view_angle == "side":
        ax.view_init(elev=0, azim=0)
    elif view_angle == "front":
        ax.view_init(elev=0, azim=90)

    # Add labels
    ax.text(0, h + 0.15, 0, "T1", fontsize=16, color=T1_color, weight="bold")
    ax.text(
        np.max(genome_x), -0.15, -0.7, "T2", fontsize=16, color=T2_color, weight="bold"
    )
    ax.text(
        np.max(genome_x), -0.15, 0.7, "T3", fontsize=16, color=T3_color, weight="bold"
    )

    # Add genome position markers along X-axis
    for i, frac in enumerate([0, 0.25, 0.5, 0.75, 1]):
        pos_idx = int(frac * (len(genome_pos_sub) - 1))
        x_pos = frac * np.max(genome_x)
        ax.text(
            x_pos,
            -0.5,
            0,
            f"{genome_pos_sub[pos_idx]:.1f}kb",
            fontsize=10,
            alpha=0.8,
            ha="center",
        )

    # Set limits and clean up
    ax.set_xlim(-3, np.max(genome_x) + 3)
    ax.set_ylim(-0.6, h + 0.4)
    ax.set_zlim(-1.0, 1.0)

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_zlabel("")

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.8, aspect=20)
    cbar.set_label("Genome Position (kb)", rotation=270, labelpad=20, fontsize=12)

    plt.title(
        f"3D Ternary Prism (Genome as X-axis) - {view_angle.title()} View", fontsize=16
    )

    # Save the plot
    title = f"{file_name}_ternary_prism_x_{view_angle}.png"
    save_figure(fig, title)
    return fig


def plot_ternary_heatmap_data(data, granularity, file_name, grid_color="#3E3E3E"):
    """
    Plot a ternary heatmap: each subtriangle is colored by the number of data points it contains.
    - Uses the global 'style' variable for colormap selection
    - grid_color: color for grid lines (default: dark grey)
    - Empty triangles (count = 0) are shown with striped pattern
    """
    import matplotlib as mpl
    from matplotlib.colors import LinearSegmentedColormap
    if granularity == "superfine":
        alpha = 0.05
    elif granularity == "fine":
        alpha = 0.1
    elif granularity == "coarse":
        alpha = 0.25
    else:
        alpha = float(granularity)

    def create_triangular_grid(alpha):
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

    def n_twisstcompare(a1, b1, a2, b2, a3, b3, data):
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

    triangles = create_triangular_grid(alpha)
    total_points = len(data)
    counts = []
    for triangle in triangles:
        count = n_twisstcompare(
            triangle['T1'][0], triangle['T1'][1],
            triangle['T2'][0], triangle['T2'][1],
            triangle['T3'][0], triangle['T3'][1],
            data
        )
        counts.append(count)
    counts = np.array(counts)
    values = counts  # Always use count, not proportion
    vmin = 0  # Include zero in the range
    vmax = np.max(values) if np.any(values > 0) else 1

    # Use the heatmap-specific style variable for colormap
    cmap = get_professional_colormap(style=style_heatmap, truncate=False)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    fig = plt.figure(figsize=(8, 6))
    ax = plt.axes()
    triangle_x = [0, -0.5, 0.5, 0]
    triangle_y = [h, 0, 0, h]
    ax.plot(triangle_x, triangle_y, color="k", linewidth=1)

    # Grid lines
    if grid_color is None:
        steps = np.arange(alpha, 1, alpha)
        for y in steps:
            x1_l, x1_r = T1_lim(y)
            ax.hlines(y=y * h, xmin=x1_l, xmax=x1_r, color=T1_color, linewidth=1)
            x2_l, x2_r = T2_lim(y)
            x2 = np.linspace(x2_l, x2_r, 100)
            ax.plot(x2, T2(y, x2), color=T2_color, linewidth=1)
            x3_l, x3_r = T3_lim(y)
            x3 = np.linspace(x3_l, x3_r, 100)
            ax.plot(x3, T3(y, x3), color=T3_color, linewidth=1)
        ax.vlines(x=0, ymin=0, ymax=h, colors="#3E3E3E", linestyles=":")
    else:
        steps = np.arange(alpha, 1, alpha)
        for y in steps:
            x1_l, x1_r = T1_lim(y)
            ax.hlines(y=y * h, xmin=x1_l, xmax=x1_r, color=grid_color, linewidth=0.8, alpha=0.6)
            x2_l, x2_r = T2_lim(y)
            x2 = np.linspace(x2_l, x2_r, 100)
            ax.plot(x2, T2(y, x2), color=grid_color, linewidth=0.8, alpha=0.6)
            x3_l, x3_r = T3_lim(y)
            x3 = np.linspace(x3_l, x3_r, 100)
            ax.plot(x3, T3(y, x3), color=grid_color, linewidth=0.8, alpha=0.6)
        ax.vlines(x=0, ymin=0, ymax=h, colors=grid_color, linestyles=":", linewidth=1.0, alpha=0.7)

    # Plot filled triangles (only for nonzero bins)
    for idx, triangle in enumerate(triangles):
        (a1, b1), (a2, b2), (a3, b3) = triangle['T1'], triangle['T2'], triangle['T3']
        trianglex, triangley, _ = return_triangle_coord(a1, b1, a2, b2, a3, b3)
        
        if values[idx] == 0:
            # Create striped pattern for empty triangles using Polygon with hatch
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
        else:
            # Use colormap for non-zero values (no restrictions on gradient)
            color = cmap(norm(values[idx]))
            ax.fill(trianglex, triangley, color=color, edgecolor='none')

    # Corner labels
    # Removed coordinate labels to avoid overlay on heatmap
    # plt.text(-0.02, 0.88, "T1", size=12, color=T1_color)
    # plt.text(-0.03, -0.01, "T2", size=12, color=T2_color)
    # plt.text(0.54, -0.01, "T3", size=12, color=T3_color)

    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Colorbar (now includes zero values)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax = inset_axes(ax, width="3%", height="80%", loc='center right',
                     bbox_to_anchor=(0.05, 0, 1, 1), bbox_transform=ax.transAxes, borderpad=1)
    cbar = plt.colorbar(sm, cax=cax)
    cbar.ax.set_title('Count', fontsize=10, pad=6)
    cbar.set_ticks([vmin, vmax])
    cbar.ax.set_yticklabels([f"{vmin:.2g}", f"{vmax:.2g}"])


    #plt.title(style_heatmap)
    title = file_name + "_heatmap_count.png"
    save_figure(fig, title)
    return fig


def get_professional_colormap(style="RdBu_r", truncate=True):
    """
    Get professional diverging colormaps suitable for scientific publication.
    
    Args:
        style (str): One of the following professional colormap styles:
            - "RdBu_r": Red-Blue diverging (reversed) - classic scientific
            - "coolwarm": Cool-warm diverging - Nature journal style
            - "seismic": Seismic diverging - Geophysical style
            - "seismic_soft": Soft seismic - lighter, more pastel version
            - "seismic_bold": Bold seismic - more vibrant version
            - "seismic_enhanced": Enhanced seismic - better white center
            - "bwr": Blue-White-Red - clean and professional
            - "PRGn": Purple-Green diverging - colorblind friendly
            - "PiYG": Pink-Yellow-Green diverging - vibrant but professional
            - "viridis": Viridis - modern and colorblind friendly
            - "plasma": Plasma - vibrant purple to yellow
            - "magma": Magma - dark to light
            - "inferno": Inferno - dark to bright
            - "inferno_r": Reversed inferno
            - "inferno_r_soft": Reversed inferno with subtle yellow start
            - "inferno_midrange": Custom inferno with mid-range center instead of white
            - "inferno_midrange_r": Reversed inferno_midrange
            - "vlag": Violet-Light Blue diverging
            - "PuOr": Purple-Orange diverging
            - "blue_grey_red": Custom Blue-Grey-Red diverging
            - "Greys": Light gray to black - clean, minimal
            - "Blues": White to dark blue - elegant, geoscience
            - "Reds": White to dark red - elegant, warm
            - "Greens": White to dark green - natural, environmental
            - "Oranges": White to dark orange - warm, energetic
            - "YlOrRd": Yellow to orange to red - vibrant sequential
            - "Purples": Light lavender to deep purple - density plots
            - "rocket_r": Seaborn sequential reversed - light red to dark maroon
            - "mako_r": Seaborn sequential reversed - light teal to deep blue
            - "crest_r": Seaborn cool gradient reversed - light aqua to dark blue/green
            - "BuGn_r": Blue-green reversed - light to dark cyan
            - "copper": White to brownish-black - physical sciences
            - "Greys_r": Black to light gray - reversed
            - "Blues_r": Dark blue to white - reversed
            - "Purples_r": Deep purple to light lavender - reversed
            - "rocket": Seaborn sequential - dark maroon to light red
            - "mako": Seaborn sequential - deep blue to light teal
            - "crest": Seaborn cool gradient - dark blue/green to light aqua
            - "BuGn": Blue-green - dark cyan to light
            - "BuGn_r": Blue-green reversed - light to dark cyan
            - "BuViolet": Blue to violet - elegant gradient
            - "BuMagenta": Blue to magenta - vibrant gradient
            - "BuPi": Blue to pink - starts exactly like BuGn, ends in pink
            - "BuPu": Blue to purple - starts exactly like BuGn, ends in dark dramatic purple
        truncate (bool): Whether to truncate the colormap to avoid extreme colors
    
    Returns:
        matplotlib.colors.LinearSegmentedColormap: Professional colormap
    """
    # Professional colormap options
    colormap_options = {
        "RdBu_r": "RdBu_r",  # Red-Blue diverging (reversed) - classic
        "coolwarm": "coolwarm",  # Cool-warm diverging - Nature style
        "seismic": "seismic",  # Seismic diverging - Geophysical
        "bwr": "bwr",  # Blue-White-Red - clean
        "PRGn": "PRGn",  # Purple-Green - colorblind friendly
        "PiYG": "PiYG",  # Pink-Yellow-Green - vibrant but professional
        "RdBu": "RdBu",  # Red-Blue (not reversed) - alternative
        "BrBG": "BrBG",  # Brown-Blue-Green - earth tones
        "viridis": "viridis",  # Modern, colorblind friendly
        "viridis_r": "viridis_r",  # Reversed viridis
        "plasma": "plasma",  # Purple to yellow
        "magma": "magma",  # Dark to light
        "inferno": "inferno",  # Dark to bright
        "inferno_r": "inferno_r",  # Reversed inferno
        "inferno_midrange": "inferno_midrange",  # Custom inferno with mid-range center
        "inferno_midrange_r": "inferno_midrange_r",  # Reversed inferno_midrange
        "vlag": "vlag",  # Violet-Light Blue diverging
        "PuOr": "PuOr",  # Purple-Orange diverging
        "blue_grey_red": "blue_grey_red",  # Custom Blue-Grey-Red diverging
        "Greys": "Greys",  # Light gray to black - clean, minimal
        "Blues": "Blues",  # White to dark blue - elegant, geoscience
        "Reds": "Reds",  # White to dark red - elegant, warm
        "Greens": "Greens",  # White to dark green - natural, environmental
        "Oranges": "Oranges",  # White to dark orange - warm, energetic
        "YlOrRd": "YlOrRd",  # Yellow to orange to red - vibrant sequential
        "Purples": "Purples",  # Light lavender to deep purple - density plots
        "rocket_r": "rocket_r",  # Seaborn sequential reversed - light red to dark maroon
        "mako_r": "mako_r",  # Seaborn sequential reversed - light teal to deep blue
        "crest_r": "crest_r",  # Seaborn cool gradient reversed - light aqua to dark blue/green
        "BuGn_r": "BuGn_r",  # Blue-green reversed - light to dark cyan
        "copper": "copper",  # White to brownish-black - physical sciences
        "Greys_r": "Greys_r",  # Black to light gray - reversed
        "Blues_r": "Blues_r",  # Dark blue to white - reversed
        "Purples_r": "Purples_r",  # Deep purple to light lavender - reversed
        "rocket": "rocket",  # Seaborn sequential - dark maroon to light red
        "mako": "mako",  # Seaborn sequential - deep blue to light teal
        "crest": "crest",  # Seaborn cool gradient - dark blue/green to light aqua
        "BuGn": "BuGn",  # Blue-green - dark cyan to light
        "BuGn_r": "BuGn_r",  # Blue-green reversed - light to dark cyan
        "BuViolet": "BuViolet",  # Custom light blue to violet colormap (starts light like BuGn)
        "BuMagenta": "BuMagenta",  # Custom light blue to magenta colormap (starts light like BuGn)
        "BuPi": "BuPi",  # Custom light blue to pink colormap (starts light like BuGn, ends in pink)
        "BuPu": "BuPu",  # Custom light blue to purple colormap (starts light like BuGn, ends in dark dramatic purple)
    }
    
    # Handle custom seismic variants
    if style == "seismic_soft":
        base_seismic = plt.cm.seismic
        soft_colors = base_seismic(np.linspace(0.15, 0.85, 256))  # More truncation
        soft_colors = np.clip(soft_colors * 1.3, 0, 1)  # Make lighter
        return plt.cm.colors.LinearSegmentedColormap.from_list('seismic_soft', soft_colors)
    
    elif style == "seismic_bold":
        base_seismic = plt.cm.seismic
        bold_colors = base_seismic(np.linspace(0.05, 0.95, 256))  # Less truncation
        bold_colors = np.clip(bold_colors * 0.8, 0, 1)  # Make more saturated
        return plt.cm.colors.LinearSegmentedColormap.from_list('seismic_bold', bold_colors)
    
    elif style == "seismic_enhanced":
        base_seismic = plt.cm.seismic
        enhanced_colors = base_seismic(np.linspace(0.1, 0.9, 256))
        # Enhance the middle (white) region
        mid_idx = len(enhanced_colors) // 2
        enhanced_colors[mid_idx-10:mid_idx+10] = [1, 1, 1, 1]  # Pure white center
        return plt.cm.colors.LinearSegmentedColormap.from_list('seismic_enhanced', enhanced_colors)
    
    elif style == "inferno_midrange":
        # Custom inferno with mid-range center instead of white
        base_inferno = plt.cm.inferno
        # Use the full range but center around orange/yellow instead of white
        colors = base_inferno(np.linspace(0.0, 1.0, 256))
        # The middle (around 0.6-0.7 in inferno) will be orange/yellow instead of white
        return plt.cm.colors.LinearSegmentedColormap.from_list('inferno_midrange', colors)
    
    elif style == "inferno_midrange_r":
        # Reversed version of inferno_midrange
        base_inferno = plt.cm.inferno
        # Use the full range but center around orange/yellow instead of white, then reverse
        colors = base_inferno(np.linspace(0.0, 1.0, 256))
        # Reverse the colors
        colors = colors[::-1]
        return plt.cm.colors.LinearSegmentedColormap.from_list('inferno_midrange_r', colors)
    
    elif style == "inferno_centered":
        # Alternative: use inferno but center it around the orange/yellow region
        base_inferno = plt.cm.inferno
        # Map D-LR [-1, 1] to inferno [0.2, 0.8] to avoid pure black/white
        colors = base_inferno(np.linspace(0.2, 0.8, 256))
        return plt.cm.colors.LinearSegmentedColormap.from_list('inferno_centered', colors)
    
    elif style == "blue_to_red":
        # Custom continuous blue to red colormap
        # Create a smooth transition from blue to red
        colors = []
        for i in range(256):
            # Blue to red transition
            if i < 128:
                # First half: blue to purple
                ratio = i / 127.0
                r = ratio * 0.5  # 0 to 0.5
                g = 0.0
                b = 1.0 - ratio * 0.5  # 1.0 to 0.5
            else:
                # Second half: purple to red
                ratio = (i - 128) / 127.0
                r = 0.5 + ratio * 0.5  # 0.5 to 1.0
                g = 0.0
                b = 0.5 - ratio * 0.5  # 0.5 to 0.0
            colors.append([r, g, b, 1.0])
        return plt.cm.colors.LinearSegmentedColormap.from_list('blue_to_red', colors)
    
    elif style == "blue_to_red_centered":
        # Blue to red with a centered approach (avoiding pure black/white)
        # Similar to inferno_centered but blue to red
        colors = []
        for i in range(256):
            # Map to avoid extremes
            t = 0.2 + 0.6 * (i / 255.0)  # Map [0,1] to [0.2, 0.8]
            
            if t < 0.5:
                # Blue to purple
                ratio = (t - 0.2) / 0.3  # Map [0.2, 0.5] to [0, 1]
                r = ratio * 0.5
                g = 0.0
                b = 1.0 - ratio * 0.3
            else:
                # Purple to red
                ratio = (t - 0.5) / 0.3  # Map [0.5, 0.8] to [0, 1]
                r = 0.5 + ratio * 0.5
                g = 0.0
                b = 0.7 - ratio * 0.7
            colors.append([r, g, b, 1.0])
        return plt.cm.colors.LinearSegmentedColormap.from_list('blue_to_red_centered', colors)
    
    elif style == "blue_grey_red":
        # Custom blue-grey-red diverging colormap
        colors = []
        for i in range(256):
            # Blue to grey to red transition
            if i < 128:
                # First half: blue to grey
                ratio = i / 127.0
                r = ratio * 0.5  # 0 to 0.5
                g = ratio * 0.5  # 0 to 0.5
                b = 1.0 - ratio * 0.5  # 1.0 to 0.5
            else:
                # Second half: grey to red
                ratio = (i - 128) / 127.0
                r = 0.5 + ratio * 0.5  # 0.5 to 1.0
                g = 0.5 - ratio * 0.5  # 0.5 to 0.0
                b = 0.5 - ratio * 0.5  # 0.5 to 0.0
            colors.append([r, g, b, 1.0])
        return plt.cm.colors.LinearSegmentedColormap.from_list('blue_grey_red', colors)
    
    elif style == "BuViolet":
        # Custom light blue to violet colormap (starts light like BuGn)
        colors = []
        for i in range(256):
            # Light blue to violet transition
            ratio = i / 255.0
            # Start with exact light blue/cyan from BuGn (approximately RGB 0.9, 0.95, 0.95)
            # End with a pretty violet (approximately RGB 0.5, 0.0, 0.5)
            r = 0.9 - ratio * 0.4  # 0.9 to 0.5
            g = 0.95 - ratio * 0.95  # 0.95 to 0.0
            b = 0.95 + ratio * 0.05  # 0.95 to 1.0
            colors.append([r, g, b, 1.0])
        return plt.cm.colors.LinearSegmentedColormap.from_list('BuViolet', colors)
    
    elif style == "BuMagenta":
        # Custom light blue to magenta colormap (starts light like BuGn)
        colors = []
        for i in range(256):
            # Light blue to magenta transition
            ratio = i / 255.0
            # Start with exact light blue/cyan from BuGn (approximately RGB 0.9, 0.95, 0.95)
            # End with dignified magenta (approximately RGB 0.8, 0.0, 0.8) - less saturated
            r = 0.9 - ratio * 0.1  # 0.9 to 0.8
            g = 0.95 - ratio * 0.95  # 0.95 to 0.0
            b = 0.95 - ratio * 0.15  # 0.95 to 0.8
            colors.append([r, g, b, 1.0])
        return plt.cm.colors.LinearSegmentedColormap.from_list('BuMagenta', colors)
    
    elif style == "BuPi":
        # Custom light blue to pink colormap (starts exactly like BuGn, ends in pink)
        colors = []
        for i in range(256):
            # Light blue to pink transition
            ratio = i / 255.0
            # Start with exact light blue/cyan from BuGn (approximately RGB 0.9, 0.95, 0.95)
            # End with a darker, sophisticated pink (approximately RGB 0.8, 0.3, 0.5)
            r = 0.9 - ratio * 0.1  # 0.9 to 0.8
            g = 0.95 - ratio * 0.65  # 0.95 to 0.3
            b = 0.95 - ratio * 0.45  # 0.95 to 0.5
            colors.append([r, g, b, 1.0])
        return plt.cm.colors.LinearSegmentedColormap.from_list('BuPi', colors)
    
    elif style == "BuPu":
        # Custom light blue to dark dramatic purple colormap (starts exactly like BuGn, ends in dark purple)
        colors = []
        for i in range(256):
            # Light blue to dark purple transition
            ratio = i / 255.0
            # Start with exact light blue/cyan from BuGn (approximately RGB 0.9, 0.95, 0.95)
            # End with dark dramatic purple (approximately RGB 0.4, 0.0, 0.6)
            r = 0.9 - ratio * 0.5  # 0.9 to 0.4
            g = 0.95 - ratio * 0.95  # 0.95 to 0.0
            b = 0.95 - ratio * 0.35  # 0.95 to 0.6
            colors.append([r, g, b, 1.0])
        return plt.cm.colors.LinearSegmentedColormap.from_list('BuPu', colors)
    
    # Custom inferno_r_soft: reversed inferno with much brighter, softer start
    if style == "inferno_r_soft":
        base = plt.cm.get_cmap("inferno_r")
        # Start at 0.7 for much brighter, and blend first 20% with soft, light yellow
        colors = base(np.linspace(0.7, 1.0, 256))
        soft_yellow = np.array([1.0, 0.97, 0.8, 1.0])  # RGBA for soft, pastel yellow
        n_blend = int(0.2 * len(colors))
        for i in range(n_blend):
            blend_ratio = 1 - (i / n_blend)
            colors[i, :3] = blend_ratio * soft_yellow[:3] + (1 - blend_ratio) * colors[i, :3]
        return plt.cm.colors.LinearSegmentedColormap.from_list('inferno_r_soft', colors)
    
    if style not in colormap_options:
        print(f"Warning: {style} not found, using RdBu_r")
        style = "RdBu_r"
    
    cmap = sns.color_palette(colormap_options[style], as_cmap=True)
    
    if truncate:
        # Truncate the colormap to get lighter colors by cutting off the extremes
        cmap = plt.cm.colors.LinearSegmentedColormap.from_list(
            f'truncated_{style}', 
            cmap(np.linspace(0.1, 0.9, 256))  # Cut off 10% from each end
        )
    
    return cmap
