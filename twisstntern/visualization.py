#!/usr/bin/env python
# coding: utf-8

import os
from ete3 import Tree  # Tree object from ete3 to parse and pretty-print Newick trees

os.environ["QT_QPA_PLATFORM"] = "xcb"  # Force XCB backend instead of Wayland
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Add 3D plotting capability
from matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import Polygon
import seaborn as sns  # For color palettes
from twisstntern.utils import (
    cartizian,
    return_triangle_coord,
    T2,
    T3,
    T1_lim,
    T2_lim,
    T3_lim,
    T3_lim_symm,
    h,
    mid_point_triangle,
    right_triangle_coordinates_list,
)
from twisstntern.analysis import fundamental_asymmetry
from sklearn.neighbors import NearestNeighbors
import matplotlib.colors as mcolors

# ============================================================================
# GLOBAL STYLE SETTINGS - TWEAK THESE FOR VISUAL CUSTOMIZATION
# ============================================================================

# Colors for the isoclines of T1, T2, T3 in the index
T1_color = "#7B1E1E"  # for index plot
T2_color = "#277DA1"  # for index plot
T3_color = "#F4A261"  # for index plot

T1_color_data = "lightgrey"  # for the data plot
T2_color_data = "lightgrey"  # for the data plot
T3_color_data = "lightgrey"  # for the data plot


style = "RdBu_r"  # for D-LR plots (=plot results + plot_fundamental_asymmetry)- the colormap
style_heatmap = "viridis_r"  # for heatmap (=plot_ternary_heatmap_data)- the colormap
truncate = False  # for heatmap -whether to chop off the edges of the gradient map


# ============================================================================
# More extensive visual customization guide:
# ============================================================================
#
# QUICK REFERENCE FOR VISUAL TWEAKING:
#
# === COLORMAPS ===
# style = "RdBu_r"              # Main diverging colormap for D-LR plots (red-blue)
# style_heatmap = "Blues"       # Sequential colormap for heatmaps (white to blue)
#
# Popular colormap options:
# - "RdBu_r", "coolwarm", "seismic" (diverging - good for D-LR values)
# - "viridis", "plasma", "inferno", "magma" (sequential - good for heatmaps)
# - "Blues", "Reds", "Greens", "Purples" (sequential - clean and professional)
#
# === GRID LINE COLORS === in the index plot
# T1_color = "#7B1E1E"          # Color for T1 isoclines and labels (crimson)
# T2_color = "#277DA1"          # Color for T2 isoclines and labels (blue)
# T3_color = "#F4A261"          # Color for T3 isoclines and labels (orange)
#
# T1_color_data = "lightgrey"   # Grid color for T1 in data plots
# T2_color_data = "lightgrey"   # Grid color for T2 in data plots
# T3_color_data = "lightgrey"   # Grid color for T3 in data plots
#
# === COMMON COLOR OPTIONS ===
# - "black", "white", "grey", "lightgrey"
# - "#7B1E1E" (crimson), "#277DA1" (blue), "#F4A261" (orange)
# - "#A2C5F2" (light blue), "#ffb347" (orange), "#22223b" (dark blue)
#
# === FIGURE SIZES (in functions) ===
# - plot(): figsize=(8, 6)
# - plot_fundamental_asymmetry(): figsize=(7, 5)
# - plot_results(): figsize=(6, 8)
# - plot_ternary_heatmap_data(): figsize=(8, 6)
#
# === DATA POINT STYLING (in plot function) ===
# - points_color = "#A2C5F2"    # Main data point color
# - alpha = 0.4                 # Transparency (0=invisible, 1=opaque)
# - s = 25                      # Point size
# - edgecolors = "gray"         # Edge color
# - linewidths = 0.2            # Edge line width
#
# === LABEL STYLING ===
# - label_color = "black"       # Corner label color
# - label_size = 12             # Corner label font size
#
# === MARKER STYLING (for significance) ===
# - marker = "*"                # Star marker for p-values
# - s = 9, 22, 25              # Small, medium, large sizes
# - alpha = 0.4, 0.9, 1.0      # Transparency levels
#
# === EMPTY TRIANGLE STYLING ===
# - facecolor = 'white'         # Background color
# - edgecolor = 'grey'          # Border color
# - hatch = '|||'               # Hatch pattern: '|||' (vertical), '///' (diagonal), '---' (horizontal), '+++' (crosses), 'xxx' (diagonal crosses), '...' (dots)
# - linewidth = 0.5             # Border width
#
# === COLORBAR STYLING ===
# - fontsize = 10               # Title font size
# - pad = 6                     # Title padding
# - width = "3%", height = "60%" # Colorbar size
#
# ============================================================================


def save_figure(fig, filename, dpi=300):
    """
    Utility function to save figures with consistent DPI and layout.

    Parameters:
        fig: matplotlib.figure.Figure — the figure object to save
        filename: str — output filename (including .png or .pdf)
        dpi: int — resolution in dots per inch (default: 300)
    """
    import warnings

    # Suppress tight_layout warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*tight_layout.*")
        try:
            fig.tight_layout()  # Try tight_layout, but don't fail if incompatible
        except:
            pass  # Skip tight_layout if not compatible (e.g., with inset axes)

    fig.savefig(filename, dpi=dpi, bbox_inches="tight")


# Initial plotting of data-points in ternary coordinates
def plot(data, granularity, file_name):
    """
    Plot ternary coordinate grid and data points with coordinate lines.

    VISUAL TWEAKING GUIDE:
    - Change the triangle border color/width: see `ax.plot(triangle_x, triangle_y, ...)`
    - Change grid line colors: T1, T2, T3 use T1_color_data, T2_color_data, T3_color_data
    - Change data point color, size, alpha: see `points_color`, `plt.scatter(...)`
    - Change label font, color, and position: see `label_color`, `label_size`, and `plt.text(...)`
    - To show/hide axis coordinate labels, uncomment the relevant block below.
    - To change figure size, edit `fig = plt.figure(figsize=(8, 6))`
    - To change output file name, edit the `title = ...` line at the end.
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

    ax.plot(triangle_x, triangle_y, color="k", linewidth=1, zorder=3)
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
    # points_color = "#A2C5F2" #"my" blue
    # Position 0.4: #2a788e (more blue-green)
    # Position 0.45: #25848e (slightly more blue)
    # Position 0.5: #21918c (the mid-range green you asked about)
    # Position 0.55: #1e9c89 (slightly more green)
    # Position 0.6: #22a884 (more vibrant green)

    points_color = "#22a884"
    plt.scatter(
        x_data,
        y_data,
        facecolors=points_color,
        edgecolors="gray",
        alpha=0.3,
        linewidths=0.4,
        s=15,
    )
    label_color = "black"  # 25.6    tweek with locations of labels
    label_size = 12
    # Label triangle corners
    plt.text(-0.01, 0.88, r"$\mathbf{T}_1$", size=label_size, color=label_color)  # T1
    plt.text(0.51, -0.005, r"$\mathbf{T}_3$", size=label_size, color=label_color)  # T3
    plt.text(
        -0.535, -0.005, r"$\mathbf{T}_2$", size=label_size, color=label_color
    )  # T2

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

    VISUAL TWEAKING GUIDE:
    - Change triangle border color/width: see `ax.plot(triangle_x, triangle_y, ...)`
    - Change colormap for D_LR values: see `get_professional_colormap(style=style, ...)`
    - Change colorbar style, ticks, and label: see the colorbar section
    - Change alpha (transparency) of filled triangles: see `ax.fill(..., alpha=0.8)`
    - Change marker style, color, and size for significance: see `ax.scatter(...)`
    - Change legend position, text, and style: see legend section
    - Change label font, color, and position: see `ax.text(...)`
    - To change figure size, edit `fig, ax = plt.subplots(figsize=(7, 5))`
    - To change output file name, edit the `title = ...` line at the end.
    """
    fig, ax = plt.subplots(figsize=(7, 5))  # Increased width to accommodate colorbar

    # Draw triangle by connecting vertices A (top), B (left), and C (right)
    triangle_x = [0, -0.5, 0.5, 0]  # A → B → C → A
    triangle_y = [h, 0, 0, h]
    ax.plot(triangle_x, triangle_y, color="k", linewidth=1, zorder=3)
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
        ax.scatter(
            x,
            y,
            color="#fde724",
            marker="*",
            alpha=1.0,
            s=25,
            edgecolors="black",
            linewidths=0.2,
        )  # Yellow with black outline
    elif 1e-5 <= main_p_value < 0.001:
        ax.scatter(
            x,
            y,
            color="#5ec961",
            marker="*",
            alpha=1.0,
            s=32,
            edgecolors="black",
            linewidths=0.2,
        )  # Green with black outline
    elif main_p_value < 1e-5:
        ax.scatter(
            x,
            y,
            color="black",
            marker="*",
            alpha=1,
            s=35,
            edgecolors="black",
            linewidths=0.8,
        )  # Black

    # Add colorbar for D_LR values - shorter and positioned lower to match other plots
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax = inset_axes(
        ax,
        width="3%",
        height="25%",
        loc="upper right",
        bbox_to_anchor=(0.05, 0.05, 1, 1),
        bbox_transform=ax.transAxes,
        borderpad=1,
    )
    cbar = plt.colorbar(sm, cax=cax)
    cbar.ax.set_title("D_LR", fontsize=10, pad=6)
    cbar.set_ticks([-1, -0.5, 0, 0.5, 1])
    cbar.ax.set_yticklabels(["-1", "-0.5", "0", "0.5", "1"])

    # Add star marker legend - moved to left side to avoid overlap with colorbar
    legend_items = [
        (0.8, "black", "$p < 10^{-5}$", 35, 1.0, "black"),  # Black
        (0.77, "#5ec961", "$p < 0.001$", 32, 1.0, "black"),  # Green with black outline
        (0.74, "#fde724", "$p < 0.05$", 28, 1.0, "black"),  # Yellow with black outline
    ]
    for y_pos, color, label, size, alpha, edge_color in legend_items:
        ax.scatter(
            -0.45,
            y_pos,
            color=color,
            marker="*",
            alpha=alpha,
            s=size,
            edgecolors=edge_color,
            linewidths=0.8,
        )
        ax.text(-0.43, y_pos, label, size=10)

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
    ax.plot(x_side_T2, T2(0, x_side_T2), color="k", linewidth=1, zorder=3)
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
    label_size = 10  # Reduced from 12 to make labels smaller
    plt.text(-0.02, 0.88, r"$\mathbf{T}_1$", size=label_size, color=T1_color)  # T1
    plt.text(-0.03, -0.01, r"$\mathbf{T}_2$", size=label_size, color=T2_color)  # T2
    plt.text(0.54, -0.01, r"$\mathbf{T}_3$", size=label_size, color=T3_color)  # T3

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

    VISUAL TWEAKING GUIDE:
    - Change the triangle border color/width: see `ax.plot(...)` and `ax.hlines(...)`
    - Change grid line colors: T1, T2, T3 use T1_color, T2_color, T3_color
    - Change the colormap for D-LR values: see `get_professional_colormap(style=style, ...)`
    - Change the colorbar style, ticks, and label: see the colorbar section
    - Change the alpha (transparency) of filled triangles: see `plt.fill(..., alpha=0.8)`
    - Change the marker style, color, and size for significant p-values: see `ax.scatter(...)`
    - Change legend position, text, and style: see legend section
    - To show/hide axis coordinate labels, uncomment the relevant block below.
    - To change figure size, edit `fig = plt.figure(figsize=(6, 8))`
    - To change output file name, edit the `title = ...` line at the end.
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
    ax.plot(x_side_T2, T2(0, x_side_T2), "k", linewidth=1, zorder=3)
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
            # Available hatch patterns: '|||' (vertical), '///' (diagonal), '---' (horizontal),
            # '+++' (crosses), 'xxx' (diagonal crosses), '...' (dots)
            triangle_coords = list(zip(trianglex, triangley))
            empty_triangle = Polygon(
                triangle_coords,
                closed=True,
                facecolor="white",
                edgecolor="grey",
                hatch="|||",  # 🔄 Change this to experiment with different patterns
                linewidth=0.5,
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
                ax.scatter(
                    x,
                    y,
                    color="#fde724",
                    marker="*",
                    alpha=1.0,
                    s=25,
                    edgecolors="black",
                    linewidths=0.2,
                )  # Yellow with black outline
            elif 0.001 > p >= 1e-5:
                ax.scatter(
                    x,
                    y,
                    color="#5ec961",
                    marker="*",
                    alpha=1.0,
                    s=32,
                    edgecolors="black",
                    linewidths=0.2,
                )  # Green with black outline
            elif p < 1e-5:
                ax.scatter(
                    x,
                    y,
                    color="black",
                    marker="*",
                    alpha=1,
                    s=35,
                    edgecolors="black",
                    linewidths=0.2,
                )  # Black

    # Legend for empty triangles - positioned in middle between triangle and colorbar
    # Create a small striped rectangle for the legend
    legend_rect = Polygon(
        [(0.3, 0.85), (0.35, 0.85), (0.35, 0.9), (0.3, 0.9)],
        closed=True,
        facecolor="white",
        edgecolor="grey",
        hatch="|||",  # Should match the pattern used in empty triangles above
        linewidth=0.5,
    )
    ax.add_patch(legend_rect)
    plt.text(0.36, 0.865, "empty triangle", size=8)

    # Legend for p-values - positioned in middle between triangle and colorbar
    ax.scatter(
        0.3,
        0.8,
        color="black",
        marker="*",
        alpha=1,
        s=35,
        edgecolors="black",
        linewidths=0.2,
    )
    plt.text(0.32, 0.8, "$p < 10^{-5}$", size=10)
    ax.scatter(
        0.3,
        0.77,
        color="#5ec961",
        marker="*",
        alpha=1.0,
        s=32,
        edgecolors="black",
        linewidths=0.2,
    )  # Green with black outline
    plt.text(0.32, 0.77, "$p < 0.001$", size=10)
    ax.scatter(
        0.3,
        0.74,
        color="#fde724",
        marker="*",
        alpha=1.0,
        s=30,
        edgecolors="black",
        linewidths=0.2,
    )  # Yellow with black outline
    plt.text(0.32, 0.74, "$p < 0.05$", size=10)

    # Add colorbar for D-LR values - shorter and positioned lower to match other plots
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax = inset_axes(
        ax,
        width="3%",
        height="25%",
        loc="upper right",
        bbox_to_anchor=(0.05, 0.05, 1, 1),
        bbox_transform=ax.transAxes,
        borderpad=1,
    )
    cbar = plt.colorbar(sm, cax=cax)
    cbar.ax.set_title("D_LR", fontsize=10, pad=6)
    cbar.set_ticks([-1, -0.5, 0, 0.5, 1])
    cbar.ax.set_yticklabels(["-1", "-0.5", "0", "0.5", "1"])

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
    # plt.title(style)

    # Save figure
    title = f"{file_name}_analysis_granularity_{alpha}.png"
    save_figure(fig, title)
    return fig


def plot_ternary_heatmap_data(
    data, granularity, file_name, grid_color="#3E3E3E", heatmap_colormap="viridis_r"
):
    """
    Plot a ternary heatmap: each subtriangle is colored by the number of data points it contains.

    VISUAL TWEAKING GUIDE:
    - Change the colormap for the heatmap: see `get_professional_colormap(style=style_heatmap, ...)`
    - Change the color and alpha of grid lines: see `grid_color` and the grid drawing section
    - Change the color, width, and hatch pattern for empty triangles: see the Polygon creation for empty triangles (currently vertical lines)
    - Change the colorbar style, ticks, and label: see the colorbar section
    - To show/hide axis coordinate labels, uncomment the relevant block below.
    - To change figure size, edit `fig = plt.figure(figsize=(8, 6))`
    - To change output file name, edit the `title = ...` line at the end.
    """
    import matplotlib as mpl

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
                    triangles.append(
                        {"T1": (a1, b1), "T2": (a2, b2), "T3": (a3_1, b3_1)}
                    )
                a3_2 = round(a3_1 - alpha, 10)
                b3_2 = round(b3_1 - alpha, 10)
                if a3_2 >= 0:
                    triangles.append(
                        {"T1": (a1, b1), "T2": (a2, b2), "T3": (a3_2, b3_2)}
                    )
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
        n = len(
            data[
                (condition_a1 & (data.T1 <= b1))
                & (condition_a2 & (data.T2 <= b2))
                & (condition_a3 & (data.T3 <= b3))
            ]
        )
        return n

    triangles = create_triangular_grid(alpha)
    total_points = len(data)
    counts = []
    for triangle in triangles:
        count = n_twisstcompare(
            triangle["T1"][0],
            triangle["T1"][1],
            triangle["T2"][0],
            triangle["T2"][1],
            triangle["T3"][0],
            triangle["T3"][1],
            data,
        )
        counts.append(count)
    counts = np.array(counts)
    values = counts  # Always use count, not proportion
    vmin = 1  # Start from 1 instead of 0
    vmax = np.max(values) if np.any(values > 0) else 1

    # Validate colormap parameter
    allowed_colormaps = ["viridis", "viridis_r", "plasma", "inferno", "Blues", "Greys"]
    if heatmap_colormap not in allowed_colormaps:
        print(f"Warning: Unknown colormap '{heatmap_colormap}', using 'viridis_r'")
        heatmap_colormap = "viridis_r"

    # Use the specified colormap for the heatmap
    cmap = get_professional_colormap(style=heatmap_colormap, truncate=False)

    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    fig = plt.figure(figsize=(8, 6))
    ax = plt.axes()
    triangle_x = [0, -0.5, 0.5, 0]
    triangle_y = [h, 0, 0, h]
    ax.plot(triangle_x, triangle_y, color="k", linewidth=1, zorder=3)
    ax.set_xticks([])
    ax.set_yticks([])

    # Add subtle median line (central vertical line) - dotted and thinner like in plot function
    ax.vlines(x=0, ymin=0, ymax=h, colors="gray", linestyles=":", linewidth=1, zorder=2)

    # No grid lines - completely removed for clean heatmap appearance

    # Plot filled triangles (only for nonzero bins)
    for idx, triangle in enumerate(triangles):
        (a1, b1), (a2, b2), (a3, b3) = triangle["T1"], triangle["T2"], triangle["T3"]
        trianglex, triangley, _ = return_triangle_coord(a1, b1, a2, b2, a3, b3)

        if values[idx] == 0:
            # Create plain white triangles for empty regions
            triangle_coords = list(zip(trianglex, triangley))
            empty_triangle = Polygon(
                triangle_coords,
                closed=True,
                facecolor="white",
                edgecolor="none",
                linewidth=0,
            )
            ax.add_patch(empty_triangle)
        else:
            # Use colormap for non-zero values (no restrictions on gradient)
            color = cmap(norm(values[idx]))
            ax.fill(trianglex, triangley, color=color, edgecolor="none")

    # === Use EXACT same labeling and cleanup code as plot() function ===
    label_color = "black"
    label_size = 12
    # Label triangle corners
    plt.text(-0.01, 0.88, r"$\mathbf{T}_1$", size=label_size, color=label_color)  # T1
    plt.text(0.51, -0.005, r"$\mathbf{T}_3$", size=label_size, color=label_color)  # T3
    plt.text(
        -0.535, -0.005, r"$\mathbf{T}_2$", size=label_size, color=label_color
    )  # T2

    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Colorbar (starts from 1) - shorter and positioned lower to match other plots
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cax = inset_axes(
        ax,
        width="3%",
        height="25%",
        loc="upper right",
        bbox_to_anchor=(0.05, 0.05, 1, 1),
        bbox_transform=ax.transAxes,
        borderpad=1,
    )
    cbar = plt.colorbar(sm, cax=cax)
    cbar.ax.set_title("Count", fontsize=10, pad=6)
    cbar.set_ticks([vmin, vmax])
    cbar.ax.set_yticklabels([f"{vmin:.2g}", f"{vmax:.2g}"])

    # plt.title(style_heatmap)
    title = f"{file_name}_heatmap.png"
    save_figure(fig, title)
    return fig


def draw_grey_grid_lines(ax, alpha=0.1):
    """
    Draw grey grid lines on ternary plot, copied from simple_density_plot.py
    Uses fixed granularity of 0.1 as requested.
    Grid lines are drawn with zorder=1 to ensure they appear under data points.
    """
    # Draw grid lines using twisstntern functions (all grey, under data points)
    for i in range(1, int(1 / alpha)):
        y = i * alpha
        # T1 lines (horizontal)
        ax.hlines(
            y=y * h,
            xmin=T1_lim(y)[0],
            xmax=T1_lim(y)[1],
            color="grey",
            linewidth=1,
            zorder=1,
        )
        # T2 lines
        x2 = np.linspace(T2_lim(y)[0], T2_lim(y)[1], 100)
        ax.plot(x2, T2(y, x2), color="grey", linewidth=1, zorder=1)
        # T3 lines
        x3 = np.linspace(T3_lim(y)[0], T3_lim(y)[1], 100)
        ax.plot(x3, T3(y, x3), color="grey", linewidth=1, zorder=1)

    # Central vertical line
    ax.vlines(x=0, ymin=0, ymax=h, colors="grey", ls=":", zorder=1)


def plot_density_colored_radcount(data, file_name, colormap="viridis_r"):
    """
    Plot ternary coordinate grid with data points colored by local density.
    Copied from simple_density_plot.py with fixed parameters as requested.

    Fixed parameters:
    - granularity: 0.1 (always)
    - grid: True (grey grid lines)
    - point_alpha: 0.8
    - density_method: "neighbors"
    - bandwidth: 0.02
    - colormap: User-specified colormap (default: "viridis_r")
    """

    # Fixed parameters as specified
    alpha = 0.1  # Fixed granularity
    grid = True
    point_alpha = 0.8
    density_method = "neighbors"
    bandwidth = 0.02
    # Use the passed colormap parameter instead of global style_heatmap

    fig = plt.figure(figsize=(8, 6))
    ax = plt.axes()

    # === STEP 1: Draw grid lines FIRST (zorder=1) ===
    if grid:
        # Use grey grid lines with fixed granularity 0.1
        draw_grey_grid_lines(ax, alpha=0.1)

    # === STEP 2: Plot scatter points SECOND (zorder=2) ===
    # Convert data points from ternary to cartesian coordinates
    x_data, y_data = cartizian(data["T1"], data["T2"], data["T3"])

    # Calculate density for each point (copied from simple_density_plot.py)
    if density_method == "neighbors":
        points = np.column_stack([x_data, y_data])
        nn = NearestNeighbors(radius=bandwidth)
        nn.fit(points)
        density = nn.radius_neighbors(points, return_distance=False)
        density = np.array([len(neighbors) for neighbors in density])

    # Create scatter plot colored by density (copied from simple_density_plot.py)
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

    # === STEP 3: Draw triangle outline THIRD (zorder=3) ===
    triangle_x = [0, -0.5, 0.5, 0]
    triangle_y = [h, 0, 0, h]
    ax.plot(triangle_x, triangle_y, color="k", linewidth=1, zorder=3)

    # === STEP 4: Draw median line LAST and FAINTLY on top (zorder=4) ===
    ax.vlines(
        x=0,
        ymin=0,
        ymax=h,
        colors="lightgrey",
        linestyles="-",
        linewidth=0.8,
        alpha=0.6,
        zorder=4,
    )
    # Add subtle median line (central vertical line) - dotted and thinner like in plot function
    ax.vlines(x=0, ymin=0, ymax=h, colors="gray", linestyles=":", linewidth=1, zorder=4)

    ax.set_xticks([])
    ax.set_yticks([])

    # === Labeling (using ax.text for proper coordinate system) ===
    label_color = "black"
    label_size = 12
    ax.text(-0.01, 0.88, r"$\mathbf{T}_1$", size=label_size, color=label_color)
    ax.text(0.51, -0.005, r"$\mathbf{T}_3$", size=label_size, color=label_color)
    ax.text(-0.535, -0.005, r"$\mathbf{T}_2$", size=label_size, color=label_color)

    # Remove plot borders (copied from simple_density_plot.py)
    for spine in ax.spines.values():
        spine.set_color("none")

    # Add colorbar - EXACT same style as other plotting functions
    sm = plt.cm.ScalarMappable(
        cmap=colormap,
        norm=plt.cm.colors.Normalize(vmin=density.min(), vmax=density.max()),
    )
    sm.set_array([])
    cax = inset_axes(
        ax,
        width="3%",
        height="25%",
        loc="upper right",
        bbox_to_anchor=(0.05, 0.05, 1, 1),
        bbox_transform=ax.transAxes,
        borderpad=1,
    )
    cbar = plt.colorbar(sm, cax=cax)
    cbar.ax.set_title("Count", fontsize=10, pad=6)

    # Save the plot
    title = f"{file_name}_radcount.png"
    save_figure(fig, title)
    return fig


def get_professional_colormap(style="RdBu_r", truncate=True):
    """
    Returns a professional-looking colormap, optionally truncated.
    """
    # ... colormap generation logic ...
    cmap = sns.color_palette(style, as_cmap=True)
    if truncate:
        # Truncate the colormap to avoid extreme colors
        # Example: use 10-90% of the original range
        return mcolors.ListedColormap(cmap(np.linspace(0.1, 0.9, 256)))
    return cmap


def plot_ternary_projection_stretched(data, expectation, file_name, new_median=0.5):
    """
    Plots a ternary projection where the T2 axis is non-linearly stretched
    around a given expectation value.

    Points with T2 > expectation are pushed further towards the T2 vertex,
    and points with T2 < expectation are compressed away from it. The
    projection is contained within the standard ternary triangle.

    Parameters:
        data (pd.DataFrame): DataFrame with T1, T2, T3 columns.
        expectation (float): The value of T2 to stretch around.
        file_name (str): The prefix for the output plot file.
        new_median (float, optional): The point on the T2 axis to map the
                                     expectation value to. Defaults to 0.5.
    """
    p1 = data["T1"].values
    p2 = data["T2"].values
    p3 = data["T3"].values

    # Define the piecewise linear stretching function for T2
    p2_stretched = np.zeros_like(p2)

    # For p2 < expectation, map [0, expectation] -> [0, new_median]
    mask_lt = p2 < expectation
    if np.any(mask_lt):
        p2_stretched[mask_lt] = new_median * (p2[mask_lt] / expectation)

    # For p2 >= expectation, map [expectation, 1] -> [new_median, 1]
    mask_ge = p2 >= expectation
    if np.any(mask_ge):
        # Avoid division by zero if expectation is 1
        if expectation == 1:
            p2_stretched[mask_ge] = new_median
        else:
            p2_stretched[mask_ge] = new_median + (
                p2[mask_ge] - expectation
            ) * (1 - new_median) / (1 - expectation)

    # To keep the points on the simplex, we re-scale p1 and p3
    # while preserving their ratio.
    p1_plus_p3 = p1 + p3
    p1_stretched = np.zeros_like(p1)
    p3_stretched = np.zeros_like(p3)

    # Avoid division by zero where p1+p3 is 0 (i.e., at the T2 vertex)
    non_zero_mask = p1_plus_p3 > 0

    scale_factor = (1 - p2_stretched[non_zero_mask]) / p1_plus_p3[non_zero_mask]
    p1_stretched[non_zero_mask] = p1[non_zero_mask] * scale_factor
    p3_stretched[non_zero_mask] = p3[non_zero_mask] * scale_factor

    # At the T2 vertex, p1=p3=0, so stretched versions are also 0.
    # This is handled by np.zeros_like initialization.

    # Now, plot the stretched coordinates
    fig = plt.figure(figsize=(8, 6))
    ax = plt.axes()

    # Draw triangle
    triangle_x = [0, -0.5, 0.5, 0]
    triangle_y = [h, 0, 0, h]
    ax.plot(triangle_x, triangle_y, color="k", linewidth=1, zorder=3)
    ax.set_xticks([])
    ax.set_yticks([])

    # Convert stretched ternary coords to cartesian and plot
    x_data, y_data = cartizian(p1_stretched, p2_stretched, p3_stretched)

    points_color = "#22a884"
    plt.scatter(
        x_data,
        y_data,
        facecolors=points_color,
        edgecolors="gray",
        alpha=0.3,
        linewidths=0.4,
        s=15,
    )

    label_color = "black"
    label_size = 12
    plt.text(-0.01, 0.88, r"$\mathbf{T}_1$", size=label_size, color=label_color)
    plt.text(0.51, -0.005, r"$\mathbf{T}_3$", size=label_size, color=label_color)
    plt.text(-0.535, -0.005, r"$\mathbf{T}_2$", size=label_size, color=label_color)

    ax.spines["right"].set_color("none")
    ax.spines["left"].set_color("none")
    ax.spines["bottom"].set_color("none")
    ax.spines["top"].set_color("none")

    title = f"{file_name}_stretched_projection_exp_{expectation}.png"
    save_figure(fig, title)
    return fig

