#!/usr/bin/env python
# coding: utf-8

import os
os.environ['QT_QPA_PLATFORM'] = 'xcb'  # Force XCB backend instead of Wayland
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
from math import sqrt
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


# Colors for the isoclines of T1, T2, T3
T1_color = "#7B1E1E"
T2_color = "#277DA1"
T3_color = "#F4A261"


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
        ax.hlines(y=y * h, xmin=x1_l, xmax=x1_r, color=T1_color, linewidth=1)

        x2_l, x2_r = T2_lim(y)
        x2 = np.linspace(x2_l, x2_r, 100)
        ax.plot(x2, T2(y, x2), color=T2_color, linewidth=1)

        x3_l, x3_r = T3_lim(y)
        x3 = np.linspace(x3_l, x3_r, 100)
        ax.plot(x3, T3(y, x3), color=T3_color, linewidth=1)

    ax.vlines(x=0, ymin=0, ymax=h, colors="#3E3E3E", linestyles=":")

    # convert data points from ternary to cartezian coordinates and plot them
    x_data, y_data = cartizian(
        data["T1"], data["T2"], data["T3"]
    )  # coordinates of the data points to be plotted
    # plt.scatter(x_data, y_data, color="lightsteelblue", alpha=0.5, s=9)
    plt.scatter(
        x_data,
        y_data,
        facecolors="#A2C5F2",
        edgecolors="gray",
        alpha=0.7,
        linewidths=0.5,
        s=10,
    )

    # Label triangle corners
    plt.text(-0.02, 0.88, "T1", size=12, color=T1_color)
    plt.text(0.54, -0.01, "T3", size=12, color=T3_color)
    plt.text(-0.58, -0.01, "T2", size=12, color=T2_color)

    # printing the coordinates
    coord = np.arange(0, 1 + alpha, alpha)
    T_1 = np.arange(0, 0.5 + alpha / 2, alpha / 2)
    T_2 = np.arange(-0.5, 0 + alpha / 2, alpha / 2)
    T_3 = np.arange(-0.5, 0.5 + alpha, alpha)

    for i in range(len(coord)):
        c = coord[i]
        label = str(round(c, 2))
        label_complement = str(round(1 - c, 2))

        # T1 label (crimson, offset right)
        x1 = T_1[i]
        y1 = T2(0, x1)
        plt.text(x1 + 0.01, y1, label_complement, size=7, color=T1_color)

        # T2 label (dodgerblue, offset left)
        x2 = T_2[i]
        y2 = T3(0, x2)
        plt.text(x2 - 0.04, y2, label_complement, size=7, color=T2_color)

        # T3 label (darkgoldenrod, on x-axis)
        x3 = T_3[i]
        plt.text(x3, -0.03, label, size=7, color=T3_color)

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
    fig, ax = plt.subplots(figsize=(6, 4))

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
    d_lr_color_score_R = (main_d_lr + 1) / 2
    d_lr_color_score_L = (
        (main_n_l - 0.5 * (main_n_l + main_n_r)) / (0.5 * (main_n_l + main_n_r)) + 1
    ) / 2

    # Color fill triangles by D_LR score
    ax.fill(
        trianglex_R,
        triangley_R,
        color=(d_lr_color_score_R, 1 - d_lr_color_score_R, 1 - d_lr_color_score_R),
    )
    ax.fill(
        trianglex_L,
        triangley_L,
        color=(d_lr_color_score_L, 1 - d_lr_color_score_L, 1 - d_lr_color_score_L),
    )

    # Annotate p-value significance stars
    x, y = 0.15, 0.4 * h
    if 0.001 <= main_p_value < 0.05:
        ax.scatter(x, y, color="darkgoldenrod", marker="*", alpha=0.4, s=9)
    elif 1e-5 <= main_p_value < 0.001:
        ax.scatter(x, y, color="darkslateblue", marker="*", alpha=0.9, s=22)
    elif main_p_value < 1e-5:
        ax.scatter(x, y, color="black", marker="*", alpha=1, s=25)

    # Add color scale legend
    for i, score in enumerate([0, 0.5, 1]):
        ax.add_patch(
            plt.Rectangle(
                (0.2 + i * 0.05, 0.85), 0.05, 0.05, color=(score, 1 - score, 1 - score)
            )
        )
    ax.text(0.135, 0.865, "$D_{lr} =-1$", size=7)
    ax.text(0.26, 0.865, "0", size=7)
    ax.text(0.31, 0.865, "1", size=7)

    # Add star marker legend
    legend_items = [
        (0.8, "black", "$p < 10^{-5}$", 25, 1.0),
        (0.77, "darkslateblue", "$p < 0.001$", 22, 0.9),
        (0.74, "darkgoldenrod", "$p < 0.05$", 9, 0.4),
    ]
    for y_pos, color, label, size, alpha in legend_items:
        ax.scatter(0.2, y_pos, color=color, marker="*", alpha=alpha, s=size)
        ax.text(0.214, y_pos, label, size=8)

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

    # Create plot canvas
    fig = plt.figure(figsize=(5, 8))
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

    # Plot subtriangles and highlight significance
    for i, row in res.iterrows():
        (a1, b1), (a2, b2), (a3, b3) = row["coord. (T1, T2, T3)"]

        trianglex, triangley, _ = return_triangle_coord(a1, b1, a2, b2, a3, b3)

        if np.isnan(row["D-LR"]):
            plt.fill(trianglex, triangley, color="black")
        else:
            score = (row["D-LR"] + 1) / 2  # normalize to [0,1]
            plt.fill(trianglex, triangley, color=(score, 1 - score, 1 - score))

            # Mark significant midpoints
            x, y = mid_point_triangle(a1, b1, a2, b2, a3, b3)
            p = row["p-value(g-test)"]
            if 0.05 > p >= 0.001:
                ax.scatter(x, y, color="yellow", marker="*", alpha=0.4, s=9)
            elif 0.001 > p >= 1e-5:
                ax.scatter(x, y, color="darkslateblue", marker="*", alpha=0.9, s=22)
            elif p < 1e-5:
                ax.scatter(x, y, color="black", marker="*", alpha=1, s=25)

    # === Add legends ===
    # Legend for D-LR color gradient
    for j, val in enumerate([0, 0.5, 1]):
        ax.add_patch(
            plt.Rectangle(
                (0.2 + 0.05 * j, 0.85), 0.05, 0.05, color=(val, 1 - val, 1 - val)
            )
        )
    ax.add_patch(plt.Rectangle((0.518, 0.85), 0.05, 0.05, color="black"))
    plt.text(0.16, 0.865, "$D_{lr} = -1$", size=8)
    plt.text(0.26, 0.865, "0", size=8)
    plt.text(0.31, 0.865, "1", size=8)
    plt.text(0.38, 0.865, "empty triangle", size=8)

    # Legend for p-values
    ax.scatter(0.2, 0.8, color="black", marker="*", alpha=1, s=25)
    plt.text(0.214, 0.8, "$p < 10^{-5}$", size=8)
    ax.scatter(0.2, 0.77, color="darkslateblue", marker="*", alpha=0.9, s=22)
    plt.text(0.214, 0.77, "$p < 0.001$", size=8)
    ax.scatter(0.2, 0.74, color="darkgoldenrod", marker="*", alpha=0.4, s=9)
    plt.text(0.214, 0.74, "$p < 0.05$", size=8)

    # === Annotate triangle corners ===
    plt.text(-0.02, 0.88, "T1", size=12, color=T1_color)
    plt.text(-0.03, -0.01, "T2", size=12, color=T2_color)
    plt.text(0.54, -0.01, "T3", size=12, color=T3_color)

    # === Add axis labels ===
    coord = np.arange(0, 1 + alpha, alpha)
    T_1 = np.arange(0, 0.5 + alpha / 2, alpha / 2)
    T_2_3 = np.arange(0, 0.5 + alpha, alpha)

    # T1 labels
    for i, x in enumerate(T_1):
        label = str(round(1 - coord[i], 2))
        y = T2(0, x)
        plt.text(x + 0.01, y, label, size=7, color=T1_color)

    # T2 labels
    for i in range(len(T_2_3) - 2):
        label = str(round(coord[i] + alpha, 2))
        y = T2(coord[i] + alpha, 0)
        plt.text(-0.033, y, label, size=7, color=T2_color)

    # T3 labels
    for i, x in enumerate(T_2_3):
        label = str(round(coord[i] + 0.5, 2))
        plt.text(x, -0.03, label, size=7, color=T3_color)

    # Remove plot box frame
    for spine in ["right", "left", "bottom", "top"]:
        ax.spines[spine].set_color("none")

    # Save figure
    title = f"{file_name}_analysis_granularity_{alpha}.png"
    save_figure(fig, title)
    return fig


def plot_genome_position_2d(data, file_name, genome_positions=None, colormap="viridis"):
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
