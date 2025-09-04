#!/usr/bin/env python
# coding: utf-8

"""
Core visualization utilities shared between twisstntern and twisstntern_simulate packages.

This module contains shared visualization helper functions and utilities that are
used by both packages for plotting ternary data and analysis results.
"""

import numpy as np
import matplotlib.pyplot as plt
from .utils import T1_lim, T2_lim, T3_lim, T1, T2, T3, h

# Default colors for isoclines
T1_color_data = "lightgrey"
T2_color_data = "lightgrey"
T3_color_data = "lightgrey"


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


def save_figure(fig, filename, dpi=300):
    """
    Save a matplotlib figure with consistent settings.
    
    Args:
        fig: Matplotlib figure object
        filename (str): Output filename (with or without extension)
        dpi (int): Resolution in dots per inch
    """
    # Ensure .png extension
    if not filename.endswith('.png'):
        filename += '.png'
    
    # Save with consistent settings
    fig.savefig(filename, dpi=dpi, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.close(fig)


def setup_ternary_axes(ax, title=None):
    """
    Configure matplotlib axes for ternary plotting.
    
    Args:
        ax: Matplotlib axes object
        title (str): Optional plot title
    """
    # Set aspect ratio and limits
    ax.set_aspect('equal')
    ax.set_xlim(-0.6, 0.6)
    ax.set_ylim(-0.1, h + 0.1)
    
    # Remove ticks and labels
    ax.set_xticks([])
    ax.set_yticks([])
    
    # Add title if provided
    if title:
        ax.set_title(title, fontsize=14, pad=20)
    
    # Remove spines
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_triangle_boundary(ax, color='black', linewidth=1.5):
    """
    Draw the boundary of the ternary triangle.
    
    Args:
        ax: Matplotlib axes object
        color (str): Line color
        linewidth (float): Line width
    """
    # Triangle vertices: A(top), B(bottom-left), C(bottom-right)
    vertices = np.array([[0, h], [-0.5, 0], [0.5, 0], [0, h]])
    ax.plot(vertices[:, 0], vertices[:, 1], color=color, linewidth=linewidth)