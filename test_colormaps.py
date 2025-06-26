#!/usr/bin/env python3
"""
Test script to visualize different professional colormap options.
Run this to see all available color schemes and choose your favorite.
"""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import Normalize, LinearSegmentedColormap

def get_professional_colormap(style="RdBu_r", truncate=True):
    """
    Get professional diverging colormaps suitable for scientific publication.
    
    Args:
        style (str): One of the following professional colormap styles:
            - "RdBu_r": Red-Blue diverging (reversed) - classic scientific
            - "coolwarm": Cool-warm diverging - Nature journal style
            - "seismic": Seismic diverging - Geophysical style
            - "seismic_soft": Soft seismic - lighter version
            - "seismic_bold": Bold seismic - more vibrant
            - "bwr": Blue-White-Red - clean and professional
            - "PRGn": Purple-Green diverging - colorblind friendly
            - "PiYG": Pink-Yellow-Green diverging - vibrant but professional
            - "viridis": Viridis - modern and colorblind friendly
            - "plasma": Plasma - vibrant purple to yellow
            - "magma": Magma - dark to light
            - "inferno": Inferno - dark to bright
            - "inferno_midrange": Inferno with mid-range center (not white)
            - "inferno_centered": Inferno centered around orange/yellow
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
        "plasma": "plasma",  # Purple to yellow
        "magma": "magma",  # Dark to light
        "inferno": "inferno",  # Dark to bright
    }
    
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

def create_tweaked_seismic_variants():
    """Create tweaked versions of seismic colormap."""
    
    # Get base seismic colormap
    base_seismic = plt.cm.seismic
    
    # Soft seismic - lighter, more pastel
    soft_colors = base_seismic(np.linspace(0.15, 0.85, 256))  # More truncation
    soft_colors = np.clip(soft_colors * 1.3, 0, 1)  # Make lighter
    seismic_soft = LinearSegmentedColormap.from_list('seismic_soft', soft_colors)
    
    # Bold seismic - more vibrant
    bold_colors = base_seismic(np.linspace(0.05, 0.95, 256))  # Less truncation
    bold_colors = np.clip(bold_colors * 0.8, 0, 1)  # Make more saturated
    seismic_bold = LinearSegmentedColormap.from_list('seismic_bold', bold_colors)
    
    # Enhanced seismic - better contrast
    enhanced_colors = base_seismic(np.linspace(0.1, 0.9, 256))
    # Enhance the middle (white) region
    mid_idx = len(enhanced_colors) // 2
    enhanced_colors[mid_idx-10:mid_idx+10] = [1, 1, 1, 1]  # Pure white center
    seismic_enhanced = LinearSegmentedColormap.from_list('seismic_enhanced', enhanced_colors)
    
    return seismic_soft, seismic_bold, seismic_enhanced

def create_inferno_variants():
    """Create tweaked versions of inferno colormap."""
    
    # Get base inferno colormap
    base_inferno = plt.cm.inferno
    
    # Inferno midrange - full range but centers around orange/yellow
    inferno_midrange = LinearSegmentedColormap.from_list('inferno_midrange', 
                                                        base_inferno(np.linspace(0.0, 1.0, 256)))
    
    # Inferno centered - avoids pure black/white
    inferno_centered = LinearSegmentedColormap.from_list('inferno_centered', 
                                                        base_inferno(np.linspace(0.2, 0.8, 256)))
    
    return inferno_midrange, inferno_centered

def create_blue_to_red_variants():
    """Create custom blue to red colormaps."""
    
    # Blue to red continuous
    colors = []
    for i in range(256):
        if i < 128:
            ratio = i / 127.0
            r = ratio * 0.5
            g = 0.0
            b = 1.0 - ratio * 0.5
        else:
            ratio = (i - 128) / 127.0
            r = 0.5 + ratio * 0.5
            g = 0.0
            b = 0.5 - ratio * 0.5
        colors.append([r, g, b, 1.0])
    blue_to_red = LinearSegmentedColormap.from_list('blue_to_red', colors)
    
    # Blue to red centered (avoiding extremes)
    colors_centered = []
    for i in range(256):
        t = 0.2 + 0.6 * (i / 255.0)
        if t < 0.5:
            ratio = (t - 0.2) / 0.3
            r = ratio * 0.5
            g = 0.0
            b = 1.0 - ratio * 0.3
        else:
            ratio = (t - 0.5) / 0.3
            r = 0.5 + ratio * 0.5
            g = 0.0
            b = 0.7 - ratio * 0.7
        colors_centered.append([r, g, b, 1.0])
    blue_to_red_centered = LinearSegmentedColormap.from_list('blue_to_red_centered', colors_centered)
    
    return blue_to_red, blue_to_red_centered

def plot_colormap_comparison():
    """Plot all available professional colormaps for comparison."""
    
    # Standard styles
    styles = ["RdBu_r", "coolwarm", "seismic", "bwr", "PRGn", "PiYG", "RdBu", "BrBG"]
    
    # Create tweaked seismic variants
    seismic_soft, seismic_bold, seismic_enhanced = create_tweaked_seismic_variants()
    
    # Create inferno variants
    inferno_midrange, inferno_centered = create_inferno_variants()
    
    # Create blue to red variants
    blue_to_red, blue_to_red_centered = create_blue_to_red_variants()
    
    # Additional modern styles
    modern_styles = ["viridis", "plasma", "magma", "inferno"]
    
    # Combine all styles
    all_styles = styles + ["seismic_soft", "seismic_bold", "seismic_enhanced"] + modern_styles + ["inferno_midrange", "inferno_centered", "blue_to_red", "blue_to_red_centered"]
    
    # Create subplot grid
    n_cols = 4
    n_rows = (len(all_styles) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(16, 4*n_rows))
    
    if n_rows == 1:
        axes = axes.reshape(1, -1)
    
    # Create sample data (D-LR values from -1 to 1)
    x = np.linspace(-1, 1, 100)
    y = np.linspace(-1, 1, 100)
    X, Y = np.meshgrid(x, y)
    Z = X + Y  # Sample data that goes from -2 to 2
    
    for i, style in enumerate(all_styles):
        row = i // n_cols
        col = i % n_cols
        ax = axes[row, col]
        
        # Get colormap
        if style == "seismic_soft":
            cmap = seismic_soft
        elif style == "seismic_bold":
            cmap = seismic_bold
        elif style == "seismic_enhanced":
            cmap = seismic_enhanced
        elif style == "inferno_midrange":
            cmap = inferno_midrange
        elif style == "inferno_centered":
            cmap = inferno_centered
        elif style == "blue_to_red":
            cmap = blue_to_red
        elif style == "blue_to_red_centered":
            cmap = blue_to_red_centered
        else:
            cmap = get_professional_colormap(style=style, truncate=True)
        
        norm = Normalize(vmin=-1, vmax=1)
        
        # Create heatmap
        im = ax.imshow(Z, cmap=cmap, norm=norm, aspect='auto')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax, shrink=0.8)
        cbar.set_label('D-LR', fontsize=10)
        
        # Set title
        ax.set_title(f'{style}', fontsize=12, fontweight='bold')
        ax.set_xticks([])
        ax.set_yticks([])
    
    # Hide empty subplots
    for i in range(len(all_styles), n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        axes[row, col].set_visible(False)
    
    plt.tight_layout()
    plt.savefig('colormap_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("Colormap comparison saved as 'colormap_comparison.png'")
    print("\n=== Professional colormap recommendations ===")
    print("\n1. SCIENTIFIC PUBLICATION FAVORITES:")
    print("   - 'coolwarm' - Very popular in Nature, Science journals")
    print("   - 'seismic' - Geophysical style, very clean")
    print("   - 'seismic_enhanced' - Seismic with better white center")
    print("   - 'RdBu_r' - Classic scientific red-blue")
    
    print("\n2. MODERN & COLORBLIND FRIENDLY:")
    print("   - 'viridis' - Modern, colorblind friendly, Nature approved")
    print("   - 'PRGn' - Purple-Green, colorblind friendly")
    print("   - 'plasma' - Vibrant purple to yellow")
    
    print("\n3. CLEAN & PROFESSIONAL:")
    print("   - 'bwr' - Simple and clean blue-white-red")
    print("   - 'seismic_soft' - Lighter, more pastel seismic")
    print("   - 'seismic_bold' - More vibrant seismic")
    
    print("\n4. ALTERNATIVE OPTIONS:")
    print("   - 'PiYG' - Pink-Yellow-Green, vibrant but professional")
    print("   - 'BrBG' - Brown-Blue-Green, earth tones")
    print("   - 'magma' - Dark to light gradient")
    print("   - 'inferno' - Dark to bright gradient")
    print("   - 'inferno_midrange' - Inferno with orange/yellow center (not white)")
    print("   - 'inferno_centered' - Inferno avoiding pure black/white")
    print("   - 'blue_to_red' - Custom continuous blue to red")
    print("   - 'blue_to_red_centered' - Blue to red avoiding extremes")

if __name__ == "__main__":
    plot_colormap_comparison() 