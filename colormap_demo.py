#!/usr/bin/env python
"""
Colormap Demonstration Script for TWISSTNTERN

This script generates heatmaps with all available colormap options 
to help users visualize and choose their preferred colormap style.

Usage:
    python colormap_demo.py data.csv

Requirements:
    - A CSV file with topology weights (T1, T2, T3 columns)
    - TWISSTNTERN package installed
"""

import sys
import os
from pathlib import Path
import argparse

# Import TWISSTNTERN functionality
try:
    import twisstntern
    from twisstntern.utils import dump_data
    from twisstntern.visualization import plot_ternary_heatmap_data
except ImportError as e:
    print(f"Error: Could not import twisstntern. Make sure it's installed.")
    print(f"Import error: {e}")
    sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Generate heatmaps with all available colormap options",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python colormap_demo.py data.csv
    python colormap_demo.py data.csv --output-dir colormap_comparison
    
Available colormaps:
    viridis, viridis_r, plasma, inferno, Blues, Greys
        """
    )
    
    parser.add_argument(
        "data_file",
        type=str,
        help="Path to CSV file with topology weights (T1, T2, T3 columns)"
    )
    
    parser.add_argument(
        "--output-dir",
        type=str,
        default="colormap_demo",
        help="Output directory for colormap comparison plots (default: colormap_demo)"
    )
    
    parser.add_argument(
        "--granularity",
        type=float,
        default=0.02,
        help="Granularity for heatmap (default: 0.02, same as TWISSTNTERN heatmaps)"
    )
    
    args = parser.parse_args()
    
    # Check if data file exists
    data_file = Path(args.data_file)
    if not data_file.exists():
        print(f"Error: Data file not found: {data_file}")
        sys.exit(1)
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    print(f"Output directory: {output_dir}")
    
    # Load data
    print(f"Loading data from: {data_file}")
    try:
        data = dump_data(str(data_file))
        print(f"Loaded {len(data)} data points")
        print(f"Data shape: {data.shape}")
        print(f"Columns: {list(data.columns)}")
        
        # Check for required columns
        required_cols = ['T1', 'T2', 'T3']
        missing_cols = [col for col in required_cols if col not in data.columns]
        if missing_cols:
            print(f"Error: Missing required columns: {missing_cols}")
            print(f"Available columns: {list(data.columns)}")
            sys.exit(1)
            
    except Exception as e:
        print(f"Error loading data: {e}")
        sys.exit(1)
    
    # Available colormap options (same as in the validation)
    colormaps = ["viridis", "viridis_r", "plasma", "inferno", "Blues", "Greys"]
    
    print(f"\nGenerating heatmaps with {len(colormaps)} different colormaps...")
    print("=" * 60)
    
    # Generate a heatmap for each colormap
    successful_plots = 0
    for i, colormap in enumerate(colormaps, 1):
        try:
            print(f"[{i}/{len(colormaps)}] Generating heatmap with '{colormap}' colormap...")
            
            # Use a descriptive output prefix that includes the colormap name
            output_prefix = str(output_dir / f"demo_{colormap}")
            
            # Generate the heatmap (will create: demo_{colormap}_heatmap.png)
            plot_ternary_heatmap_data(
                data=data,
                granularity=args.granularity,
                file_name=output_prefix,
                heatmap_colormap=colormap
            )
            
            successful_plots += 1
            print(f"    ✓ Saved: {output_prefix}_heatmap.png")
            
        except Exception as e:
            print(f"    ✗ Error with {colormap}: {e}")
    
    print("=" * 60)
    print(f"Colormap demonstration complete!")
    print(f"Successfully generated {successful_plots}/{len(colormaps)} heatmaps")
    print(f"\nResults saved to: {output_dir}")
    print("\nGenerated files:")
    
    # List the generated files
    for colormap in colormaps:
        heatmap_file = output_dir / f"demo_{colormap}_heatmap.png"
        if heatmap_file.exists():
            print(f"  ✓ demo_{colormap}_heatmap.png")
        else:
            print(f"  ✗ demo_{colormap}_heatmap.png (failed)")
    
    print(f"\nNow you can compare the different colormap styles and choose your favorite!")
    print(f"To use a specific colormap in your analysis:")
    print(f'  results = twisstntern.run_analysis("data.csv", heatmap_colormap="plasma")')


if __name__ == "__main__":
    main() 