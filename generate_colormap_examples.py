#!/usr/bin/env python3
"""
Script to generate example plots for all available colormaps using the real data and real plotting functions.
Each plot is saved to 'colorOptions/' with the colormap name in the filename.
"""
import os
import sys
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import importlib

from twisstntern.visualization import plot_results, plot_fundamental_asymmetry, get_professional_colormap
from twisstntern.pipeline import run_analysis

# Path to your real data
DATA_PATH = "/home/hlifchit/projects/twissting_baby/dataCompare/NoMigration_topology_weights.csv"
OUTPUT_DIR = Path("colorOptions")
OUTPUT_DIR.mkdir(exist_ok=True)

# List of colormaps to test
COLORMAPS = [
    "RdBu_r", "coolwarm", "seismic", "bwr", "PRGn", "PiYG", 
    "RdBu", "BrBG", "viridis", "plasma", "magma", "inferno",
    "seismic_soft", "seismic_bold", "seismic_enhanced",
    "inferno_midrange", "inferno_centered", 
    "blue_to_red", "blue_to_red_centered"
]

def main():
    # Run the analysis pipeline to get results and data
    print("Running analysis pipeline on real data...")
    # We want the results and data, but not to save all the default plots
    # So we will run the pipeline, but patch the plotting functions to no-ops during the run
    import twisstntern.visualization as viz
    import types
    
    # Patch plotting functions to no-ops for the pipeline run
    orig_plot_results = viz.plot_results
    orig_plot_fundamental_asymmetry = viz.plot_fundamental_asymmetry
    viz.plot_results = lambda *a, **k: None
    viz.plot_fundamental_asymmetry = lambda *a, **k: None
    
    # Run pipeline
    results, fundamental_results, _ = run_analysis(
        DATA_PATH,
        granularity=0.1,
        output_dir="colorOptions"
    )
    # Load the data as DataFrame
    data = pd.read_csv(DATA_PATH)
    
    # Restore plotting functions
    viz.plot_results = orig_plot_results
    viz.plot_fundamental_asymmetry = orig_plot_fundamental_asymmetry
    
    # Now, for each colormap, patch get_professional_colormap and generate plots
    for cmap_name in COLORMAPS:
        print(f"Generating plots for colormap: {cmap_name}")
        # Patch get_professional_colormap to always return this colormap
        orig_get_cmap = viz.get_professional_colormap
        def patched_get_cmap(style=None, truncate=True):
            return orig_get_cmap(style=cmap_name, truncate=truncate)
        viz.get_professional_colormap = patched_get_cmap
        
        # Plot results
        try:
            plot_results(results, 0.1, str(OUTPUT_DIR / f"results_{cmap_name}"))
        except Exception as e:
            print(f"Error plotting results for {cmap_name}: {e}")
        try:
            plot_fundamental_asymmetry(data, str(OUTPUT_DIR / f"fundamental_{cmap_name}"))
        except Exception as e:
            print(f"Error plotting fundamental asymmetry for {cmap_name}: {e}")
        
        # Restore get_professional_colormap
        viz.get_professional_colormap = orig_get_cmap
    print("All plots generated in colorOptions/")

if __name__ == "__main__":
    main() 