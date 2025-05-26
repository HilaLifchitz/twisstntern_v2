"""
Twisstntern - A package for ternary coordinate analysis and visualization.

This package provides tools for analyzing and visualizing data in ternary coordinates,
including statistical tests, asymmetry analysis, and various plotting functions.

Main Features:
- Data loading and preprocessing for ternary coordinate data
- Statistical analysis of ternary coordinate distributions
- Visualization of ternary coordinate data with various granularity levels
- Fundamental asymmetry analysis
- Triangle-based analysis with customizable granularity

Dependencies:
- numpy
- pandas
- matplotlib
- scipy
- sympy
"""

__version__ = '1.0.0'
__author__ = 'Twisstntern Team'

# Import core functions
from .core import (
    dump_data,
    T1, T2, T3,
    ternary_coord,
    cartizian,
    T1_lim, T2_lim, T3_lim, T3_lim_symm,
    return_triangle_coord,
    ref,
    mid_point_triangle
)

# Import analysis functions
from .analysis import (
    fundemental_asymmetry,
    n,
    D_LR,
    log_likelihood_ratio_test,
    triangles_analysis
)

# Import visualization functions
from .visualization import (
    plot,
    plot_results,
    plot_fundemental_asymmetry,
    plotting_triangle_index
)

# Define what should be imported with "from twisstntern import *"
__all__ = [
    # Core functions
    'dump_data',
    'T1', 'T2', 'T3',
    'ternary_coord',
    'cartizian',
    'T1_lim', 'T2_lim', 'T3_lim', 'T3_lim_symm',
    'return_triangle_coord',
    'ref',
    'mid_point_triangle',
    
    # Analysis functions
    'fundemental_asymmetry',
    'n',
    'D_LR',
    'log_likelihood_ratio_test',
    'triangles_analysis',
    
    # Visualization functions
    'plot',
    'plot_results',
    'plot_fundemental_asymmetry',
    'plotting_triangle_index'
]

# Check for required dependencies
required_packages = {
    'numpy': 'numpy',
    'pandas': 'pandas',
    'matplotlib': 'matplotlib',
    'scipy': 'scipy',
    'sympy': 'sympy'
}

missing_packages = []
for package, import_name in required_packages.items():
    try:
        __import__(import_name)
    except ImportError:
        missing_packages.append(package)

if missing_packages:
    print("Warning: The following required packages are not installed:")
    for package in missing_packages:
        print(f"  - {package}")
    print("\nPlease install them using pip:")
    print(f"pip install {' '.join(missing_packages)}")
