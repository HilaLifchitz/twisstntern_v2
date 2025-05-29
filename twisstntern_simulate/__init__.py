"""
Twisstntern Simulate - A Python package for demographic simulation and tree sequence generation.

This package provides tools for:
- Loading and validating demographic configurations
- Running msprime simulations
- Generating tree sequences
- Saving simulation results

Dependencies:
- msprime
- numpy
- pandas
- yaml
"""

__version__ = '1.0.0'
__author__ = 'Hila Lifchitz'

from .config import Config
from .simulation import (
    run_simulation,
    simulate_locus,
    simulate_chromosome,
    save_tree_sequences
)

# Define what gets imported with "from twisstntern_simulate import *"
__all__ = [
    'Config',
    'run_simulation',
    'simulate_locus',
    'simulate_chromosome',
    'save_tree_sequences'
]

# Check for required dependencies
required_packages = {
    'msprime': 'msprime',
    'numpy': 'numpy',
    'pandas': 'pandas',
    'yaml': 'pyyaml'
}

missing_packages = []
for package, pip_name in required_packages.items():
    try:
        __import__(package)
    except ImportError:
        missing_packages.append(pip_name)

if missing_packages:
    print("Warning: The following required packages are missing:")
    for package in missing_packages:
        print(f"  - {package}")
    print("\nTo install missing packages, run:")
    print(f"pip install {' '.join(missing_packages)}") 