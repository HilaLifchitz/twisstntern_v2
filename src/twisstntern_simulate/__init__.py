"""
Twisstntern Simulate - A Python package for demographic simulation and tree sequence generation.

This package provides tools for:
- Loading and validating demographic configurations
- Running msprime simulations
- Generating tree sequences
- Processing tree sequences with twisst
- Running complete twisstntern analysis pipelines

Dependencies:
- msprime
- numpy
- pandas
- yaml
- tskit
- ete3
"""

__version__ = "1.0.0"
__author__ = "Hila Lifchitz"

# Configure logging once for the entire package
import logging
logging.basicConfig(level=logging.INFO)

# Suppress msprime INFO messages
logging.getLogger("msprime").setLevel(logging.WARNING)
logging.getLogger("msprime.ancestry").setLevel(logging.WARNING)

from .hydra_config import TwisstnternSimulateConfig
# Lazy import of pipeline to avoid early dependency loading
# from .pipeline import run_pipeline
from .simulation import run_simulation, simulate_locus, simulate_chromosome

# Define what gets imported with "from twisstntern_simulate import *"
__all__ = [
    "TwisstnternSimulateConfig",
    "run_simulation",
    "simulate_locus", 
    "simulate_chromosome",
]

# Check for required dependencies
required_packages = {
    "msprime": "msprime",
    "numpy": "numpy",
    "pandas": "pandas",
    "yaml": "pyyaml",
    "tskit": "tskit",
    "ete3": "ete3",
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

