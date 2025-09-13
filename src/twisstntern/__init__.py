#!/usr/bin/env python
# coding: utf-8

"""
twisstntern - A package for analyzing ternary data from topology weights

This package supports:
1. Direct analysis of CSV files containing topology weights
2. Processing tree files (TreeSequence, Newick, Nexus) to generate topology weights
3. Comprehensive ternary analysis and visualization of the results

Tree file formats supported:
- TreeSequence (.trees, .ts)
- Newick (.newick, .nwk, .tree)
- Nexus (.nexus)

The package automatically detects file format and processes accordingly.
"""

from .pipeline import run_analysis, detect_file_type, process_tree_file
# Import from core instead of local modules
from ..core import (
    cartizian,
    return_triangle_coord,
    dump_data,
    number_triangles,
    right_triangle_coordinates_list,
    fundamental_asymmetry,
    triangles_analysis,
)
from .visualization import (
    plot,
    plot_results,
    plotting_triangle_index,
    plot_fundamental_asymmetry,
)
from .tree_processing import (
    detect_and_read_trees,
    simplify_topologies,
    ts_chromosome_to_twisst_weights,
    ts_to_twisst_weights,
    newick_to_twisst_weights,
    trees_to_twisst_weights_unified,
    parse_topology_mapping,
    normalize_topology_string,
    compare_topologies,
    reorder_weights_by_topology_preference,
    print_topology_mapping_with_trees,
)

# Versioning follows Semantic Versioning (semver.org):
# MAJOR.MINOR.PATCH — breaking | feature | fix
__version__ = "0.1.0"