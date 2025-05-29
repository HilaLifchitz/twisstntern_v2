#!/usr/bin/env python
# coding: utf-8

"""
twisstntern - A package for analyzing ternary data
"""

from .core import (cartizian, return_triangle_coord, dump_data)
from .analysis import (run_analysis, fundemental_asymmetry, 
                      triangles_analysis, number_triangles)
from .visualization import (plot, plot_results, plotting_triangle_index,
                          plot_fundemental_asymmetry)

__version__ = "0.1.0" 