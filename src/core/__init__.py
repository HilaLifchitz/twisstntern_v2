"""
Core utilities shared between twisstntern and twisstntern_simulate packages.

This module contains common functionality used by both packages:
- Ternary coordinate system mathematics
- Statistical analysis functions  
- Shared visualization utilities
- External dependencies (twisst)
"""

from .utils import (
    # Ternary coordinate functions
    T1, T2, T3,
    T1_lim, T2_lim, T3_lim, T3_lim_symm,
    ternary_coord, cartizian,
    
    # Triangle mathematics
    return_triangle_coord, ref, mid_point_triangle,
    build_conditions, n,
    
    # Data processing
    dump_data,
    
    # Statistical functions
    D_LR, log_likelihood_ratio_test,
    
    # Triangle utilities
    number_triangles, right_triangle_coordinates_list,
)

from .analysis import fundamental_asymmetry, triangles_analysis
from .visualization import get_professional_colormap, draw_isoclines, save_figure, setup_ternary_axes, draw_triangle_boundary

# Logging utilities
from .logger import setup_logging, get_logger, log_system_info, log_analysis_start, log_analysis_complete, log_error, log_simulation_config, log_topologies

# External dependencies
from . import external

__version__ = "0.1.0"
__all__ = [
    # Utils
    "T1", "T2", "T3",
    "T1_lim", "T2_lim", "T3_lim", "T3_lim_symm", 
    "ternary_coord", "cartizian",
    "return_triangle_coord", "ref", "mid_point_triangle",
    "build_conditions", "n",
    "dump_data",
    "D_LR", "log_likelihood_ratio_test",
    "number_triangles", "right_triangle_coordinates_list",
    
    # Analysis
    "fundamental_asymmetry", "triangles_analysis",
    
    # Visualization
    "get_professional_colormap", "draw_isoclines", "save_figure", 
    "setup_ternary_axes", "draw_triangle_boundary",
    
    # Logging
    "setup_logging", "get_logger", "log_system_info", "log_analysis_start", 
    "log_analysis_complete", "log_error", "log_simulation_config", "log_topologies",
    
    # External
    "external",
]