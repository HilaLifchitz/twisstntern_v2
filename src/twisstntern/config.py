#!/usr/bin/env python
# coding: utf-8

from dataclasses import dataclass
from typing import Optional, List
from omegaconf import MISSING


@dataclass
class AnalysisConfig:
    """Configuration for analysis parameters."""
    granularity: float = 0.1
    superfine_granularity: float = 0.05
    fine_granularity: float = 0.1
    coarse_granularity: float = 0.25
    heatmap_granularity: float = 0.02  # Fixed granularity for heatmap


@dataclass
class VisualizationConfig:
    """Configuration for visualization parameters."""
    # Main colormaps
    style: str = "RdBu_r"  # for D-LR plots
    style_heatmap: str = "viridis_r"  # for heatmap
    
    # Colors for isoclines
    t1_color: str = "#7B1E1E"  # for index plot
    t2_color: str = "#277DA1"  # for index plot  
    t3_color: str = "#F4A261"  # for index plot
    
    t1_color_data: str = "lightgrey"  # for data plot
    t2_color_data: str = "lightgrey"  # for data plot
    t3_color_data: str = "lightgrey"  # for data plot
    
    # Plot styling
    figure_size: List[int] = None  # Will default to [8, 7]
    dpi: int = 300
    point_size: float = 2
    point_alpha: float = 0.3
    line_width: float = 0.4
    label_size: int = 14
    
    def __post_init__(self):
        if self.figure_size is None:
            self.figure_size = [8, 7]


@dataclass
class OutputConfig:
    """Configuration for output settings."""
    output_dir: str = "Results"
    log_file: bool = True
    verbose: bool = False


@dataclass
class ProcessingConfig:
    """Configuration for data processing parameters."""
    downsample_n: Optional[int] = None
    downsample_i: Optional[int] = None
    normalize_data: bool = True
    remove_equal_t2_t3: bool = True
    axis_order: Optional[List[str]] = None
    

@dataclass
class TreeProcessingConfig:
    """Configuration for tree file processing."""
    taxon_names: Optional[List[str]] = None
    outgroup: Optional[str] = None
    topology_mapping: Optional[str] = None


@dataclass
class TwisstnternConfig:
    """Main configuration class combining all sub-configs."""
    # Required parameters
    file: str = MISSING
    
    # Sub-configurations
    analysis: AnalysisConfig = AnalysisConfig()
    visualization: VisualizationConfig = VisualizationConfig()
    output: OutputConfig = OutputConfig()
    processing: ProcessingConfig = ProcessingConfig()
    tree_processing: TreeProcessingConfig = TreeProcessingConfig()
    
    # Note: Hydra configuration is handled in config files
