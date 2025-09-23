"""
Hydra structured configuration dataclasses for twisstntern_simulate.

These dataclasses define the structure and types for all configuration
parameters, providing type safety and validation.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Union, Any


@dataclass 
class PopulationConfig:
    """Configuration for a single population."""
    name: str
    Ne: float
    growth_rate: float = 0.0
    sample_size: Optional[int] = None


@dataclass
class SplitConfig:
    """Configuration for a population split event."""
    time: float
    derived_pop1: str
    derived_pop2: str
    ancestral_pop: str


@dataclass
class DownsamplingConfig:
    """Configuration for downsampling options."""
    tree_interval: Optional[int] = None
    tree_start_index: int = 0
    position_interval_kb: Optional[int] = None
    position_start_kb: int = 0


@dataclass
class FundamentalAsymmetryConfig:
    """Configuration for fundamental asymmetry analysis."""
    filter_t2_equals_t3: bool = True


@dataclass
class AnalysisConfig:
    """Configuration for analysis parameters."""
    granularity: Union[float, str] = 0.1
    downsampling: DownsamplingConfig = field(default_factory=DownsamplingConfig)
    topology_mapping: Optional[str] = None
    fundamental_asymmetry: FundamentalAsymmetryConfig = field(default_factory=FundamentalAsymmetryConfig)


@dataclass
class PlotSettingsConfig:
    """Configuration for general plot settings."""
    figure_size: List[float] = field(default_factory=lambda: [10, 8])
    dpi: int = 300
    alpha: float = 0.7


@dataclass
class ColorsConfig:
    """Configuration for color schemes."""
    T1_color: str = "#7B1E1E"
    T2_color: str = "#277DA1"
    T3_color: str = "#F4A261"
    T1_color_data: str = "lightgrey"
    T2_color_data: str = "lightgrey"
    T3_color_data: str = "lightgrey"


@dataclass
class HeatmapConfig:
    """Configuration for heatmap plots."""
    colormap: str = "viridis_r"
    granularity: float = 0.02


@dataclass
class DensityConfig:
    """Configuration for density plots."""
    colormap: str = "viridis"


@dataclass
class OutputConfig:
    """Configuration for output formats."""
    formats: List[str] = field(default_factory=lambda: ["png"])
    transparent: bool = False


@dataclass
class VisualizationConfig:
    """Configuration for visualization parameters."""
    plot_settings: PlotSettingsConfig = field(default_factory=PlotSettingsConfig)
    colors: ColorsConfig = field(default_factory=ColorsConfig)
    heatmap: HeatmapConfig = field(default_factory=HeatmapConfig)
    density: DensityConfig = field(default_factory=DensityConfig)
    output: OutputConfig = field(default_factory=OutputConfig)


@dataclass
class SimulationConfig:
    """Configuration for simulation parameters."""
    mode: str = "locus"
    
    # Locus mode parameters
    n_loci: Optional[int] = None
    locus_length: Optional[int] = None
    
    # Chromosome mode parameters
    chromosome_length: Optional[int] = None
    rec_rate: Optional[float] = None
    
    # Common parameters
    mutation_rate: float = 1e-8
    ploidy: int = 1
    
    # Demographic parameters
    populations: List[PopulationConfig] = field(default_factory=list)
    splits: List[SplitConfig] = field(default_factory=list)
    migration: Dict[str, float] = field(default_factory=dict)


@dataclass
class TwisstnternSimulateConfig:
    """Main configuration class for twisstntern_simulate."""
    
    # Core settings
    output_dir: str = "Results"
    verbose: bool = False
    quiet: bool = False
    seed: Optional[int] = None
    population_labels: Dict[str, str] = field(default_factory=dict)
    
    # Sub-configurations
    simulation: SimulationConfig = field(default_factory=SimulationConfig)
    analysis: AnalysisConfig = field(default_factory=AnalysisConfig)
    visualization: VisualizationConfig = field(default_factory=VisualizationConfig)
    
    # Hydra-specific
    hydra: Dict[str, Any] = field(default_factory=dict)
