# TWISSTNTERN Hydra Usage Guide

## Overview

Both `twisstntern` and `twisstntern_simulate` now use a simplified Hydra configuration system with integrated simulation templates and automatic configuration logging.

## Basic Usage

### TWISSTNTERN Analysis

```bash
# Basic analysis with CSV data
python -m twisstntern file=data.csv

# With custom parameters
python -m twisstntern file=trees.newick granularity=0.05 output=MyAnalysis verbose=true

# With taxon names and outgroup
python -m twisstntern file=data.newick taxon_names=[A,B,C,D] outgroup=D verbose=true

# With downsampling
python -m twisstntern file=big_dataset.csv downsample=10 output=DownsampledResults
```

### TWISSTNTERN_SIMULATE (Using Embedded Configuration)

```bash
# Use embedded simulation config (no external file needed!)
python -m twisstntern_simulate

# With parameter overrides
python -m twisstntern_simulate simulation.seed=12345 output=Sim_Results verbose=true

# Override population parameters
python -m twisstntern_simulate simulation.populations.0.Ne=5000 simulation.migration.p1>p2=0.05

# Override simulation mode
python -m twisstntern_simulate simulation.simulation_mode=chromosome simulation.chromosome_length=1e7
```

### TWISSTNTERN_SIMULATE (Using External Config)

```bash
# Use external config file (old way still supported)
python -m twisstntern_simulate config_file=my_simulation.yaml
```

## Results and Configuration Logging

All runs automatically save:

1. **Analysis results**: PNGs, CSV files, plots in the output directory
2. **Final configuration**: `analysis_config.yaml` or `simulation_config.yaml` in results folder
3. **Log files**: Detailed execution logs

### Example Result Structure
```
Results/
├── analysis_config.yaml          # Final configuration used (with all overrides)
├── data_analysis.png              # Analysis plots
├── data_fundamental_asymmetry.png # Asymmetry analysis  
├── data_triangle_analysis.csv     # Raw results data
├── analysis.log                   # Execution log
└── ...                            # Other output files
```

## Hydra Sweeps (Parameter Studies)

Hydra makes it easy to run parameter sweeps for systematic studies.

### Example: Granularity Sweep

```bash
# Test different granularity values
python -m twisstntern file=data.csv granularity=0.05,0.1,0.15,0.2 output=GranularitySweep --multirun
```

This creates:
```
multirun/2024-01-01/12-00-00/
├── 0/                           # granularity=0.05
│   ├── analysis_config.yaml
│   └── [results...]
├── 1/                           # granularity=0.1  
│   ├── analysis_config.yaml
│   └── [results...]
├── 2/                           # granularity=0.15
│   └── [results...]
└── 3/                           # granularity=0.2
    └── [results...]
```

### Example: Simulation Parameter Sweep

```bash
# Sweep over migration rates
python -m twisstntern_simulate simulation.migration.p1>p2=0.01,0.02,0.05,0.1 --multirun

# Sweep over population sizes
python -m twisstntern_simulate simulation.populations.0.Ne=1000,5000,10000 simulation.populations.1.Ne=1000,5000,10000 --multirun

# Complex sweep: migration rate vs population size
python -m twisstntern_simulate \
  simulation.migration.p1>p2=0.01,0.05 \
  simulation.populations.0.Ne=1000,5000 \
  --multirun
```

### Advanced Sweeps with Hydra

```bash
# Sweep with custom output naming
python -m twisstntern \
  file=data.csv \
  granularity=0.05,0.1,0.2 \
  hydra.sweep.dir=GranularityStudy \
  hydra.sweep.subdir='gran_${granularity}' \
  --multirun

# Parallel execution (if you have joblib installed)
python -m twisstntern \
  file=data.csv \
  granularity=0.05,0.1,0.15,0.2 \
  hydra/launcher=joblib \
  --multirun
```

## Configuration File Examples

### Custom Analysis Configuration

Create `my_analysis.yaml`:
```yaml
# Custom analysis configuration
file: "data/population_trees.newick"
output: "PopulationAnalysis"
granularity: 0.05
verbose: true
taxon_names: ["Population_A", "Population_B", "Population_C", "Population_D"]
outgroup: "Population_D"
downsample: "5"
```

Run with: `python -m twisstntern --config-path=. --config-name=my_analysis`

### Custom Simulation Configuration

Create `custom_sim.yaml`:
```yaml
# Override defaults from base simulate_config.yaml
defaults:
  - simulate_config

# Custom simulation parameters
simulation:
  seed: 42
  n_loci: 5000
  populations:
    - name: "O"
      Ne: 2000
      sample_size: 15
    - name: "p1" 
      Ne: 1500
      sample_size: 15
  migration:
    "p1>p2": 0.02
    "p2>p1": 0.02

output: "CustomSimulation"
granularity: 0.08
verbose: true
```

Run with: `python -m twisstntern_simulate --config-path=. --config-name=custom_sim`

## Advanced Features

### Using Hydra Variable Interpolation

```yaml
# Use variables and interpolation
defaults:
  - simulate_config

base_ne: 1000
migration_rate: 0.01

simulation:
  populations:
    - name: "p1"
      Ne: ${base_ne}
    - name: "p2"  
      Ne: ${oc.decode:${base_ne}*2}  # 2000
  migration:
    "p1>p2": ${migration_rate}
    "p2>p1": ${migration_rate}
```

### Environment-Specific Configurations

```bash
# Development mode
python -m twisstntern file=test_data.csv verbose=true granularity=0.2

# Production mode  
python -m twisstntern file=production_data.csv verbose=false granularity=0.05 output=ProductionResults
```

## Quick Reference

### Key Parameters - TWISSTNTERN

| Parameter | Type | Default | Example |
|-----------|------|---------|---------|
| `file` | str | null | `file=data.csv` |
| `output` | str | "Results" | `output=MyResults` |
| `granularity` | float | 0.1 | `granularity=0.05` |
| `verbose` | bool | false | `verbose=true` |
| `taxon_names` | list | null | `taxon_names=[A,B,C,D]` |
| `downsample` | str | null | `downsample=10` |

### Key Parameters - TWISSTNTERN_SIMULATE

| Parameter | Type | Default | Example |
|-----------|------|---------|---------|
| `config_file` | str | null | `config_file=sim.yaml` |
| `simulation.seed` | int | 4576 | `simulation.seed=12345` |
| `simulation.n_loci` | int | 10000 | `simulation.n_loci=5000` |
| `simulation.migration.p1>p2` | float | 0.01 | `simulation.migration.p1>p2=0.05` |
| `granularity` | float | 0.1 | `granularity=0.08` |
| `output` | str | "Results" | `output=SimResults` |

## Troubleshooting

### Common Issues

1. **File not found**: Use absolute paths or ensure files exist relative to execution directory
2. **Configuration errors**: Check the saved `*_config.yaml` file in results to see final configuration
3. **Sweep issues**: Use `--multirun` flag for parameter sweeps

### Getting Help

```bash
# Show configuration options
python -m twisstntern --help
python -m twisstntern_simulate --help

# Show current configuration
python -m twisstntern --cfg job
python -m twisstntern_simulate --cfg job

# Validate configuration without running
python -m twisstntern file=data.csv --cfg job
```

## Migration from Old System

### Old Way (External config files required)
```bash
python -m twisstntern_simulate -c config_template.yaml -o Results
```

### New Way (Embedded config, cleaner interface)
```bash
python -m twisstntern_simulate output=Results  # Uses embedded config!
```

The new system is more flexible, logs configurations automatically, and supports powerful sweep functionality for parameter studies. 