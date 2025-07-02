# TWISSTNTERN Hydra Configuration Guide

This guide explains how to use the new Hydra-based configuration system for both `twisstntern` and `twisstntern_simulate` modules.

## Overview

Both modules have been refactored to use [Hydra](https://hydra.cc/) for configuration management instead of argparse. This provides:
- Hierarchical configuration management
- Easy parameter overrides from command line
- Configuration composition and templating
- Better parameter validation
- Consistent interface across modules

## Basic Usage

### TWISSTNTERN Analysis

```bash
# Basic usage with file parameter override
python -m twisstntern file=path/to/your/data.csv

# With additional parameters
python -m twisstntern file=data.csv output=MyResults granularity=0.05 verbose=true

# Using a specific configuration
python -m twisstntern --config-name=batch file=data.csv

# Using configuration from different path
python -m twisstntern --config-path=examples --config-name=basic_analysis
```

### TWISSTNTERN_SIMULATE

```bash
# Basic simulation usage
python -m twisstntern_simulate config_file=config_template.yaml

# With parameter overrides
python -m twisstntern_simulate config_file=config_template.yaml output=SimResults verbose=true seed=12345

# With advanced options
python -m twisstntern_simulate config_file=sim.yaml granularity=0.1 downsample=10 density_colormap=plasma
```

## Configuration Files

### Available CLI Configurations

The system includes several pre-configured CLI modes:

- `base`: Default configuration for standard analysis
- `batch`: Optimized for batch processing (verbose enabled)
- `compare`: For comparison analysis
- `plot`: For plotting/visualization focus
- `convert`: For file conversion operations
- `custom`: Template for custom configurations
- `simulate`: For simulation mode

### Configuration Structure

```yaml
# Example configuration file
file: "path/to/input.csv"           # Input file path
output: "Results"                   # Output directory
granularity: 0.1                    # Analysis granularity
verbose: true                       # Enable verbose logging
taxon_names: ["A", "B", "C", "D"]  # Taxon names (optional)
outgroup: "D"                       # Outgroup specification
topology_mapping: null             # Custom topology mapping
downsample: "10"                    # Downsampling configuration
```

### Simulation Configuration

```yaml
config_file: "simulation_config.yaml"  # Path to simulation config
output: "SimResults"                    # Output directory
verbose: false                          # Logging verbosity
quiet: false                            # Suppress output
seed: 12345                            # Random seed override
granularity: 0.1                       # Analysis granularity
topology_mapping: null                 # Custom topology mapping
override: []                           # Config value overrides
downsample: null                       # Tree/locus downsampling
downsampleKB: null                     # Kilobase downsampling
density_colormap: "viridis"            # Visualization colormap
```

## Command Line Overrides

Hydra allows easy parameter overrides from the command line:

```bash
# Override single parameters
python -m twisstntern file=data.csv granularity=0.05

# Override nested parameters (for simulation)
python -m twisstntern_simulate config_file=sim.yaml override=["migration.p1>p2=0.05"]

# Multiple overrides
python -m twisstntern file=data.csv output=Results verbose=true taxon_names=["Pop1","Pop2","Pop3","Pop4"]
```

## Advanced Usage

### Using Different Configuration Groups

```bash
# Use batch processing configuration
python -m twisstntern --config-name=config cli=batch file=data.csv

# Use plotting configuration
python -m twisstntern --config-name=config cli=plot file=data.csv
```

### Configuration Composition

You can create custom configuration files that extend base configurations:

```yaml
# custom_analysis.yaml
defaults:
  - cli: base
  - _self_

# Override specific parameters
granularity: 0.05
verbose: true
output: "CustomResults"
```

### Environment-Specific Configurations

```bash
# Development environment
python -m twisstntern file=test_data.csv verbose=true

# Production environment  
python -m twisstntern file=production_data.csv verbose=false output=ProductionResults
```

## Migration from Argparse

### Old CLI Usage (argparse)
```bash
# Old way
python -m twisstntern input_file.csv --granularity 0.1 --output Results --verbose
python -m twisstntern_simulate -c config.yaml -o Results --verbose --seed 123
```

### New CLI Usage (Hydra)
```bash
# New way
python -m twisstntern file=input_file.csv granularity=0.1 output=Results verbose=true
python -m twisstntern_simulate config_file=config.yaml output=Results verbose=true seed=123
```

## Parameter Reference

### TWISSTNTERN Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `file` | str | null | Input file path (required) |
| `output` | str | "Results" | Output directory |
| `granularity` | float/str | 0.1 | Analysis granularity |
| `taxon_names` | list | null | List of taxon names |
| `outgroup` | str | null | Outgroup taxon name |
| `topology_mapping` | str | null | Custom topology ordering |
| `downsample` | str | null | Downsampling format (N or N+i) |
| `verbose` | bool | false | Enable verbose logging |

### TWISSTNTERN_SIMULATE Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `config_file` | str | null | Simulation config file (required) |
| `output` | str | "Results" | Output directory |
| `verbose` | bool | false | Enable verbose output |
| `quiet` | bool | false | Suppress output |
| `seed` | int | null | Random seed override |
| `granularity` | float | 0.1 | Analysis granularity |
| `topology_mapping` | str | null | Custom topology mapping |
| `override` | list | [] | Config value overrides |
| `downsample` | str | null | Tree/locus downsampling |
| `downsampleKB` | str | null | Kilobase-based downsampling |
| `density_colormap` | str | "viridis" | Visualization colormap |

## Tips and Best Practices

1. **Use configuration files** for complex analyses with many parameters
2. **Use command line overrides** for quick parameter changes
3. **Leverage configuration groups** (batch, plot, etc.) for common use cases
4. **Create custom configurations** for repeated analysis workflows
5. **Use environment variables** with Hydra's variable interpolation for flexible deployments

## Troubleshooting

### Common Issues

1. **Missing required parameters**: Ensure `file` (for twisstntern) or `config_file` (for twisstntern_simulate) is specified
2. **Configuration file not found**: Check file paths and ensure config files exist
3. **Parameter type mismatches**: Hydra will validate parameter types based on defaults

### Getting Help

```bash
# Show available options for twisstntern
python -m twisstntern --help

# Show available options for twisstntern_simulate  
python -m twisstntern_simulate --help

# Show configuration structure
python -m twisstntern --cfg job
``` 