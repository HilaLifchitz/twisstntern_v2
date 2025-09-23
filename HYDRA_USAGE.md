# TWISSTNTERN Hydra Configuration Guide

Simple guide for configuring and running TWISSTNTERN packages with Hydra.

This guide covers both:
- **TWISSTNTERN** - Analysis of existing tree data (CSV files, Newick trees)
- **TWISSTNTERN Simulate** - Demographic simulation and tree sequence generation

## Quick Start

### TWISSTNTERN (Analysis)

**Single Analysis:**
```bash
conda activate tnt3
python -m src.twisstntern.main file=DATA.csv analysis.granularity=VALUE
```

**Parameter Sweep:**
```bash
python -m src.twisstntern.main --multirun hydra/launcher=joblib analysis.granularity=VALUES file=DATA.csv
```

### TWISSTNTERN Simulate (Simulation)

**Single Simulation:**
```bash
conda activate tnt3
python -m src.twisstntern_simulate.main analysis.granularity=VALUE simulation.n_loci=N
```

**Parameter Sweep:**
```bash
python -m src.twisstntern_simulate.main --multirun hydra/launcher=joblib analysis.granularity=VALUES simulation.n_loci=N
```

## Configuration System

Both packages use **unified configuration files** that handle single runs and parameter sweeps automatically:
- **TWISSTNTERN**: `src/twisstntern/configs/default.yaml`
- **TWISSTNTERN Simulate**: `src/twisstntern_simulate/configs/config.yaml`

## Usage Examples

### TWISSTNTERN (Analysis)

```bash
# Single analysis (your original command still works)
python -m src.twisstntern.main file="./examples/data_files/Littorina_data.csv" analysis.granularity=0.05

# Parameter sweep (2 values, all CPU cores)
python -m src.twisstntern.main --multirun hydra/launcher=joblib analysis.granularity=0.05,0.1 file="./examples/data_files/Littorina_data.csv"

# Comprehensive sweep (8 values)
python -m src.twisstntern.main --multirun hydra/launcher=joblib analysis.granularity=0.01,0.025,0.05,0.075,0.1,0.15,0.2,0.25 file="./examples/data_files/Littorina_data.csv"
```

### TWISSTNTERN Simulate (Simulation)

```bash
# Single simulation
python -m src.twisstntern_simulate.main analysis.granularity=0.1 simulation.n_loci=1000

# Parameter sweep (2 granularity values)
python -m src.twisstntern_simulate.main --multirun hydra/launcher=joblib analysis.granularity=0.05,0.1 simulation.n_loci=100

# Multi-parameter sweep (granularity + simulation mode)
python -m src.twisstntern_simulate.main --multirun hydra/launcher=joblib analysis.granularity=0.05,0.1 simulation=locus,chromosome simulation.n_loci=50

# Chromosome mode simulation
python -m src.twisstntern_simulate.main simulation=chromosome analysis.granularity=0.1 simulation.chromosome_length=1000000
```

#### Mode-Specific Sweeps

**Locus sweep (granularity × n_loci)** – explores how changing locus counts affects topology weights while keeping locus length fixed. Results are grouped under `outputs/simulate_sweeps/locus/` with the Hydra configs that produced them.
```bash
python -m src.twisstntern_simulate.main --multirun \
  hydra/launcher=joblib hydra.launcher.n_jobs=4 \
  hydra.sweep.dir=./outputs/simulate_sweeps \
  'hydra.sweep.subdir=locus/granularity_${analysis.granularity}' \
  simulation=locus \
  analysis.granularity=0.05,0.1 \
  simulation.n_loci=1000,5000 \
  simulation.locus_length=5000 \
  analysis.downsampling.tree_interval=10
```

**Chromosome sweep (granularity × chromosome length)** – runs the recombination-aware mode with longer chromosomes to reveal how genome length influences recovered trees. Outputs land in `outputs/simulate_sweeps/chromosome/` alongside per-run Hydra metadata.
```bash
python -m src.twisstntern_simulate.main --multirun \
  hydra/launcher=joblib hydra.launcher.n_jobs=4 \
  hydra.sweep.dir=./outputs/simulate_sweeps \
  'hydra.sweep.subdir=chromosome/granularity_${analysis.granularity}' \
  simulation=chromosome \
  analysis.granularity=0.05,0.1 \
  simulation.chromosome_length=1000000 \
  analysis.downsampling.tree_interval=10
```

### Resource Management

**For TWISSTNTERN (Analysis):**
```bash
# Limit to 4 cores (memory-constrained systems)
python -m src.twisstntern.main --multirun hydra/launcher=joblib hydra.launcher.n_jobs=4 analysis.granularity=0.01,0.05,0.1,0.2 file="./examples/data_files/Littorina_data.csv"

# Use different backend with 2 cores
python -m src.twisstntern.main --multirun hydra/launcher=joblib hydra.launcher.n_jobs=2 hydra.launcher.backend=multiprocessing analysis.granularity=0.05,0.1 file="./examples/data_files/Littorina_data.csv"

# Enable verbose output for debugging
python -m src.twisstntern.main --multirun hydra/launcher=joblib hydra.launcher.verbose=10 analysis.granularity=0.05,0.1 file="./examples/data_files/Littorina_data.csv"
```

**For TWISSTNTERN Simulate:**
```bash
# Limit to 4 cores for simulation sweeps
python -m src.twisstntern_simulate.main --multirun hydra/launcher=joblib hydra.launcher.n_jobs=4 analysis.granularity=0.05,0.1 simulation.n_loci=100

# Use multiprocessing backend for stability
python -m src.twisstntern_simulate.main --multirun hydra/launcher=joblib hydra.launcher.backend=multiprocessing analysis.granularity=0.05,0.1 simulation.n_loci=50
```

### Custom Output Directories
```bash
# Single run with custom output
python -m src.twisstntern.main file="./examples/data_files/Littorina_data.csv" analysis.granularity=0.1 output.output_dir=./my_results

# Parameter sweep with custom output
python -m src.twisstntern.main --multirun hydra/launcher=joblib analysis.granularity=0.05,0.1,0.2 output.output_dir="./my_sweep" file="./examples/data_files/Littorina_data.csv"
```

### Multi-Parameter Sweeps
```bash
# Sweep multiple parameters simultaneously (creates 3×2 = 6 combinations)
python -m src.twisstntern.main --multirun hydra/launcher=joblib analysis.granularity=0.05,0.1,0.2 visualization.style="RdBu_r","viridis" file="./examples/data_files/Littorina_data.csv"
```

## Configuration Options

### Launcher Settings
| Parameter | Default | Override | Description |
|-----------|---------|----------|-------------|
| `n_jobs` | -1 (all cores) | `hydra.launcher.n_jobs=4` | Number of parallel jobs |
| `backend` | "loky" | `hydra.launcher.backend=multiprocessing` | Parallel backend ("loky", "multiprocessing", "threading") |
| `verbose` | 0 (quiet) | `hydra.launcher.verbose=10` | Debug output level (0-10) |

### Analysis Parameters

**Common to both packages:**
| Parameter | Default | Override | Description |
|-----------|---------|----------|-------------|
| `granularity` | 0.1 | `analysis.granularity=0.05` | Main analysis granularity |

**TWISSTNTERN specific:**
| Parameter | Default | Override | Description |
|-----------|---------|----------|-------------|
| `heatmap_granularity` | 0.02 | `analysis.heatmap_granularity=0.01` | Fixed granularity for heatmaps |

**TWISSTNTERN Simulate specific:**
| Parameter | Default | Override | Description |
|-----------|---------|----------|-------------|
| `n_loci` | 10000 | `simulation.n_loci=1000` | Number of loci (locus mode) |
| `locus_length` | 10000 | `simulation.locus_length=5000` | Length of each locus |
| `chromosome_length` | - | `simulation.chromosome_length=1000000` | Chromosome length (chromosome mode) |
| `mutation_rate` | 1e-8 | `simulation.mutation_rate=2e-8` | Mutation rate per site per generation |

### Common Overrides

**TWISSTNTERN (Analysis):**
```bash
# Analysis parameters
analysis.granularity=0.05
analysis.heatmap_granularity=0.01

# Output settings
output.output_dir=./custom_results
output.verbose=true

# Visualization
visualization.style=viridis
visualization.dpi=600
```

**TWISSTNTERN Simulate:**
```bash
# Simulation parameters
simulation.n_loci=1000
simulation.locus_length=5000
simulation.mutation_rate=2e-8

# Analysis parameters
analysis.granularity=0.05

# Output settings
output_dir=./custom_results
verbose=true

# Switch simulation modes
simulation=chromosome  # or simulation=locus
```

## Output Structure

**Single run (both packages):**
```
outputs/2025-09-09/07-18-49/
├── *.png (plots)
├── *.csv (results)
├── *.log (logs)
└── .hydra/ (config)
```

**Parameter sweep (both packages):**
```
outputs/sweep/2025-09-09/07-18-57/
├── granularity_0.05/
│   ├── *.png (plots)
│   ├── *.csv (results)
│   └── *.log (logs)
├── granularity_0.1/
│   ├── *.png (plots)
│   ├── *.csv (results)
│   └── *.log (logs)
└── multirun.yaml (sweep metadata)
```

**TWISSTNTERN Simulate also generates:**
```
outputs/2025-09-09/07-18-49/
├── locus_trees.newick (or chromosome_trees.newick)
├── locus_topology_weights.csv
└── ... (analysis outputs)
```

## Performance

- **Sequential execution**: ~2s per granularity value
- **Parallel execution**: All values simultaneously (~3x faster)
- **Resource usage**: Scales with CPU cores (monitor with `htop`)

## Best Practices

1. **Start small**: Test with 2-3 parameter values first
2. **Monitor resources**: Use `htop` during execution
3. **Limit cores**: Use `hydra.launcher.n_jobs=4` on memory-constrained systems
4. **Debug issues**: Enable verbose output with `hydra.launcher.verbose=10`
5. **File paths**: Use absolute paths or ensure correct working directory

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Memory errors | Reduce `hydra.launcher.n_jobs` |
| Process hangs | Try `hydra.launcher.backend=multiprocessing` |
| Too many open files | Reduce parallel jobs |
| Slow execution | Check if parallel execution is enabled |

## File Structure

Your configuration system uses minimal files:

```
src/twisstntern/configs/
├── default.yaml                    # Single config for everything
├── analysis/                       # Analysis presets
├── visualization/                   # Visualization presets
└── hydra/launcher/
    └── joblib.yaml                 # Single launcher config
```

## Key Benefits

- ✅ **One config file**: No complex multi-file setup
- ✅ **CLI configurable**: Override any parameter from command line
- ✅ **True parallel**: Uses all CPU cores by default
- ✅ **Backward compatible**: Original commands still work
- ✅ **Pure Hydra**: No external scripts or dependencies

---

**Most common commands:**

**TWISSTNTERN (Analysis):**
```bash
# Single analysis
python -m src.twisstntern.main file="data.csv" analysis.granularity=0.05

# Parameter sweep (all cores)
python -m src.twisstntern.main --multirun hydra/launcher=joblib analysis.granularity=0.01,0.05,0.1,0.2 file="data.csv"

# Parameter sweep (limited cores)
python -m src.twisstntern.main --multirun hydra/launcher=joblib hydra.launcher.n_jobs=4 analysis.granularity=0.01,0.05,0.1,0.2 file="data.csv"
```

**TWISSTNTERN Simulate:**
```bash
# Single simulation
python -m src.twisstntern_simulate.main analysis.granularity=0.1 simulation.n_loci=1000

# Parameter sweep (all cores)
python -m src.twisstntern_simulate.main --multirun hydra/launcher=joblib analysis.granularity=0.05,0.1,0.2 simulation.n_loci=100

# Parameter sweep (limited cores)
python -m src.twisstntern_simulate.main --multirun hydra/launcher=joblib hydra.launcher.n_jobs=4 analysis.granularity=0.05,0.1 simulation.n_loci=50
```
