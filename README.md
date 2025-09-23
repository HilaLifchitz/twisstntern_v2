# TWISSTNTERN

<img src="images/logo.png" height="270pt" align="bottom">

## Table of Contents
- [What does TWISSTNTERN do?](#what-does-twisstntern-do)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input](#input)
- [Command-Line Usage](#-command-line-usage)
- [Python Interface](#-python-interface)
- [Heatmap Colormap Customization](#-heatmap-colormap-customization)
- [Output](#output)
- [Simulating Data with twisstntern_simulate](#-simulating-data-with-twisstntern_simulate)
- [twisstntern_simulate Usage](#-twisstntern_simulate-command-line-usage)
- [twisstntern_simulate Python Interface](#-twisstntern_simulate-python-interface)
- [Advanced Features](#-advanced-features)
- [Citation](#citation)
- [Dependencies](#dependencies)
- [Contributing](#contributing)

---

## A method for analysing topology weights in a ternary framework

**TWISSTNTERN** is a modern analytical framework for visualizing and analyzing topology weights in ternary plots. This command-line package can calculate topology weights directly from tree files or analyze pre-computed weights from CSV files.

## What does TWISSTNTERN do?

In a tree with four populations, there are only 3 possible unrooted subtrees that can be observed for each sampled subtree. This makes a ternary plot a natural framework for analyzing the joint distribution of weights, as it is possible to graphically represent each genomic window as a single point based on the three topology weights.

The three corners of the ternary plot, [1,0,0], [0,1,0], [0,0,1], correspond with genomic windows that show taxon-level relationships that are consistent with one of the three possible subtrees; that is, 100% of the sampled subtrees perfectly match one of the three alternative trees, implying that samples from each of the four groups are monophyletic. In contrast, the very center of the ternary plot—[0.33,0.33,0.33]—corresponds with a genomic window where all three of the possible subtrees were sampled at equal frequency. Any other location in the ternary plot indicates a bias toward one of the subtrees, but with some resemblance to at least one of the other alternative topologies.

In an idealized four population model (3 splits with no migration), we expect the distribution of weights for many loci to be biased toward the top of the triangle which represents the subtree that matches the demographic history. Incomplete lineage sorting will generate a symmetrical distribution of weights between the left and right sides of the plot. This is because there is an equal chance that any discordant gene tree will more-closely resemble either alternative topology. However, processes like gene flow result in an asymmetrical pattern of discordance between the left and right halves of the triangle. The strength of the genealogical asymmetry (and its significance) can be quantified using the D_LR statistic, which is similar to the site-based statistic, Patterson's D. D_LR can be calculated at the genome-wide scale (two full half-triangles), or between left-right sub-triangles so that one can see how the strength of the asymmetry varies among loci that show different levels of discordance.

<img src="images/method_overview.png" height="450pt" align="bottom">

---

## Installation

### 🛠️ Option 1: Install from GitHub (Recommended)

```bash
pip install git+https://github.com/HilaLifchitz/twisstntern_v2
```

### 🛠️ Option 2: Development Mode (for contributors)

```bash
git clone https://github.com/HilaLifchitz/twisstntern_v2.git
cd twisstntern_v2
pip install -e .[dev]
```

After installation, you can use any of these equivalent commands:

```bash
# Download to current directory  
python -m src.twisstntern_simulate --get-config

# Download to specific location
python -m src.twisstntern_simulate --get-config /path/to/my_config.yaml
```

**Or use the Python interface:**

```python
from src.twisstntern_simulate.utils import download_config_template

# Download to current directory
config_path = download_config_template()

# Download to specific location  
config_path = download_config_template("my_simulation_config.yaml")
```

This automatically downloads the latest `config_template.yaml` from GitHub, ensuring you always have the most up-to-date configuration options.
# Main analysis commands
twisstntern your_data.csv              # Hydra CLI (recommended)
twisst-analyze your_data.csv           # Alias for Hydra CLI
python -m legacy.twisstntern your_data.csv    # Legacy argparse module

# Simulation commands  
twisstntern-simulate -c config.yaml    # Hydra CLI (recommended)
twisst-simulate -c config.yaml         # Alias for Hydra CLI
python -m legacy.twisstntern_simulate -c config.yaml  # Legacy argparse module
```

### 🧬 Automatic twisst.py Handling

You **do not need to manually install or download `twisst.py`**! Both packages (`twisstntern` and `twisstntern_simulate`) include the required `twisst.py` script automatically. No external dependencies or manual setup required.

---

## Quick Start

### 🚀 Analyze Existing Data

```bash
# Analyze CSV topology weights (creates timestamped Results_YYYY-MM-DD_HH-MM-SS/ directory)
twisstntern your_data.csv

# Analyze tree sequence file with custom output directory
twisstntern your_data.trees -o custom_analysis

# Analyze Newick trees with custom settings
twisstntern trees.newick --taxon-names O P1 P2 P3 --outgroup O --granularity fine
```

### 🧬 Simulate and Analyze Data

```bash
# Get configuration template (for external configs)
python -m src.twisstntern_simulate --get-config

# New Hydra approach - no config file needed
python -m src.twisstntern_simulate simulation=chromosome analysis.granularity=0.05

# Or use external config with overrides
python -m src.twisstntern_simulate config_file=config_template.yaml simulation.mode=chromosome
```

## 📦 **What Changed in v2**

### **🔄 Refactored `twisstntern_simulate`**

The simulation package has been completely refactored:

**Before (v1):**
- Argparse-based CLI with limited override options
- Single monolithic configuration approach
- Now preserved under `legacy/twisstntern_simulate/`

**After (v2):**
- **Hydra-based configuration system** with powerful composition
- **Modular config structure** in `src/twisstntern_simulate/configs/`
- **Advanced parameter overrides** without editing files
- **Same functionality, better interface**

### **🎯 Benefits of Hydra System**

- **Config composition**: `simulation=chromosome` switches entire mode
- **Deep overrides**: `simulation.populations[0].Ne=5000` 
- **Parameter sweeps**: Easy batch experiments with different parameters
- **Type safety**: Structured configs with validation
- **Backwards compatible**: Existing YAML files still work via `config_file=`

---
# TWISSTNTERN
## Input

### 📊 **Tree File Formats**
- **TreeSequence** (`.trees`, `.ts`): TSKit tree sequence files
- **Newick** (`.newick`, `.nwk`, `.tree`): Single or multiple Newick format trees  
- **Nexus** (`.nexus`): Nexus format files

### 📈 **Pre-computed Data**
- **CSV** (`.csv`): Topology weights data (T1, T2, T3 columns)

---

## 🔧 Command-Line Usage

```bash
twisstntern INPUT [OPTIONS]
```

### **Essential Parameters**

- `INPUT`: **(Required)** Input file path - tree file or CSV weights
- `-o`, `--output`: Output directory (default: auto-generated `Results_YYYY-MM-DD_HH-MM-SS/`)
- `--granularity`: Analysis resolution - `coarse` (0.25), `fine` (0.1), `superfine` (0.05), or custom float
- `--verbose`: Enable detailed logging


### **📊 Granularity Control**

| Keyword     | Value  | Use Case |
|-------------|--------|----------|
| `coarse`    | `0.25` | Quick overview, small datasets |
| `fine`      | `0.1`  | Standard analysis (default) |
| `superfine` | `0.05` | High-resolution, large datasets |
| Custom      | e.g., `0.02` | Specific research needs |

  
 **Constraint:** `1 / granularity` must be an **even integer**
- ✅ `--granularity 0.05` → 1 / 0.05 = 20 → even ✔️  
- ❌ `--granularity 0.2` → 1 / 0.2 = 5 → odd ✘  


### **🕐 Automatic Timestamped Output**

When no `-o` flag is provided, TWISSTNTERN automatically creates timestamped directories:
- **Format**: `Results_YYYY-MM-DD_HH-MM-SS/` (e.g., `Results_2025-07-03_14-30-25/`)
- **Benefits**: Prevents overwrites, organizes multiple runs, chronological sorting
- **Custom output**: Use `-o custom_name` to specify your own directory name

### **Tree-Specific Parameters**

- `--taxon-names`: Space-separated taxon names (e.g., `O P1 P2 P3`)
- `--outgroup`: Outgroup taxon name
- `--topology-mapping`: Custom topology assignment (see [Advanced Features](#-advanced-features))

### **Data Trimming**

- `--downsample N`: Keep every Nth data point
- `--downsample "N+i"`: Keep every Nth starting from index i, i<N

### **Examples**

```bash
# Basic analysis (creates Results_2025-07-03_14-30-25/)
twisstntern data.csv

# Custom output directory
twisstntern data.csv -o my_analysis

# Tree analysis with custom granularity
twisstntern trees.newick --taxon-names O P1 P2 P3 --outgroup O --granularity 0.05

# Downsampled analysis with verbose output
twisstntern large_dataset.csv --downsample "100+5" --verbose
```

---

## 🐍 Python Interface

### **Basic Usage**

```python
from legacy.twisstntern import run_analysis
```
This outputs:
 - results: pandas.DataFrame (grid of subtriangle analysis)
 - fundamental_results: tuple (n_right, n_left, D_LR, G_test, p_value)
 - csv_file: str (path to topology weights CSV)
)
```python
# Analyze CSV data
results, fundamental_results, csv_file = run_analysis(
    file="data.csv",
    granularity=0.1,
    output_dir="my_results",
)

# Analyze tree data
results, fundamental_results, csv_file = run_analysis(
    file="trees.newick",
    taxon_names=["O", "P1", "P2", "P3"],
    outgroup="O",
    granularity="superfine"
)
```

### **Advanced Usage**

```python
# Custom colormap and topology mapping
results, fundamental_results, csv_file = run_analysis(
    file="data.trees",
    granularity=0.1,
    colormap="plasma",
    topology_mapping='T1="(0,(3,(1,2)))"; T2="(0,(1,(2,3)))"; T3="(0,(2,(1,3)))";'
    downsample_N=100,
    downsample_i=20
)
```
---

## 🎨 Heatmap & radcount Colormap Customization

Customize ternary heatmap and radcount colors through the Python interface:

### **Available Colormaps**

| Colormap      |  Best For |
|---------------|---------|
| `"viridis_r"` |  Reversed viridis (default) | General use, colorblind-friendly |
| `"viridis"`   | Modern purple-to-yellow | Scientific publications |
| `"plasma"`    | Vibrant purple-pink-yellow | Eye-catching presentations |
| `"inferno"`   | Dark dramatic gradient | High contrast needed |
| `"Blues"`     |  Classic sequential blue | Publication-ready, minimal |
| `"Greys"`     |  Clean grayscale | Print-friendly, minimal |


**Note**: Colormap customization is available only through the Python interface to keep the CLI focused on core parameters.

---

## Output

### **Generated Files**

**📊 Analysis Data:**
- `[prefix]_topology_weights.csv` - Raw topology weight data
- `[prefix]_triangle_analysis_[granularity].csv` - Triangle-based statistics
- `twisstntern_YYYYMMDD_HHMMSS.log` - Detailed analysis log

**📈 Visualizations:**
- `[prefix]_fundamental_asymmetry.png` - Left vs right asymmetry bar chart
- `[prefix]_analysis_granularity_[value].png` - Ternary plot with triangle coloring
- `[prefix]_granuality_[value].png` - Main ternary plot with density visualization
- `[prefix]_index_granularity_[value].png` - Triangle index boundaries
- `[prefix]_heatmap.png` - **Ternary heatmap with customizable colormap**

### **File Naming Convention**

- `[prefix]` = Input filename (e.g., `data.trees` → `data_`)
- `[granularity]` = Analysis resolution (e.g., `0.1`, `fine`, `superfine`)

### **Example Output Structure**

```
Results_2025-07-03_14-30-25/
├── data_topology_weights.csv
├── data_triangle_analysis_0.1.csv
├── data_fundamental_asymmetry.png
├── data_analysis_granularity_0.1.png
├── data_granuality_0.1.png
├── data_index_granularity_0.1.png
├── data_radcount.png       # 🎨 Customizable colormap
├── data_heatmap.png        # 🎨 Customizable colormap
└── twisstntern_20250703_143025.log
```

---

# 🧬 Simulating Data with `twisstntern_simulate`

Generate and analyze simulated tree sequence data with demographic modeling using `msprime`. **Now uses Hydra configuration system for powerful parameter management.**

**Key Features:**
- 🔬 **Demographic modeling**: Population sizes, migration, growth rates
- 🧬 **Flexible simulation modes**: Independent loci or chromosome with recombination
- 🎯 **Hydra configuration**: Advanced parameter composition and override system
- 🎨 **Integrated analysis**: Automatic topology weight calculation and ternary plotting
- 📊 **Consistent output**: Same visualization pipeline as main package

## 🔧 twisstntern_simulate Command-Line Usage

### **⚡ New Hydra Configuration System**

The package now uses [Hydra](https://hydra.cc) for configuration management, providing powerful config composition and override capabilities.

```bash
# Basic usage with config composition
python -m src.twisstntern_simulate simulation=chromosome

# Parameter overrides (no config file needed)
python -m src.twisstntern_simulate simulation.n_loci=1000 seed=12345 verbose=true

# External config file with overrides
python -m src.twisstntern_simulate config_file=my_config.yaml simulation.mode=chromosome
```

### **🗂️ Configuration Structure**

Configs are stored in `src/twisstntern_simulate/configs/`:

```
configs/
├── config.yaml              # Main config with defaults
├── simulation/
│   ├── locus.yaml           # Locus mode parameters
│   └── chromosome.yaml      # Chromosome mode parameters
├── analysis/
│   └── default.yaml         # Analysis settings (granularity, downsampling)
└── visualization/
    └── default.yaml         # Plot settings and colormaps
```

### **🎛️ Parameter Override Examples**

```bash
# Simulation mode switching
python -m src.twisstntern_simulate simulation=chromosome
python -m src.twisstntern_simulate simulation=locus

# Population parameters
python -m src.twisstntern_simulate simulation.populations[0].Ne=10000
python -m src.twisstntern_simulate simulation.populations[0].sample_size=20

# Analysis parameters
python -m src.twisstntern_simulate analysis.granularity=0.05
python -m src.twisstntern_simulate analysis.downsampling.tree_interval=10

# Visualization settings
python -m src.twisstntern_simulate visualization.heatmap.colormap=plasma
python -m src.twisstntern_simulate visualization.colors.T1_color="#FF5733"

# Multiple overrides
python -m src.twisstntern_simulate \
  simulation=chromosome \
  simulation.chromosome_length=2000000 \
  analysis.granularity=0.02 \
  seed=12345 \
  verbose=true
```

### **📁 External Config Files**

Use existing YAML configs with Hydra overrides:

```bash
# Legacy config support
python -m src.twisstntern_simulate config_file=my_demographic_model.yaml

# With additional overrides
python -m src.twisstntern_simulate \
  config_file=my_config.yaml \
  simulation.mode=chromosome \
  seed=42
```

### **🔍 Configuration Help**

```bash
# View current config structure
python -m src.twisstntern_simulate --help

# Print resolved config without running
python -m src.twisstntern_simulate simulation=chromosome --cfg job
```

### **📂 Output Directory Management**

**Default behavior:**
- Outputs go to `Results/` directory (configurable via `output_dir=path`)
- Hydra does NOT change working directories (unlike default Hydra behavior)
- All file paths work as expected

```bash
# Custom output directory
python -m src.twisstntern_simulate output_dir=MyResults simulation=chromosome

# Outputs will be in: MyResults/
#   ├── locus_topology_weights.csv
#   ├── locus_trees.newick  
#   ├── locus_analysis_granularity_0.1.png
#   └── ... (same files as before)
```

**Hydra integration notes:**
- Configuration is managed by Hydra but outputs remain in your specified directory
- No `.hydra/` directories created in your workspace
- Same output structure as original version

---

## 🐍 twisstntern_simulate Python Interface

```python
from src.twisstntern_simulate.pipeline import run_pipeline
```
Like with twisstntern's run_analysis, returns:
 - results: pandas.DataFrame (triangle analysis)
 - fundamental_results: tuple (n_right, n_left, D_LR, G_test, p_value)
 - csv_file: str (path to topology weights CSV) -->

```python
# Basic simulation and analysis (same interface)
results, fundamental_results, csv_file = run_pipeline(
    config_path="config.yaml",
    output_dir="simulation_results",
    granularity=0.1
)

# With parameter overrides and custom colormap
results, fundamental_results, csv_file = run_pipeline(
    config_path="config.yaml",
    output_dir="simulation_results",
    granularity=0.1,
    seed_override=12345,
    config_overrides=["populations.p1.Ne=2000", "migration.p1>p2=0.05"],
    heatmap_colormap="plasma"
)

# Using Hydra structured configs
from src.twisstntern_simulate.hydra_config import TwisstnternSimulateConfig
from omegaconf import OmegaConf

# Create programmatic config
cfg = OmegaConf.structured(TwisstnternSimulateConfig)
cfg.simulation.mode = "chromosome"
cfg.simulation.chromosome_length = 2000000
cfg.analysis.granularity = 0.05
cfg.seed = 12345
```

---

## 🔬 Advanced Features

### **🌳 Topology Mapping**

 Defines how each topology is assigned to the axes of the ternary plot.  
  This allows you to explicitly control which topology appears on each vertex of the triangle. i.e. which topology appears as T1 (on top),which is T2 (bottom-left) and which T3 (bottom-right).


```bash
# For TreeSequence files (population IDs: 0,1,2,3)
twisstntern data.trees \
  --topology-mapping 'T1="(0,(3,(1,2)))"; T2="(0,(1,(2,3)))"; T3="(0,(2,(1,3)))";'

# For Newick files (population names: O,P1,P2,P3)  
twisstntern data.newick \
  --taxon-names O P1 P2 P3 --outgroup O \
  --topology-mapping 'T1="(O,(P3,(P1,P2)))"; T2="(O,(P1,(P2,P3)))"; T3="(O,(P2,(P1,P3)))";'
```


### **📝 Comprehensive Logging**

Every analysis generates detailed logs with:
- 🖥️ **System information**: Python version, platform, package versions
- ⚙️ **Parameters**: All settings used for the analysis
- 🌳 **Topology details**: ASCII tree diagrams and topology strings
- 📊 **Results summary**: Statistical outcomes and file generation
- 🕐 **Performance**: Timing and memory usage information

Access logs:
- **File**: `twisstntern_YYYYMMDD_HHMMSS.log` in output directory  
- **Console**: Use `--verbose` for real-time detailed output

---

## Axis Order for CSV Input (`--axis` flag)

By default, TWISSTNTERN expects CSV files to have columns in the order T1, T2, T3. If your CSV columns are in a different order, you can specify the mapping using the `--axis` flag.

**Usage:**

```
twisstntern data.csv --axis T2 T1 T3
```

This will interpret the first column as T2, the second as T1, and the third as T3. If you do not provide the flag, the default is `T1 T2 T3` (first column is T1, second is T2, third is T3).

**Example:**

If your CSV looks like this:

```
A,B,C
0.1,0.6,0.3
0.2,0.3,0.5
```

and you want to map A to T2, B to T1, and C to T3, run:

```
twisstntern mydata.csv --axis T2 T1 T3
```

## Output Directory (`--output` flag)

You can control where TWISSTNTERN writes all results and plots using the `--output` (or `-o`) flag. By default, if you do not specify this flag, a timestamped directory (e.g., `Results_2025-07-03_17-11-06`) will be created for each run. If you specify a directory name, all output files will be placed there. The directory will be created automatically if it does not exist.

**Usage:**

```
twisstntern data.csv --output my_results_folder
```

This will save all results and plots in the `my_results_folder` directory.

If you do not provide `--output`, a new directory named like `Results_YYYY-MM-DD_HH-MM-SS` will be created for you.


## Citation

### 📚 Citation

The **TwisstNTern** method was first introduced and applied in:

> **Stankowski et al. (2023)**  
> Stankowski, S., Zagrodzka, Z. B., Garlovsky, M., Pal, A., Shipilina, D., Garcia Castillo, D., Lifchitz, H. Broquet, T., Leader, E., Reeve, J., Johannesson, K., Westram, A. M., & Butlin, R. K. (2023).  
> *The genetic basis of a recent transition to live-bearing in marine snails.*  
> **Science**. https://doi.org/10.1126/science.adi2982

If you use **TwisstNTern**, please cite the article above.

---

### 🧰 Underlying Tools

**TwisstNTern** builds upon the **TWISST** framework for topology weighting.  
If you use this method, please also cite:

> **Martin & Van Belleghem (2017)**  
> Martin, S. H., & Van Belleghem, S. M. (2017).  
> *Exploring evolutionary relationships across the genome using topology weighting.*  
> **Genetics**, 206(1), 429–438. https://doi.org/10.1534/genetics.116.194720

And acknowledge the **TWISST** software available at:  
👉 https://github.com/simonhmartin/twisst

---

## Dependencies

### **Core Requirements**
- Python ≥ 3.8
- NumPy ≥ 1.21.0
- Pandas ≥ 1.3.0
- SciPy ≥ 1.7.0
- Matplotlib ≥ 3.4.0
- scikit-learn ≥ 1.0.0

### **Tree Processing**
- tskit ≥ 0.4.0
- msprime ≥ 1.0.0 *(for simulations)*
- ete3 ≥ 3.1.0

### **Configuration & Data**
- PyYAML ≥ 6.0.0
- requests ≥ 2.25.0
- hydra-core ≥ 1.3.0 *(for twisstntern_simulate)*
- omegaconf ≥ 2.3.0 *(for twisstntern_simulate)*

**Note**: `twisst.py` is included automatically - no manual installation required.

---

## Contributing

We welcome contributions! Please:

1. **Fork** the repository
2. **Create** a feature branch (`git checkout -b feature/amazing-feature`)
3. **Commit** your changes (`git commit -m 'Add amazing feature'`)
4. **Push** to the branch (`git push origin feature/amazing-feature`)
5. **Open** a Pull Request

### **Development Setup**

```bash
git clone https://github.com/HilaLifchitz/twisstntern_v2.git
cd twisstntern_v2
pip install -r requirements.txt
pip install -e .[dev]

# Run tests
pytest

# Format code
black .
```

---

<div align="center">

**🧬 Happy analyzing with TWISSTNTERN! 🔬**

*For issues, feature requests, or questions, please visit our [GitHub repository](https://github.com/HilaLifchitz/twisstntern_v2).*

</div>
