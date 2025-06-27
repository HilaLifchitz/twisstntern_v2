# TWISSTNTERN

<img src="images/logo.png" height="270pt" align="bottom">

## Table of Contents
- [What does TWISSTNTERN do?](#what-does-twisstntern-do)
- [Installation](#installation)
- [Automatic twisst.py Handling](#automatic-twisstpy-handling)
- [Input](#input)
- [Command-Line Usage](#command-line-usage)
- [Topology Mapping](#topology-mapping)
- [Granularity Settings](#granularity-settings)
- [Output](#output)
- [Simulating Data with twisstntern_simulate](#simulating-data-with-twisstntern_simulate)
- [twisstntern_simulate Input](#input-1)
- [twisstntern_simulate Command-Line Usage](#command-line-usage-1)
- [twisstntern_simulate Configuration Overrides](#configuration-overrides)
- [twisstntern_simulate Output](#output-1)
- [Citation](#citation)
- [Dependencies](#dependencies)
- [Contributing](#contributing)

## A method for analysing topology weights in a ternary framework

TWISSTNTERN is an analytical framework for visualising and analysing topology weights in a ternary plot. This is a modernized, command-line version that can calculate topology weights directly from tree files or analyze pre-computed weights from CSV files.

## What does TWISSTNTERN do?

In a tree with four populations, there are only 3 possible unrooted subtrees that can be observed for each sampled subtree. This makes a ternary plot a natural framework for analyzing the joint distribution of weights, as it is possible to graphically represent each genomic window as a single point based on the three topology weights.

The three corners of the ternary plot, [1,0,0], [0,1,0], [0,0,1], correspond with genomic windows that show taxon-level relationships that are consistent with one of the three possible subtrees; that is, 100% of the sampled subtrees perfectly match one of the three alternative trees, implying that samples from each of the four groups are monophyletic. In contrast, the very center of the ternary plot—[0.33,0.33,0.33]—corresponds with a genomic window where all three of the possible subtrees were sampled at equal frequency. Any other location in the ternary plot indicates a bias toward one of the subtrees, but with some resemblance to at least one of the other alternative topologies.

In an idealized four population model (3 splits with no migration), we expect the distribution of weights for many loci to be biased toward the top of the triangle which represents the subtree that matches the demographic history. Incomplete lineage sorting will generate a symmetrical distribution of weights between the left and right sides of the plot. This is because there is an equal chance that any discordant gene tree will more-closely resemble either alternative topology. However, processes like gene flow result in an asymmetrical pattern of discordance between the left and right halves of the triangle. The strength of the genealogical asymmetry (and its significance) can be quantified using the D_LR statistic, which is similar to the site-based statistic, Patterson's D. D_LR can be calculated at the genome-wide scale (two full half-triangles), or between left-right sub-triangles so that one can see how the strength of the asymmetry varies among loci that show different levels of discordance.

<img src="images/method_overview.png" height="450pt" align="bottom">

## Installation

You can install the package in one of the following ways:

### 🛠️ Option 1: Install from GitHub

```bash
pip install git+https://github.com/HilaLifchitz/twisstntern_v2
```

### 🛠️ Option 2: Development Mode (recommended for contributors)

```bash
pip install -r requirements.txt
pip install -e .[dev]
```

---

### 🧬 Automatic twisst.py Handling

You **do not need to manually install or download `twisst.py`** to use this package!

Both the main analysis (`twisstntern`) and simulation (`twisstntern_simulate`) pipelines will:
- **Automatically check for the required `twisst.py` script** (used for topology weighting).
- **Download and patch it if missing or outdated**—no manual steps required.
- **Work out-of-the-box** in any fresh environment, as long as you install this package and its dependencies.

### How it works
- When you run:
  ```bash
  python -m twisstntern treeSfiles/CHROM.trees
  # or
  python -m twisstntern_simulate -c config_template.yaml
  ```
  The pipeline will ensure `twisst.py` is present and compatible before starting the analysis.

- If you want to **force a fresh download** (e.g., if you suspect a corrupted file), use:
  ```bash
  python -m twisstntern_simulate -c config_template.yaml --force-download
  ```

- If you want to **skip the check** (not recommended unless you know what you're doing), use:
  ```bash
  python -m twisstntern_simulate -c config_template.yaml --skip-twisst-check
  ```

### Why is this needed?
- `twisst.py` is a third-party script not available on PyPI.
- This package bundles and manages it for you, so you never have to worry about missing dependencies.

---

## Input

**Tree File Options:**

- **TreeSequence** (`.trees`, `.ts`): TSKit tree sequence files
- **Newick** (`.newick`, `.nwk`, `.tree`): Single or multiple Newick format trees
- **Nexus** (`.nexus`): Nexus format files

**Weights Data Files:**

- **CSV** (`.csv`): Pre-computed topology weights (no normalization required)

---

## 🔧 Command-Line Usage
**TWISSTNTERN** accepts both tree files and pre-computed topology weights. The input file can be specified either as a positional argument or using the `--input` flag.

```bash
python -m twisstntern INPUT [GRANULARITY] [OPTIONS]
```

#### **Parameters**

- `INPUT`: **(Required)** Input file path - tree file (`.trees`, `.newick`, `.nwk`, `.tree`, `.nexus`) or topology weights CSV `.csv`.
- `--granularity`: _(Optional, default: `0.1`)_  
  Sets the resolution of the ternary triangle analysis. Accepts either a float (e.g., `0.05`) or a keyword (`coarse`, `fine`, `superfine`).  
  Smaller values produce finer subdivisions of the triangle and more detailed results.
- `-i`, `--input`: Alternative way to specify input file
- `-o`, `--output`: Output directory (default: `Results/`)
- `--taxon-names`: Space-separated taxon names for Newick/Nexus files (e.g., `O P1 P2 P3`)
- `--outgroup`: Outgroup taxon name for Newick/Nexus files
- `--topology-mapping`: _(Optional)_  
  Manually specify which topology corresponds to each axis label (T1, T2, T3) in the ternary plot.  
  Useful for ensuring consistency across runs or datasets. Format:  
  `'T1="(0,(3,(1,2)))"; T2="(0,(1,(2,3)))"; T3="(0,(2,(1,3)))";'`
- `--downsample N` or `--downsample "N+i"`: Downsample the data by keeping only every Nth tree/locus. 
  - **Format**: `N` (every Nth starting from index 0) or `"N+i"` (every Nth starting from index i< N)
- `--verbose`: Enable verbose logging (DEBUG level)
- `--help`: Show help message

---

#### 🌳 Topology Mapping

**TWISSTNTERN** now supports **custom topology ordering** for phylogenetic analyses. By default, topologies are labeled T1, T2, T3 based on the order they are discovered by twisst. However, you can specify which topology should be assigned to each label using the `--topology-mapping` argument.

#### **Format**

```bash
--topology-mapping 'T1="(topology1)"; T2="(topology2)"; T3="(topology3)";'
```

#### **Examples**

**For TreeSequence files** (population IDs: 0, 1, 2, 3):

```bash
python -m twisstntern data.trees \
  --topology-mapping 'T1="(0,(3,(1,2)))"; T2="(0,(1,(2,3)))"; T3="(0,(2,(1,3)))";'
```

**For Newick files** (population names: O, P1, P2, P3):

```bash
python -m twisstntern data.newick \
  --taxon-names O P1 P2 P3 --outgroup O \
  --topology-mapping 'T1="(O,(P3,(P1,P2)))"; T2="(O,(P1,(P2,P3)))"; T3="(O,(P2,(P1,P3)))";'
```

#### 🧪 Granularity Settings

Control the resolution of the triangle-based analysis using the `--granularity` argument. Supported predefined values:

| Keyword     | Value  |
| ----------- | ------ |
| `coarse`    | `0.25` |
| `fine`      | `0.1`  |
| `superfine` | `0.05` |

You can also provide a **custom float value** for granularity (alpha), e.g. `--granularity 0.01`.  
**Note:** `1 / alpha` must be an **even integer**.

- ✅ Valid example: `0.01` → `1 / 0.01 = 100`
- ❌ Invalid example: `0.2` → `1 / 0.2 = 5`

---

## **Output**

TWISSTNTERN generates a comprehensive set of output files saved to the specified output directory (default: `Results/`):

**Analysis Files:**

- `[prefix]_topology_weights.csv` - Raw topology weight data for each genomic window. In chromosome mode, this file will always include a `position` column indicating the start position of each tree along the chromosome (in base pairs), which is used for KB-based downsampling and genome position-aware analyses.
- `[prefix]_triangle_analysis_[granularity].csv` - Triangle-based sub-analysis results with statistics
- `twisstntern_YYYYMMDD_HHMMSS.log` - Detailed log file with complete analysis record

**Visualization Files:**

- `[prefix]_fundamental_asymmetry.png` - Fundamental asymmetry bar chart showing left vs right bias
- `[prefix]_analysis_granularity_[value].png` - Ternary plot with data points colored by triangle regions
- `[prefix]_granuality_[value].png` - Main ternary plot with density visualization and statistical overlays
- `[prefix]_index_granularity_[value].png` - Triangle index visualization showing region boundaries
- `[prefix]_heatmap_count_granularity_[value].png` - Ternary heatmap showing data point density in each subtriangle

**File Naming:**

- `[prefix]` is derived from the input filename (e.g., `data.trees` → `data_`)
- `[value]` and `[granularity]` represent the granularity setting used (e.g., `0.1`, `0.05`, `superfine`, `fine`, `coarse`)
  - **Granularity values**: `superfine` = 0.05, `fine` = 0.1, `coarse` = 0.25, or custom float values
  - **All output files** now include the granularity value in their names for easy identification

**Example Output (granularity 0.1):**

```
Results/
├── data_topology_weights.csv
├── data_triangle_analysis_0.1.csv
├── data_fundamental_asymmetry.png
├── data_analysis_granularity_0.1.png
├── data_granuality_0.1.png
├── data_index_granularity_0.1.png
├── data_heatmap_count_granularity_0.1.png
└── twisstntern_20250618_151932.log
```

---

#### 📝 Logging

**TWISSTNTERN** includes comprehensive logging to track analysis progress and results. Every analysis automatically generates a detailed log file saved to the output directory.

#### **Log File Creation**

- **Automatic**: Log files are created for every analysis run
- **Location**: Saved in the output directory (default: `Results/`)
- **Format**: `twisstntern_YYYYMMDD_HHMMSS.log`

#### **Logging Levels**

```bash
# Standard logging (INFO level)
python -m twisstntern input.csv

# Verbose logging (DEBUG level) for detailed technical information
python -m twisstntern input.newick --taxon-names O P1 P2 P3 --outgroup O --verbose
```

#### **What Gets Logged**

- **System Information**: Python version, platform, package versions
- **Analysis Parameters**: Input file, granularity, taxon names, topology mapping
- **Processing Steps**: File format detection, tree processing, statistical analysis
- **Topology Information**: Complete topology details with both string representations and ASCII tree diagrams
- **Results Summary**: Fundamental asymmetry values, file generation, timing
- **Error Context**: Detailed error messages with helpful suggestions
<!--

#### **Console vs File Output**

- **Console**: Clean, colored progress messages for real-time feedback
- **Log File**: Complete technical details with timestamps and module context, including detailed topology logging -->

<!-- #### **Topology Logging**

When processing tree files (Newick, TreeSequence), TWISSTNTERN automatically logs detailed topology information to the log file, including:

- **Topology Strings**: Simplified Newick format for each topology (e.g., `(O,((P1,P2),P3));`)
- **ASCII Tree Diagrams**: Beautiful visual representations of each topology structure
- **Topology Labels**: Clear identification of T1, T2, T3 assignments -->

<!--
**Example Topology Log Content:**

```
2025-06-17 17:12:42,554 - INFO - Newick file topologies (default order)
2025-06-17 17:12:42,555 - INFO - ==================================================
2025-06-17 17:12:42,555 - INFO - T1:
2025-06-17 17:12:42,555 - INFO -   String: (O,((P1,P2),P3));
2025-06-17 17:12:42,555 - INFO -   ASCII Tree:
2025-06-17 17:12:42,555 - INFO -        /-O
2025-06-17 17:12:42,555 - INFO -       |
2025-06-17 17:12:42,555 - INFO -     --|      /-P1
2025-06-17 17:12:42,555 - INFO -       |   /-|
2025-06-17 17:12:42,555 - INFO -        \-|   \-P2
2025-06-17 17:12:42,555 - INFO -          |
2025-06-17 17:12:42,555 - INFO -           \-P3
2025-06-17 17:12:42,556 - INFO - T2:
2025-06-17 17:12:42,556 - INFO -   String: (O,((P1,P3),P2));
2025-06-17 17:12:42,556 - INFO -   ASCII Tree:
...
``` -->

---

### 🐍 Python Interface

You can also use TWISSTNTERN as a Python module:

```python
from twisstntern import run_analysis

# Basic usage
results, fundamental_results, csv_file = run_analysis(
    file="your_file.trees",
    granularity=0.1
)

# With topology mapping
results, fundamental_results, csv_file = run_analysis(
    file="your_file.trees",
    granularity=0.1,
    topology_mapping='T1="(0,(3,(1,2)))"; T2="(0,(1,(2,3)))"; T3="(0,(2,(1,3)))";'
)

# For Newick files
results, fundamental_results, csv_file = run_analysis(
    file="your_file.newick",
    granularity=0.1,
    taxon_names=["O", "P1", "P2", "P3"],
    outgroup="O",
    topology_mapping='T1="(O,(P3,(P1,P2)))"; T2="(O,(P1,(P2,P3)))"; T3="(O,(P2,(P1,P3)))";'
)
```

---

# 🧬 Simulating Data with `twisstntern_simulate`

The `twisstntern_simulate` module lets you generate simulated tree sequence data and analyze it using the same ternary pipeline as the main package.  
It is ideal for testing, benchmarking, and exploring different demographic scenarios.

- Runs simulations using `msprime`
- Automatically Downloads Twisst from https://github.com/simonhmartin/twisst
- Saves trees in Newick format with standardized taxon names
- Automatically passes simulated data to the analysis pipeline
- Outputs CSVs and plots in a structured results directory

---

## Input

Simulation parameters (demography, sample sizes, sequence length, etc.) are specified in a YAML file.  
See `config_template.yaml` for a template and documentation of available options.

---

## 🔧 Command-Line Usage

The only required input is a configuration file provided via the `--config` flag.

```bash
python -m twisstntern_simulate -c CONFIG [-o OUTPUT] [--downsample N] [--downsampleKB N] [--skip-twisst-check] [--force-download] [--verbose] [--quiet] [--log-file LOG_FILE] [--seed SEED] [--granularity GRANULARITY]
```

- `-c`, `--config`: **(Required)** Path to a YAML configuration file specifying simulation parameters (see `twisstntern_simulate/config_template.yaml` for an example).
- `-o`, `--output`: Output directory for results. Defaults to `Results/`.
- `--downsample N` or `--downsample "N+i"`: Downsample the data by keeping only every Nth tree/locus.Or every Nth starting from index i< N. Works for both locus and chromosome mode.
- `--downsampleKB Nkb` or `--downsampleKB Nkb+ikb`: **(chromosome mode only)** Downsample by keeping one tree every N kilobases along the simulated chromosome.
  - **Format**: `Nkb` (every N kb starting from position 0) or `Nkb+ikb` (every N kb starting from position i)
  - **Units supported**: kb (kilobases), mb (megabases), gb (gigabases), or bp (base pairs)
  - **Examples**: 
    - `--downsampleKB "100kb"` → every 100kb starting from position 0
    - `--downsampleKB "100kb+50kb"` → every 100kb starting from 50kb
    - `--downsampleKB "50kb+25000"` → every 50kb starting from 25,000 bp
  - **Constraint**: i < N (starting position must be less than sampling interval)
  - **This option is ignored in locus mode**
- `--granularity`: Set analysis granularity (default = 0.1).
- `--topology-mapping`: Manually specify which topology corresponds to each axis label (T1, T2, T3) in the ternary plot.
  Format: `'T1="(0,(3,(1,2)))"; T2="(0,(1,(2,3)))"; T3="(0,(2,(1,3)))";'`
  Useful for ensuring consistency across runs or datasets.
- `--override`: Override configuration values using the format 'key=value' or 'nested.key=value'. Examples:
  - `--override "migration.p2>p3=0.1"` → Set migration rate from population p2 to p3 to 0.1
  - `--override "populations.p1.Ne=2000"` → Set effective population size of p1 to 2000
  - `--override "samplesize=20"` → Set sample size to 20 for all non-ancestral populations
  - `--override "seed=12345"` → Set random seed to 12345
- `--skip-twisst-check`: Skip checking for the TWISST executable.
- `--force-download`: Force re-download of the TWISST executable.
- `--verbose`: Enable verbose logging.
- `--quiet`: Suppress most output.
- `--log-file LOG_FILE`: Write logs to a file.
- `--seed SEED`: Set a random seed for reproducibility.

**Examples:**

```bash
# Basic downsampling - every 10th locus/tree starting from index 0
python -m twisstntern_simulate -c config_template.yaml --downsample 10

# Enhanced downsampling - every 10th tree starting from index 3
python -m twisstntern_simulate -c config_template.yaml --downsample "10+3"

# KB-based downsampling - every 100kb starting from position 0
python -m twisstntern_simulate -c config_template.yaml --downsampleKB "100kb"

# Enhanced KB-based downsampling - every 100kb starting from 50kb
python -m twisstntern_simulate -c config_template.yaml --downsampleKB "100kb+50kb"

# KB-based downsampling with different units - every 1000kb starting from 30mb
python -m twisstntern_simulate -c config_template.yaml --downsampleKB "1000kb+30mb"

# Override sample size for all populations
python -m twisstntern_simulate -c config_template.yaml --override "samplesize=20"

# Multiple overrides including sample size
python -m twisstntern_simulate -c config_template.yaml --override "samplesize=15" --override "migration.p2>p3=0.1" --override "populations.p1.Ne=2000"
```

---

### ⚙️ Configuration Overrides

You can override any configuration parameter from the YAML file directly via the command line using the `--override` argument. This is extremely useful for testing different parameter values, running parameter sweeps, or exploring sensitivity analyses without editing configuration files.

#### **Format**

```bash
--override 'parameter_path=value'
```

#### **Supported Parameter Types**

**Top-level Parameters:**

```bash
--override "seed=1234"                    # Random seed (int)
--override "ploidy=2"                     # Ploidy level (int: 1=haploid, 2=diploid)
--override "mutation_rate=1e-7"           # Mutation rate (float, supports scientific notation)
--override "chromosome_length=5e7"        # Chromosome length (float)
--override "rec_rate=2e-8"               # Recombination rate (float)
--override "n_loci=500"                  # Number of loci for locus mode (int)
--override "locus_length=10000"          # Locus length for locus mode (int)
```

**Population Parameters** (format: `populations.{pop_name}.{parameter}=value`):

```bash
--override "populations.p1.Ne=2000"           # Effective population size (float)
--override "populations.p1.sample_size=15"    # Number of samples (int)
--override "populations.p2.growth_rate=0.01"  # Population growth rate (float)
```

**Migration Parameters** (format: `migration.{source}>{dest}=rate`):

```bash
--override "migration.p1>p2=0.05"        # Migration from p1 to p2 (float)
--override "migration.p2>p3=0.8"         # Migration from p2 to p3 (float)
--override "migration.O>p1=0.001"        # Migration from outgroup to p1 (float)
```

#### **Value Type Conversion**

The override system automatically converts values to appropriate types:

- **Integers**: `seed=1234` → `1234` (int)
- **Floats**: `Ne=1000.5` → `1000.5` (float)
- **Scientific notation**: `mutation_rate=1e-7` → `0.0000001`