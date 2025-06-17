# TWISSTNTERN

<img src="images/logo.png" height="270pt" align="bottom">

## A method for analysing topology weights in a ternary framework

TWISSTNTERN is an analytical framework for visualising and analysing topology weights in a ternary plot. This is a modernized, command-line version that can calculate topology weights directly from tree files or analyze pre-computed weights from CSV files.

## What does TWISSTNTERN do?

In a tree with four populations, there are only 3 possible unrooted subtrees that can be observed for each sampled subtree. This makes a ternary plot a natural framework for analyzing the joint distribution of weights, as it is possible to graphically represent each genomic window as a single point based on the three topology weights.

The three corners of the ternary plot, [1,0,0], [0,1,0], [0,0,1], correspond with genomic windows that show taxon-level relationships that are consistent with one of the three possible subtrees; that is, 100% of the sampled subtrees perfectly match one of the three alternative trees, implying that samples from each of the four groups are monophyletic. In contrast, the very center of the ternary plot—[0.33,0.33,0.33]—corresponds with a genomic window where all three of the possible subtrees were sampled at equal frequency. Any other location in the ternary plot indicates a bias toward one of the subtrees, but with some resemblance to at least one of the other alternative topologies.

In an idealized four population model (3 splits with no migration), we expect the distribution of weights for many loci to be biased toward the top of the triangle which represents the subtree that matches the demographic history. Incomplete lineage sorting will generate a symmetrical distribution of weights between the left and right sides of the plot. This is because there is an equal chance that any discordant gene tree will more-closely resemble either alternative topology. However, processes like gene flow result in an asymmetrical pattern of discordance between the left and right halves of the triangle. The strength of the genealogical asymmetry (and its significance) can be quantified using the D_LR statistic, which is similar to the site-based statistic, Patterson's D. D_LR can be calculated at the genome-wide scale (two full half-triangles), or between left-right sub-triangles so that one can see how the strength of the asymmetry varies among loci that show different levels of discordance.

<img src="images/method_overview.png" height="450pt" align="bottom">

## Supported File Formats

**Tree Files:**
- **TreeSequence** (`.trees`, `.ts`): TSKit tree sequence files  
- **Newick** (`.newick`, `.nwk`, `.tree`): Single or multiple Newick format trees  
- **Nexus** (`.nexus`): Nexus format files

**Data Files:**
- **CSV** (`.csv`): Pre-computed topology weights (no normalization required)

---

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

## Usage

### 🔧 Command-Line Interface

```bash
python -m twisstntern [input_file] [granularity or keyword] [options]
```

**Input File Options:**
```bash
# Method 1: Positional argument
python -m twisstntern input_file.trees

# Method 2: Using -i or --input flag
python -m twisstntern -i input_file.trees
python -m twisstntern --input input_file.trees
```

**Examples:**

```bash
# Analyze a .trees file with default granularity (0.1)
python -m twisstntern tree_file.trees
# OR using flag syntax:
python -m twisstntern -i tree_file.trees

# Use a custom granularity value
python -m twisstntern tree_file.trees 0.25

# Analyze a Newick tree with specified taxon names and outgroup.
# Granularity can be specified using a keyword (e.g. coarse, fine, superfine)
python -m twisstntern -i tree_file.newick --granularity superfine --taxon-names O P1 P2 P3 --outgroup O

# Enable verbose logging for detailed output
python -m twisstntern tree_file.newick --taxon-names O P1 P2 P3 --outgroup O --verbose

# Specify an output directory (otherwise, 'Results/' will be created)
# * Note: files ending in `.trees` are msprime TreeSequence objects, which do not require specifying --taxon-names or --outgroup.
python -m twisstntern --input tree_file.tree --granularity 0.1 --taxon-names O P1 P2 P3 --outgroup O --output /your/custom/output_dir

# Use custom topology mapping to define T1, T2, T3 ordering
python -m twisstntern -i tree_file.trees --topology-mapping 'T1="(0,(3,(1,2)))"; T2="(0,(1,(2,3)))"; T3="(0,(2,(1,3)))";'

# Analyze a precomputed CSV file
python -m twisstntern weights.csv fine
```

---

### 🌳 Topology Mapping

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

#### **Benefits**
- **Consistent Analysis**: Ensure the same biological topology is always labeled as T1 across different datasets
- **Comparative Studies**: Compare results across analyses with consistent topology labels
- **Clear Visualization**: See ASCII trees showing exactly which topology corresponds to which label

#### **Output**
When topology mapping is applied, you'll see the ASCII tree visualizations printed to screen, e.g:
```
T1:
   /-0
--|
  |   /-3
   \-|
     |   /-1
      \-|
         \-2

T2:
   /-0
--|
  |   /-1
   \-|
     |   /-2
      \-|
         \-3

T3:
   /-0
--|
  |   /-2
   \-|
     |   /-1
      \-|
         \-3
```

---

### 🧪 Granularity Settings

Control the resolution of the triangle-based analysis using the `--granularity` argument. Supported predefined values:

| Keyword     | Value   |
|-------------|---------|
| `coarse`    | `0.25`  |
| `fine`      | `0.1`   |
| `superfine` | `0.05`  |

You can also provide a **custom float value** for granularity (alpha), e.g. `--granularity 0.01`.  
**Note:** `1 / alpha` must be an **even integer**.

- ✅ Valid example: `0.01` → `1 / 0.01 = 100`
- ❌ Invalid example: `0.2` → `1 / 0.2 = 5`

---

### 📝 Logging and Output

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

#### **Console vs File Output**
- **Console**: Clean, colored progress messages for real-time feedback
- **Log File**: Complete technical details with timestamps and module context, including detailed topology logging

#### **Topology Logging**
When processing tree files (Newick, TreeSequence), TWISSTNTERN automatically logs detailed topology information to the log file, including:

- **Topology Strings**: Simplified Newick format for each topology (e.g., `(O,((P1,P2),P3));`)
- **ASCII Tree Diagrams**: Beautiful visual representations of each topology structure
- **Topology Labels**: Clear identification of T1, T2, T3 assignments
- **Mapping Information**: When custom topology mapping is used, both original and reordered topologies are logged

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
```

This topology logging provides:
- **Reproducibility**: Exact topology structures used in analysis
- **Debugging**: Visual confirmation of topology assignments
- **Documentation**: Permanent record of tree structures for publications

#### **Benefits**
- **Full Traceability**: Complete record of every analysis for reproducibility
- **Debugging Support**: Detailed technical information when issues occur
- **Performance Monitoring**: Timing and file size tracking
- **Error Diagnosis**: Clear context when problems arise

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

## 🧬 Simulating Data with `twisstntern_simulate`

The `twisstntern_simulate` module allows you to generate simulated tree sequence data and analyze it using the same ternary pipeline as the main package. This is useful for testing, benchmarking, and exploring the behavior of the analysis under different demographic scenarios.

### 🔧 Command-Line Usage

In addition to the required configuration file (`--config`), you can optionally specify an output directory using `--output`. If no output directory is provided, a default `Results/` folder will be created.

You may also override individual parameters from the YAML config directly via additional command-line options (e.g., to test alternative values without editing the file), see below


```bash
python -m twisstntern_simulate -c CONFIG[-o OUTPUT][--skip-twisst-check][--force-download][--verbose][--quiet][--log-file LOG_FILE][--seed SEED][--mode {locus,chromosome}][--granularity GRANULARITY]
```

- `-c`, `--config`: **(Required)** Path to a YAML configuration file specifying simulation parameters (see `twisstntern_simulate/config_template.yaml` for an example).
- `-o`, `--output`: (Optional) Output directory for results. Defaults to `Results/`.
- `--skip-twisst-check`: Skip checking for the TWISST executable.
- `--force-download`: Force re-download of the TWISST executable.
- `--verbose`: Enable verbose logging.
- `--quiet`: Suppress most output.
- `--log-file LOG_FILE`: Write logs to a file.
- `--seed SEED`: Set a random seed for reproducibility.
- `--mode {locus,chromosome}`: Simulation mode (default: locus).
- `--granularity GRANULARITY`: Set analysis granularity (overrides config).

**Example:**
```bash
python -m twisstntern_simulate -c config_template.yaml -o SimResults/
```

### 📄 What It Does

- Runs msprime-based simulations according to your config.
- Saves simulated trees in Newick format (with consistent taxon naming for downstream analysis).
- Passes the simulated data directly to the ternary analysis pipeline.
- Outputs results (CSV, plots, etc.) in the specified output directory, formatted identically to the main package.

### ⚙️ Configuration

Simulation parameters (demography, sample sizes, sequence length, etc.) are specified in a YAML file.  
See `twisstntern_simulate/config_template.yaml` for a template and documentation of available options.

### 📝 Output

- All results (trees, weights, CSVs, plots) are saved in the output directory.

---

**Integration:**  
You can use the outputs from `twisstntern_simulate` directly with the main `twisstntern` analysis tools for further exploration.

---

## Citation

If you use TWISSTNTERN, please cite:

**Stankowski, S., Z. B. Zagrodzka, M. Garlovsky, A. Pal, D. Shipilina, D Garcia Castillo, T. Broquet, E. Leader, J. Reeve, K. Johannesson, A. M. Westram, R. K. Butlin. 2023. Selection on many loci drove the origin and spread of a key innovation. _bioRxiv_ doi: https://doi.org/10.1101/2023.02.13.528213**

Stankowski et al 2023 is where we first used the TWISSTNTERN method to study patterns of tree discordance in _Littorina_.

For the underlying topology weighting method, please also cite:
**Martin, S. H., & Van Belleghem, S. M. (2017). Exploring evolutionary relationships across the genome using topology weighting. _Genetics_, 206(1), 429-438. https://doi.org/10.1534/genetics.116.194720**

---

## Dependencies

- `numpy>=1.21.0`
- `pandas>=1.3.0`  
- `scipy>=1.7.0`
- `matplotlib>=3.4.0`
- `tskit>=0.4.0`
- `msprime>=1.0.0`
- `ete3>=3.1.0`
- `requests>=2.25.0`

---

## Features

- Ternary data loading and preprocessing
- Visualization with adjustable granularity
- Fundamental asymmetry detection
- Triangle-based statistical analysis
- TWISST-style topology weighting integration
- **Custom topology mapping** with beautiful ASCII tree visualization
- **Comprehensive logging** with automatic log file generation
- Tree sequence generation and demographic simulation via `msprime`

---

## Contributing

Feedback, issues, and enhancement suggestions are welcome!  
To install in development mode and test the CLI:

```bash
pip install -r requirements.txt
pip install -e .
python -m twisstntern --help
python -m twisstntern_simulate --help
```
