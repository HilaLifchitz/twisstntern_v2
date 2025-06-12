# TWISSTNTERN

A Python package for analyzing ternary data derived from topology weights in phylogenetic trees.

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

**Examples:**

```bash
# Analyze a .trees file with default granularity (0.1)
python -m twisstntern tree_file.trees

# Use a custom granularity value
python -m twisstntern tree_file.trees 0.25

# Analyze a Newick tree with specified taxon names and outgroup.
# Granularity can be specified using a keyword (e.g. coarse, fine, superfine)
python -m twisstntern tree_file.newick --granularity superfine --taxon-names O P1 P2 P3 --outgroup O

# Specify an output directory (otherwise, 'Results/' will be created)
python -m twisstntern tree_file.tree --granularity 0.1 --taxon-names O P1 P2 P3 --outgroup O --output /your/custom/output_dir

# Analyze a precomputed CSV file
python -m twisstntern weights.csv fine
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

### 🐍 Python Interface

You can also use TWISSTNTERN as a Python module:

```python
from twisstntern import run_analysis

results, fundamental_results, csv_file = run_analysis(
    file="your_file.trees",
    granularity=0.1
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
