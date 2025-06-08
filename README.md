# TWISSTNTERN

A package for analyzing ternary data from topology weights.

## Supported File Formats

**Tree Files:**
- TreeSequence (.trees, .ts): TSKit tree sequence files
- Newick (.newick, .nwk, .tree): Single or multiple Newick format trees  
- Nexus (.nexus): Nexus format files

**Data Files:**
- CSV (.csv): Pre-computed topology weights (normalization not required)

## Installation

1. Install the package and its dependencies:

```bash
pip install git+https://github.com/HilaLifchitz/twisstntern_v2
pip install -r requirements.txt
pip install -e .
```

2. Or install directly with all dependencies:

```bash
pip install -e .[dev]
```

## Usage

### Command-Line Interface

granularities can be specified with the following keywords:

                   coarse = 0.25
                     fine = 0.1
                superfine = 0.05

```bash
# Analyze a .trees file with custom granularity (0.25)
# Taxon names and outgroup are automatically detected
python -m twisstntern tree_file.trees 0.25

# default granularity is 0.1 
python -m twisstntern tree_file.trees 

# Analyze Newick files — requires --taxon-names and --outgroup
# Also explicitly provide --granularity
python -m twisstntern phylogeny.tree --granularity 0.1 --taxon-names O P1 P2 P3 --outgroup O
python -m twisstntern trees.newick --granularity superfine --taxon-names O P1 P2 P3 --outgroup O

# Analyze a CSV file  
python -m twisstntern weights.csv fine



### Python interface:

```python
from twisstntern import run_analysis

# Analyze any supported file
results, fundamental_results, csv_file = run_analysis(
    file="your_file.trees",
    granularity=0.1
)
```

## Dependencies

- numpy>=1.21.0
- pandas>=1.3.0  
- scipy>=1.7.0
- matplotlib>=3.4.0
- tskit>=0.4.0
- msprime>=1.0.0
- ete3>=3.1.0
- requests>=2.25.0

## Features

- Data loading and preprocessing for ternary coordinate data
- Statistical analysis of ternary coordinate distributions
- Visualization of ternary coordinate data with various granularity levels
- Fundamental asymmetry analysis
- Triangle-based analysis with customizable granularity
- Tree sequence processing with twisst integration
- Demographic simulation and tree sequence generation

## Contributing

Feel free to submit issues or enhancement requests! 

## Installation Instructions

1. Install dependencies

```bash
pip install -r requirements.txt
```

2. Install the package in development mode

```bash
pip install -e .
```

3. Test it works

```bash
python -m twisstntern --help
``` 