# TWISSTNTERN

A package for analyzing ternary data from topology weights.

## Installation

1. Install the package and its dependencies:

```bash
pip install -r requirements.txt
pip install -e .
```

2. Or install directly with all dependencies:

```bash
pip install -e .[dev]
```

## Usage

### Command line interface:

```bash
# Analyze a tree file with custom granuality (0.25)
python -m twisstntern tree_file.trees 0.25

# granualities can be specified by words:

#                   coarse = 0.25
#                     fine = 0.1
#                superfine = 0.05

# default granuality is 0.1 
python -m twisstntern tree_file.trees 


# Analyze a CSV file  
python -m twisstntern weights.csv fine

# Analyze Newick trees with specific taxa
python -m twisstntern trees.newick 0.05 --taxon-names A B C D --outgroup A
```


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

Feel free to submit issues and enhancement requests! 

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