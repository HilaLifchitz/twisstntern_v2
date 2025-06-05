# Twisstntern v2

A Python package for ternary coordinate analysis and visualization.

## Recent Updates

### Code Refactoring
- **Modular Structure:** The codebase is now organized into several main modules:
  - `core.py`: Core mathematical and helper functions
  - `analysis.py`: Statistical and triangle-based analysis functions
  - `visualization.py`: All plotting and visualization functions
  - `tree_processing.py`: Tree sequence processing and twisst integration
  - `simulation.py`: Demographic simulation and tree sequence generation
  - `__init__.py`: Exposes all main functions at the package level

### Enhanced Features
- **Flexible Data Loading:** The `dump_data` function now:
  - Automatically detects CSV file format and delimiter
  - Handles various CSV formats (comma, tab, semicolon, space-separated)
  - Provides detailed feedback about file reading process
  - Skips empty lines and invalid data
  - Normalizes data automatically

- **Tree Sequence Processing:**
  - Direct integration with twisst for topology analysis
  - Support for newick format tree sequences
  - Automatic conversion of topology weights to ternary coordinates
  - Efficient processing of large tree sequences

- **Demographic Simulation:**
  - Support for neutral demographic simulations
  - Configurable population sizes and migration rates
  - Flexible recombination and mutation rate settings
  - Output in both tree sequence and CSV formats

- **Improved User Interface:**
  - Interactive CSV file selection
  - Excludes system files and checkpoints
  - Clear granularity selection options
  - Detailed error reporting

## Installation

Required dependencies:
```bash
pip install numpy pandas matplotlib scipy sympy ete3 msprime requests
```

## Usage

### Basic CSV Input
```python
import twisstntern

# Load and process data
data = twisstntern.dump_data("myfile.csv")

# Run analysis
results = twisstntern.triangles_analysis(data, "fine", "myfile")

# Visualize results
twisstntern.plot(data, 0.1, "myfile")
```

### Tree Sequence Processing
```python
import twisstntern

# Define taxon groups
groups = {
    "pop1": ["sample1", "sample2"],
    "pop2": ["sample3", "sample4"],
    "pop3": ["sample5", "sample6"]
}

# Process tree sequence
ternary_data = twisstntern.process_tree_sequence(
    "trees.trees",
    groups,
    window_size=1000
)
```

### Demographic Simulation
```python
import twisstntern

# Create demographic model
model = twisstntern.DemographicModel("config.yaml")

# Simulate a single locus
ts = twisstntern.simulate_locus(
    model,
    n_ind=10,
    rec_rate=1e-8,
    mutation_rate=1e-8
)

# Simulate a chromosome
ts = twisstntern.simulate_chromosome(
    model,
    n_ind=10,
    n_loci=100,
    rec_rate=1e-8,
    mutation_rate=1e-8
)

# Save tree sequences
twisstntern.save_tree_sequences(ts, "output_dir")
```

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