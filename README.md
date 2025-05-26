# Twisstntern v2

A Python package for ternary coordinate analysis and visualization.

## Recent Updates

### Code Refactoring
- **Modular Structure:** The codebase is now organized into three main modules:
  - `core.py`: Core mathematical and helper functions
  - `analysis.py`: Statistical and triangle-based analysis functions
  - `visualization.py`: All plotting and visualization functions
  - `__init__.py`: Exposes all main functions at the package level

### Enhanced Features
- **Flexible Data Loading:** The `dump_data` function now:
  - Automatically detects CSV file format and delimiter
  - Handles various CSV formats (comma, tab, semicolon, space-separated)
  - Provides detailed feedback about file reading process
  - Skips empty lines and invalid data
  - Normalizes data automatically

- **Improved User Interface:**
  - Interactive CSV file selection
  - Excludes system files and checkpoints
  - Clear granularity selection options
  - Detailed error reporting

## Installation

Required dependencies:
```bash
pip install numpy pandas matplotlib scipy sympy
```

## Usage

```python
import twisstntern

# Load and process data
data = twisstntern.dump_data("myfile.csv")

# Run analysis
results = twisstntern.triangles_analysis(data, "fine", "myfile")

# Visualize results
twisstntern.plot(data, 0.1, "myfile")
```

## Features

- Data loading and preprocessing for ternary coordinate data
- Statistical analysis of ternary coordinate distributions
- Visualization of ternary coordinate data with various granularity levels
- Fundamental asymmetry analysis
- Triangle-based analysis with customizable granularity

## Contributing

Feel free to submit issues and enhancement requests! 