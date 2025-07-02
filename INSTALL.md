# Installation Guide

This package has been modernized to support installation via `uv` and standard pip with `pyproject.toml`.

## Quick Installation

### Using uv (recommended)

```bash
# Install uv if you haven't already
curl -LsSf https://astral.sh/uv/install.sh | sh

# Install the package in development mode
uv pip install -e .

# Or install from a wheel
uv pip install .
```

### Using pip

```bash
# Install in development mode
pip install -e .

# Or install from wheel
pip install .
```

## Development Installation

For development with all optional dependencies:

```bash
# Using uv
uv pip install -e ".[dev,test,docs]"

# Using pip
pip install -e ".[dev,test,docs]"
```

## Available CLI Commands

After installation, the following commands are available:

### Main Commands
- `twisstntern` - Main analysis command
- `twisst-analyze` - Alias for main analysis
- `twisstntern-simulate` - Simulation package
- `twisst-simulate` - Alias for simulation

### Custom Commands (Placeholders)

The following commands are available as placeholders. To enable them:

1. Edit `twisstntern/cli/base.py`
2. Uncomment the desired command class
3. Register it in the CLIRegistry
4. Uncomment the corresponding entry in `pyproject.toml`
5. Reinstall the package

Available placeholders:
- `twisst-custom` - Custom analysis pipeline
- `twisst-batch` - Batch processing
- `twisst-compare` - Compare analysis results
- `twisst-plot` - Standalone plotting
- `twisst-convert` - Format conversion

### Example: Enabling a Custom Command

1. In `twisstntern/cli/base.py`, uncomment the `CustomCommand` class:
```python
class CustomCommand(BaseCommand):
    # ... implementation
```

2. Register it:
```python
CLIRegistry.register(CustomCommand())
```

3. In `pyproject.toml`, uncomment:
```toml
twisst-custom = "twisstntern.cli.custom:main"
```

4. Reinstall:
```bash
uv pip install -e .
```

## Development Workflow

```bash
# Clone the repository
git clone <repository-url>
cd twisstntern_v2

# Install development dependencies
uv pip install -e ".[dev]"

# Run tests
pytest

# Format code
black twisstntern/ twisstntern_simulate/

# Type checking
mypy twisstntern/

# Run the main command
twisstntern --help
```

## Package Structure

```
twisstntern_v2/
├── pyproject.toml          # Modern package configuration
├── twisstntern/            # Main package
│   ├── cli/               # CLI framework
│   │   ├── __init__.py
│   │   ├── base.py        # Base CLI classes
│   │   ├── custom.py      # Custom command entry point
│   │   ├── batch.py       # Batch command entry point
│   │   ├── compare.py     # Compare command entry point
│   │   ├── plot.py        # Plot command entry point
│   │   └── convert.py     # Convert command entry point
│   └── ...
└── twisstntern_simulate/   # Simulation package
    └── ...
```

## Legacy Files

The following files are now replaced by `pyproject.toml`:
- `setup.py` (kept for compatibility)
- `requirements.txt` (dependencies now in pyproject.toml)

You can safely remove these files if you want to fully modernize the package. 