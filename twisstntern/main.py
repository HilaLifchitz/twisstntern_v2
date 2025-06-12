#!/usr/bin/env python
# coding: utf-8

"""
Main entry point for TWISSTNTERN analysis pipeline.
This file allows running the package directly with: python main.py
"""

import sys
from pathlib import Path

# Add the package directory to the Python path
package_dir = Path(__file__).parent / "twisstntern"
sys.path.insert(0, str(package_dir.parent))

# Import and run the main function from the package
from twisstntern.__main__ import main

if __name__ == "__main__":
    main() 