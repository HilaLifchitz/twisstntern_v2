#!/usr/bin/env python
"""
Entry point for running twisstntern_simulate as a module.

This allows users to run:
    python -m twisstntern_simulate [args]
    
instead of:
    python -m twisstntern_simulate.main [args]
"""

from twisstntern_simulate.main import main

if __name__ == "__main__":
    main()
