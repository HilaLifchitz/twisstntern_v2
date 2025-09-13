"""
External dependencies for the core module.

This module provides access to external tools and libraries used by both
twisstntern and twisstntern_simulate packages.
"""

# Make twisst module available at package level
from . import twisst

__all__ = ["twisst"]