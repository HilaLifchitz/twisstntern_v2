#!/usr/bin/env python3
"""
Debug Environment Script
Run this to diagnose environment issues with twisstntern_simulate
"""

import sys
import os
from pathlib import Path

print("=" * 60)
print("TWISSTNTERN ENVIRONMENT DEBUGGING")
print("=" * 60)

# Basic environment info
print("\n1. BASIC ENVIRONMENT:")
print(f"   Current directory: {os.getcwd()}")
print(f"   Python executable: {sys.executable}")
print(f"   Python version: {sys.version}")
print(f"   CONDA_DEFAULT_ENV: {os.environ.get('CONDA_DEFAULT_ENV', 'NOT SET')}")

# Check if we're in the right directory
print("\n2. DIRECTORY CHECK:")
current_dir = Path.cwd()
project_files = ['setup.py', 'twisstntern_simulate', 'config_template.yaml']
missing_files = []
for file in project_files:
    if (current_dir / file).exists():
        print(f"   ✓ {file} found")
    else:
        print(f"   ✗ {file} MISSING")
        missing_files.append(file)

if missing_files:
    print(f"\n   WARNING: You might not be in the twissting_baby directory!")
    print(f"   Missing files: {missing_files}")
    print(f"   Try: cd /home/hlifchit/projects/twissting_baby")

# Check Python path
print("\n3. PYTHON PATH:")
for i, path in enumerate(sys.path[:5]):
    print(f"   [{i}] {path}")
if len(sys.path) > 5:
    print(f"   ... and {len(sys.path) - 5} more paths")

# Check imports
print("\n4. IMPORT TESTS:")
import_tests = [
    ("pandas", "pandas"),
    ("numpy", "numpy"),
    ("msprime", "msprime"),
    ("twisstntern.utils", "twisstntern main package"),
    ("twisstntern_simulate", "twisstntern_simulate package"),
]

for module, description in import_tests:
    try:
        __import__(module)
        print(f"   ✓ {description}")
    except ImportError as e:
        print(f"   ✗ {description} FAILED: {e}")

# Check package installation
print("\n5. PACKAGE INSTALLATION:")
try:
    import pkg_resources
    pkg = pkg_resources.get_distribution("twisstntern")
    print(f"   ✓ twisstntern version: {pkg.version}")
    print(f"   ✓ installed at: {pkg.location}")
except Exception as e:
    print(f"   ✗ twisstntern not properly installed: {e}")

# Test specific problematic import
print("\n6. SPECIFIC IMPORT TEST:")
try:
    from twisstntern_simulate.simulation import run_simulation
    print("   ✓ twisstntern_simulate.simulation imports successfully")
except ImportError as e:
    print(f"   ✗ twisstntern_simulate.simulation import FAILED: {e}")
    print("   This is likely the source of your problem!")

# Check if we can run the module
print("\n7. MODULE EXECUTION TEST:")
try:
    import subprocess
    result = subprocess.run([sys.executable, "-m", "twisstntern_simulate", "--help"], 
                          capture_output=True, text=True, timeout=10)
    if result.returncode == 0:
        print("   ✓ Module execution works")
    else:
        print("   ✗ Module execution FAILED")
        print(f"   Error: {result.stderr}")
except Exception as e:
    print(f"   ✗ Module execution test failed: {e}")

print("\n" + "=" * 60)
print("DIAGNOSIS COMPLETE")
print("=" * 60)
print("\nIf you see any ✗ marks above, those are the issues to fix!")
print("Common solutions:")
print("1. Activate conda environment: conda activate twisstntern_env")
print("2. Navigate to project directory: cd /home/hlifchit/projects/twissting_baby") 
print("3. Reinstall package: pip install -e .")
print("4. If pandas import fails: conda install pandas") 