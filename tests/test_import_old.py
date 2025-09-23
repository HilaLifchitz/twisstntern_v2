#!/usr/bin/env python3
"""
Simple import test to isolate the pandas/simulation import issue
"""

print("Testing imports step by step...")

# Ensure repository root is on sys.path for legacy imports
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

try:
    print("1. Testing pandas...")
    import pandas as pd
    print("   ✅ pandas imported successfully")
except Exception as e:
    print(f"   ❌ pandas failed: {e}")
    exit(1)

try:
    print("2. Testing legacy.twisstntern_simulate base...")
    import legacy.twisstntern_simulate as legacy_sim
    print("   ✅ legacy.twisstntern_simulate imported successfully")
except Exception as e:
    print(f"   ❌ legacy.twisstntern_simulate failed: {e}")
    exit(1)

try:
    print("3. Testing simulation module specifically...")
    from legacy.twisstntern_simulate.simulation import run_simulation
    print("   ✅ simulation module imported successfully")
except Exception as e:
    print(f"   ❌ simulation module failed: {e}")
    print(f"   This is the problem! Error details: {e}")
    exit(1)

try:
    print("4. Testing module execution...")
    import subprocess
    result = subprocess.run([sys.executable, "-m", "legacy.twisstntern_simulate", "--help"], 
                          capture_output=True, text=True, timeout=5)
    if result.returncode == 0:
        print("   ✅ Module execution works!")
        print("\n🎉 EVERYTHING WORKS! You should be able to run:")
        print("   python -m legacy.twisstntern_simulate -c config_template.yaml")
    else:
        print("   ❌ Module execution failed")
        print(f"   Error: {result.stderr}")
except Exception as e:
    print(f"   ❌ Module execution test failed: {e}")

print("\nDone!") 
