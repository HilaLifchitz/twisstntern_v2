#!/usr/bin/env python3
"""
Simple import test to isolate the pandas/simulation import issue
"""

print("Testing imports step by step...")

try:
    print("1. Testing pandas...")
    import pandas as pd
    print("   ✅ pandas imported successfully")
except Exception as e:
    print(f"   ❌ pandas failed: {e}")
    exit(1)

try:
    print("2. Testing twisstntern_simulate base...")
    import twisstntern_simulate
    print("   ✅ twisstntern_simulate imported successfully")
except Exception as e:
    print(f"   ❌ twisstntern_simulate failed: {e}")
    exit(1)

try:
    print("3. Testing simulation module specifically...")
    from twisstntern_simulate.simulation import run_simulation
    print("   ✅ simulation module imported successfully")
except Exception as e:
    print(f"   ❌ simulation module failed: {e}")
    print(f"   This is the problem! Error details: {e}")
    exit(1)

try:
    print("4. Testing module execution...")
    import subprocess
    import sys
    result = subprocess.run([sys.executable, "-m", "twisstntern_simulate", "--help"], 
                          capture_output=True, text=True, timeout=5)
    if result.returncode == 0:
        print("   ✅ Module execution works!")
        print("\n🎉 EVERYTHING WORKS! You should be able to run:")
        print("   python -m twisstntern_simulate -c config_template.yaml")
    else:
        print("   ❌ Module execution failed")
        print(f"   Error: {result.stderr}")
except Exception as e:
    print(f"   ❌ Module execution test failed: {e}")

print("\nDone!") 