#!/usr/bin/env python3
"""
This script is used to compare the processing of the generator object (locus mode) performance.
It is used to compare the new and old ts_processing algorithms.
"""
import time



# First, let's try importing just what we need manually

# import yaml

# from dataclasses import dataclass


# from typing import List, Optional, Dict, Any


# from twisstntern_simulate.config import Config

from twisstntern_simulate.ts_processing import ts_to_twisst_weights as new
from twisstntern_simulate.simulation import run_simulation

from twisstntern.tree_processing import ts_to_twisst_weights as old

what = "locus"

res=run_simulation("/home/hlifchit/twisstntern2.0/twisstntern_simulate/config_template.yaml")

print("=== Comparing ts_processing algorithms  ===")
print(f"implemating {what} mode")
print("=== Testing new ts_processing alg ===")
start = time.perf_counter()

s1=new(res[what])

end = time.perf_counter()
print(f"New took {end - start:.4f} seconds")
print("--------------------------------")
# res=run_simulation("/home/hlifchit/twisstntern2.0/twisstntern_simulate/config_template.yaml")
# print("=== Testing old ts_processing alg ===")
# start = time.perf_counter()

# t1=new(res[what])

# end = time.perf_counter()
# print(f"Old took {end - start:.4f} seconds")

# print("=== Test complete ===") 
