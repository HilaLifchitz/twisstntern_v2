#!/usr/bin/env python
# coding: utf-8

import os
from pathlib import Path
import numpy as np
import random
from math import *
import scipy
from scipy.stats import binom
from scipy.stats.distributions import chi2
import pandas as pd
from sympy import Eq, Symbol as sym, solve
import sys
import csv
import warnings

# Constants
h = sqrt(3)/2  # height of an equilateral triangle with side length 1

# Ternary coordinate functions
def T1(y, x):
    """Calculate T1 coordinate"""
    return y*h

def T2(y, x):
    """Calculate T2 coordinate"""
    a = 2*h
    b = sqrt(3)*(0.5-y) 
    return -a*x + b

def T3(y, x):
    """Calculate T3 coordinate"""
    a = 2*h
    b = sqrt(3)*(0.5-y) 
    return a*x + b

def ternary_coord(a, b):
    """Convert Cartesian coordinates to ternary coordinates"""
    t1 = b/h
    
    y = sym('y')
    eqa_2 = Eq(T2(y,a), b)
    t2 = solve(eqa_2,y)[0]
    
    z = sym('z')
    eqa_3 = Eq(T3(z,a), b)
    t3 = solve(eqa_3,z)[0]
    
    return (t1,t2,t3)

def cartizian(x, y, z):
    """Convert ternary coordinates to Cartesian coordinates"""
    return ((z-y)/2, x*h)

# Triangle limit functions
def T1_lim(y):
    """Return x-axis limits for T1(y,.)"""
    y = y*h
    x_l = sym('x_l')
    eqa_l = Eq(T3(0,x_l),y)
    x_L = float(solve(eqa_l,x_l)[0])

    x_r = sym('x_r')
    eqa_r = Eq(T2(0,x_r),y)
    x_R = float(solve(eqa_r,x_r)[0])
    return (x_L,x_R)

def T2_lim(y):
    """Return x-axis limits for T2(y,.)"""
    x_l = sym('x_l')
    eqa_l = Eq(T3(0,x_l),T2(y, x_l))
    x_L = float(solve(eqa_l,x_l)[0])

    x_r = sym('x_r')
    eqa_r = Eq(T2(y,x_r),0)
    x_R = float(solve(eqa_r,x_r)[0])
    return (x_L,x_R)

def T3_lim(y):
    """Return x-axis limits for T3(y,.)"""
    x_l = sym('x_l')
    eqa_l = Eq(T3(y,x_l),0)
    x_L = float(solve(eqa_l,x_l)[0])

    x_r = sym('x_r')
    eqa_r = Eq(T3(y,x_r),T2(0, x_r))
    x_R = float(solve(eqa_r,x_r)[0])
    return (x_L,x_R)

def T3_lim_symm(y):
    """Return x-axis limits for T3(y,.) in right half of triangle"""
    x_l = sym('x_l')
    eqa_l = Eq(T3(y,x_l),0)
    x_L = float(solve(eqa_l,x_l)[0])

    x_r = sym('x_r')
    eqa_r = Eq(T3(y,x_r),T2(0, x_r))
    x_R = float(solve(eqa_r,x_r)[0])
    
    return (max(x_L,0),x_R)

def return_triangle_coord(a1, b1, a2, b2, a3, b3):
    """Return triangle coordinates and orientation"""
    x_s = sym('x_s')
    x_l = sym('x_l')
    x_r = sym('x_r')

    eqa_r = Eq(T2(a2, x_r), T3(b3,x_r))
    x_r = float(solve(eqa_r,x_r)[0])

    eqa_l = Eq(T2(b2, x_l), T3(a3,x_l))
    x_l = float(solve(eqa_l,x_l)[0])
    
    if abs(T2(b2,x_l)-a1*h) <= 10**-14:
        direction = "up"
        eqa_s = Eq(T2(a2, x_s), T3(a3,x_s))
        x_s = float(solve(eqa_s,x_s)[0])
        triangley = [a1*h, b1*h, a1*h, a1*h]
    
    if abs(T2(b2,x_l)-b1*h) <= 10**-14:
        direction = "down"
        eqa_s = Eq(T2(b2, x_s), T3(b3,x_s))
        x_s = float(solve(eqa_s,x_s)[0])
        triangley = [b1*h, a1*h, b1*h, b1*h]
    
    trianglex = [x_l, x_s, x_r, x_l]
    return (trianglex, triangley, direction)

def ref(a1, b1, a2, b2, a3, b3):
    """Return reflected coordinates"""
    return (a1, b1, a3, b3, a2, b2)

def mid_point_triangle(a1, b1, a2, b2, a3, b3):
    """Calculate midpoint of triangle"""
    trianglex, triangley, direction = return_triangle_coord(a1, b1, a2, b2, a3, b3)
    mid_x = trianglex[1]
    
    if direction == "up":
        mid_y = triangley[0] + (trianglex[1]-trianglex[0])/2
    else:
        mid_y = triangley[0] + 2*(trianglex[0]-trianglex[1])/3
        
    return (mid_x, round(mid_y, 4))

def dump_data(file):
    """Load and process data from CSV file"""
    warnings.filterwarnings('ignore', category=FutureWarning)
    
    print(f"\nAttempting to read file: {file}")
    
    with open(file, 'r', encoding='utf-8') as f:
        sample = f.read(1024)
        f.seek(0)
        
        try:
            dialect = csv.Sniffer().sniff(sample)
            print(f"Detected format: delimiter='{dialect.delimiter}', quoting={dialect.quoting}")
        except:
            print("Could not detect format automatically, will try multiple formats")
            dialect = None
        
        lines = f.readlines()
    
    print("\nFirst few lines of the file:")
    for i, line in enumerate(lines[:5]):
        print(f"Line {i}: {line.strip()}")
    
    data = None
    delimiters = [',', '\t', ';', ' ', '|']
    
    if dialect:
        try:
            data = pd.read_csv(file, dialect=dialect, header=None, names=['T1', 'T2', 'T3'])
            print(f"Successfully read with detected dialect: {dialect.delimiter}")
        except Exception as e:
            print(f"Failed to read with detected dialect: {str(e)}")
    
    if data is None or data.shape[1] != 3:
        for sep in delimiters:
            try:
                data = pd.read_csv(file, sep=sep, header=None, names=['T1', 'T2', 'T3'])
                if data.shape[1] == 3:
                    print(f"Successfully read with delimiter: {sep}")
                    break
                else:
                    print(f"Found {data.shape[1]} columns with delimiter {sep}, trying next...")
            except Exception as e:
                print(f"Failed to read with delimiter {sep}: {str(e)}")
                continue
    
    if data is None or data.shape[1] != 3:
        print("\nTrying to read line by line...")
        rows = []
        for line in lines:
            line = line.strip()
            if not line:
                continue
                
            for sep in delimiters:
                cols = [col.strip() for col in line.split(sep)]
                if len(cols) == 3:
                    try:
                        row = [float(col) for col in cols]
                        rows.append(row)
                        break
                    except ValueError:
                        continue
        
        if rows:
            data = pd.DataFrame(rows, columns=['T1', 'T2', 'T3'])
            print("Successfully read data line by line")
    
    if data is None or data.shape[1] != 3:
        raise ValueError(f"Could not read file with 3 columns. Please check the file format.")
    
    print(f"\nSuccessfully read {len(data)} rows of data")
    
    for col in data.columns:
        data[col] = pd.to_numeric(data[col], errors='coerce')
    
    data = data.dropna()
    
    if len(data) == 0:
        raise ValueError("No valid numeric data found after processing")

    n_rows = data.shape[0]
    for i in range(n_rows):
        s = sum(data.iloc[i,:])
        if s == 0:
            continue
        data.iloc[i,:] = (data.iloc[i,:])/s
    
    data = data.loc[data.iloc[:,1] != data.iloc[:,2]]
   
    return data 