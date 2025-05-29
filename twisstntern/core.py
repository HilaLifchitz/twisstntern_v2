#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import warnings
from math import sqrt, log
from sympy import Eq, Symbol as sym, solve
from scipy.stats.distributions import chi2

# Constants
h = sqrt(3)/2  # Height of the equilateral triangle
a = 2*h  # slope multiplier used in T2 and T3


# Define the isocline equations for the ternary simplex.
# Each function T_i(tau_i, x) returns the Cartesian y-coordinate of the line
# along which the barycentric coordinate tau_i is constant.
#
# In a ternary plot (i.e., an equilateral triangle where each point represents a triplet of weights
# (tau_1, tau_2, tau_3) summing to 1), every point inside the triangle corresponds to a unique
# convex combination of the three vertices.
#
# These functions describe, in Cartesian coordinates, the geometric level sets of the barycentric
# coordinates. That is: given a fixed value tau_i in [0, 1], T_i(tau_i, x) gives the y-position
# of the point(s) in the triangle where the i-th barycentric coordinate equals tau_i.
#
# These lines — also called isoclines — are used to visualize or compute slices of the simplex
# where one barycentric coordinate is held constant. For instance:
#   - T1 traces horizontal slices (tau_1 constant),
#   - T2 traces left-leaning slices (tau_2 constant),
#   - T3 traces right-leaning slices (tau_3 constant).
#
# This is useful for plotting
def T1(y, x):
    """Calculate T1 coordinate in ternary space.
    
    Args:
        y, x: Input coordinates
        
    Returns:
        float: T1 coordinate (y*h)
    """
    return y*h

def T2(y, x):
    """Calculate T2 coordinate in ternary space.
    
    Args:
        y, x: Input coordinates
        
    Returns:
        float: T2 coordinate (-2h*x + sqrt(3)*(0.5-y))
    """
    b = sqrt(3)*(0.5-y) 
    return -a*x + b

def T3(y, x):
    """Calculate T3 coordinate in ternary space.
    
    Args:
        y, x: Input coordinates
        
    Returns:
        float: T3 coordinate (2h*x + sqrt(3)*(0.5-y))
    """
    b = sqrt(3)*(0.5-y) 
    return a*x + b

def ternary_coord(a, b):
    """Convert Cartesian coordinates to ternary coordinates.
    
    Args:
        a, b: Cartesian coordinates
        
    Returns:
        tuple: (t1, t2, t3) ternary coordinates
    """
    t1 = b/h
    
    y = sym('y')
    eqa_2 = Eq(T2(y,a), b)
    t2 = solve(eqa_2,y)[0]
    
    z = sym('z')
    eqa_3 = Eq(T3(z,a), b)
    t3 = solve(eqa_3,z)[0]
    
    return (t1,t2,t3)

def cartizian(x, y, z):
    """Convert ternary coordinates to Cartesian coordinates.
    
    Args:
        x, y, z: Ternary coordinates
        
    Returns:
        tuple: (x, y) Cartesian coordinates
    """
    return ((z-y)/2, x*h)

def T1_lim(y):
    """Calculate T1 limit in ternary space.
    
    Args:
        y: y-coordinate
        
    Returns:
        tuple: (x_L, x_R) x-axis limits
    """
    y = y*h  # the real height of coordinate y in T1
    x_l = sym('x_l')  # x_left
    eqa_l = Eq(T3(0,x_l),y)
    x_L = float(solve(eqa_l,x_l)[0])

    x_r = sym('x_r')  # x_right
    eqa_r = Eq(T2(0,x_r),y)
    x_R = float(solve(eqa_r,x_r)[0])
    return (x_L,x_R)

def T2_lim(y):
    """Calculate T2 limit in ternary space.
    
    Args:
        y: y-coordinate
        
    Returns:
        tuple: (x_L, x_R) x-axis limits
    """
    x_l = sym('x_l')  # x_left
    eqa_l = Eq(T3(0,x_l),T2(y, x_l))
    x_L = float(solve(eqa_l,x_l)[0])

    x_r = sym('x_r')  # x_right
    eqa_r = Eq(T2(y,x_r),0)
    x_R = float(solve(eqa_r,x_r)[0])
    return (x_L,x_R)

def T3_lim(y):
    """Calculate T3 limit in ternary space.
    
    Args:
        y: y-coordinate
        
    Returns:
        tuple: (x_L, x_R) x-axis limits
    """
    x_l = sym('x_l')  # x_left
    eqa_l = Eq(T3(y,x_l),0)
    x_L = float(solve(eqa_l,x_l)[0])

    x_r = sym('x_r')  # x_right
    eqa_r = Eq(T3(y,x_r),T2(0, x_r))
    x_R = float(solve(eqa_r,x_r)[0])
    return (x_L,x_R)

def T3_lim_symm(y):
    """Calculate symmetric T3 limit in ternary space.
    
    Args:
        y: y-coordinate
        
    Returns:
        tuple: (x_L, x_R) x-axis limits for right half
    """
    x_l = sym('x_l')  # x_left
    eqa_l = Eq(T3(y,x_l),0)
    x_L = float(solve(eqa_l,x_l)[0])

    x_r = sym('x_r')  # x_right
    eqa_r = Eq(T3(y,x_r),T2(0, x_r))
    x_R = float(solve(eqa_r,x_r)[0])
    
    return (max(x_L,0),x_R)

def return_triangle_coord(a1, b1, a2, b2, a3, b3):
    """Calculate triangle coordinates and direction.
    
    Args:
        a1, b1: First vertex coordinates
        a2, b2: Second vertex coordinates
        a3, b3: Third vertex coordinates
        
    Returns:
        tuple: (trianglex, triangley, direction) where:
            - trianglex: List of x coordinates
            - triangley: List of y coordinates
            - direction: "up" or "down"
    """
    direction = None
    triangley = None
    
    x_s = sym('x_s')  # s- spitz node
    x_l = sym('x_l')  # l- left node
    x_r = sym('x_r')  # r- right node

    eqa_r = Eq(T2(a2, x_r), T3(b3,x_r))  # true for both up & down triangles
    x_r = float(solve(eqa_r,x_r)[0])

    eqa_l = Eq(T2(b2, x_l), T3(a3,x_l))  # true for both up & down triangles
    x_l = float(solve(eqa_l,x_l)[0])
    
    # A triangle is classified as an "up-triangle" if and only if the y-component
    # of any of the base nodes intersects with b1/h
    if abs(T2(b2,x_l)-a1*h) <= 10**-14:  # T2(b2,x_l) == (a1*h)
        direction = "up"
        
        eqa_s = Eq(T2(a2, x_s), T3(a3,x_s))
        x_s = float(solve(eqa_s,x_s)[0])

        triangley = [a1*h, b1*h, a1*h, a1*h]  # y coordinates of [x_l,x_s,x_r,x_l]
    
    elif abs(T2(b2,x_l)-b1*h) <= 10**-14:  # T2(b2,x_l) == (b1*h)
        direction = "down"
    
        eqa_s = Eq(T2(b2, x_s), T3(b3,x_s))
        x_s = float(solve(eqa_s,x_s)[0])
    
        triangley = [b1*h, a1*h, b1*h, b1*h]  # y coordinates of [x_l,x_s,x_r,x_l]
    
    else:
        raise ValueError("Couldn't classify triangle as 'up' or 'down'. Check coordinates.")

    trianglex = [x_l, x_s, x_r, x_l]
    return (trianglex, triangley, direction)

def dump_data(file):
    """Load and process data from CSV file.
    
    Args:
        file: Path to the input CSV file
        
    Returns:
        pandas.DataFrame: Processed data with normalized coordinates
        
    Raises:
        ValueError: If no valid numeric data is found
    """
    warnings.filterwarnings('ignore', category=FutureWarning)
    
    print(f"\nAttempting to read file: {file}")
    
    # Read the file and skip comment lines
    with open(file, 'r', encoding='utf-8') as f:
        lines = []
        for line in f:
            if not line.strip().startswith('#'):
                lines.append(line)
    
    print("\nFirst few lines of the file (after removing comments):")
    for i, line in enumerate(lines[:5]):
        print(f"Line {i}: {line.strip()}")
    
    # Try to read with pandas, skipping the header row
    try:
        data = pd.read_csv(file, sep='\t', skiprows=3, names=['T1', 'T2', 'T3'])
        print(f"Successfully read with tab delimiter, skipping comments and header")
    except Exception as e:
        print(f"Failed to read with tab delimiter: {str(e)}")
        raise
    
    print(f"\nSuccessfully read {len(data)} rows of data")
    
    # Convert to numeric and handle any non-numeric values
    for col in data.columns:
        data[col] = pd.to_numeric(data[col], errors='coerce')
    
    data = data.dropna()
    
    if len(data) == 0:
        raise ValueError("No valid numeric data found after processing")

    # Normalize the data
    n_rows = data.shape[0]
    for i in range(n_rows):
        s = sum(data.iloc[i,:])
        if s == 0:
            continue
        data.iloc[i,:] = (data.iloc[i,:])/s
    
    data = data.loc[data.iloc[:,1] != data.iloc[:,2]]
   
    return data 

def n(a1, b1, a2, b2, a3, b3, data):
    """Count points in triangle"""
    n_points = 0
    for i in range(data.shape[0]):
        x = cartizian(data.iloc[i,0], data.iloc[i,1], data.iloc[i,2])[0]
        y = cartizian(data.iloc[i,0], data.iloc[i,1], data.iloc[i,2])[1]
        
        trianglex, triangley, direction = return_triangle_coord(a1, b1, a2, b2, a3, b3)
        
        if direction == "up":
            if (y >= triangley[0] and y <= triangley[1] and 
                x >= trianglex[0] and x <= trianglex[2]):
                n_points += 1
        else:
            if (y >= triangley[0] and y <= triangley[1] and 
                x >= trianglex[0] and x <= trianglex[2]):
                n_points += 1
                
    return n_points

def D_LR(n_r, n_l):
    """Calculate D_LR statistic"""
    if n_r + n_l == 0:
        return np.nan
    return (n_r - n_l)/(n_r + n_l)

def log_likelihood_ratio_test(n_r, n_l):
    """Perform log-likelihood ratio test"""
    if n_r + n_l == 0:
        return np.nan, np.nan
    
    n_total = n_r + n_l
    p = 0.5
    
    expected_r = n_total * p
    expected_l = n_total * p
    
    if expected_r == 0 or expected_l == 0:
        return np.nan, np.nan
    
    # Handle zero counts in observed values
    if n_r == 0:
        n_r_term = 0
    else:
        n_r_term = n_r * log(n_r/expected_r)
        
    if n_l == 0:
        n_l_term = 0
    else:
        n_l_term = n_l * log(n_l/expected_l)
    
    g_test = 2 * (n_r_term + n_l_term)
    p_value = 1 - chi2.cdf(g_test, 1)
    
    return g_test, p_value

def fundemental_asymmetry(data):
    """Calculate fundamental asymmetry statistics"""
    n_r = n(0, 0.5, 0, 0.5, 0, 0.5, data)
    n_l = n(-0.5, 0, 0, 0.5, 0, 0.5, data)
    d_lr = D_LR(n_r, n_l)
    g_test, p_value = log_likelihood_ratio_test(n_r, n_l)
    return n_r, n_l, d_lr, g_test, p_value 