#!/usr/bin/env python
# coding: utf-8

"""
Core utilities for ternary coordinate system mathematics and statistical analysis.

This module contains shared functions used by both twisstntern and twisstntern_simulate
packages for working with ternary (barycentric) coordinate systems and performing
statistical analysis on topology weight data.
"""

import numpy as np
import pandas as pd
import warnings
from math import sqrt, log
from scipy.stats import binom, chi2
from pathlib import Path

# Constants
h = sqrt(3) / 2  # Height of the equilateral triangle
a = 2 * h  # slope multiplier used in functions T2 and T3


# Functions T1,T2,T3: Define the isocline equations for the ternary simplex.
# Each function T_i(tau_i, x) returns the Cartesian y-coordinate of the line
# along which the barycentric coordinate tau_i is constant and a given cartezian x is given.
# (if a symbolic 'x' is given it exactly gives the line in the triangle for which T_i is constant at tau_i)
#
# These lines — also called isoclines — are used to visualize or compute slices of the simplex
# where one barycentric coordinate is held constant. For instance:
#   - T1 traces horizontal slices (tau_1 constant),
#   - T2 traces left-leaning slices (tau_2 constant),
#   - T3 traces right-leaning slices (tau_3 constant).


def T1(tau_1, x):
    return tau_1 * h


def T2(tau_2, x):
    b = sqrt(3) * (0.5 - tau_2)
    return -a * x + b


def T3(tau_3, x):
    b = sqrt(3) * (0.5 - tau_3)
    return a * x + b


# T1_lim, T2_lim, T3_lim, T3_lim_symm: Functions for computing x-axis limits of isoclines in a ternary triangle.
#
# Each function T*_lim(tau_i) computes the [x_left, x_right] interval
# over which the isocline corresponding to a constant barycentric coordinate tau_i
# lies within the triangle (i.e., is visible and valid for plotting).
#
# These delimiters are useful for rendering the level set lines of tau_1, tau_2, or tau_3
# in Cartesian coordinates when the triangle is embedded as:
#   A = (0, h), B = (-0.5, 0), C = (+0.5, 0)
#   where h = sqrt(3)/2 (height of equilateral triangle)


def T1_lim(t1):
    """
    Return the x-limits [x_L, x_R] of the horizontal isocline where T1 coordinate tau_1 = t1.

    This corresponds to the valid x-range of the line y = t1 * h,
    intersected with the triangle's edges (T2 and T3).
    """
    y = t1 * h
    x_L = (y - h) / (2 * h)
    x_R = (h - y) / (2 * h)
    return (x_L, x_R)


def T2_lim(t2):
    """
    Calculate the x-limits for isocline T2 = t2.
    The range is where T2(t2, x) ≥ 0 within the triangle boundaries.
    """
    if t2 < 0 or t2 > 1:
        return None, None
#   A = (0, h), B = (-0.5, 0), C = (+0.5, 0)
    
    x_L = -0.5 * t2
    x_R = 0.5 - t2
    
    return x_L, x_R


def T3_lim(t3):
    """
    Calculate the x-limits for isocline T3 = t3.
    The range is where T3(t3, x) ≥ 0 within the triangle boundaries.
    """
    if t3 < 0 or t3 > 1:
        return None, None
    
    x_L = 0.5 * (2 * t3 - 1)
    x_R = 0.5 * t3
    
    return x_L, x_R


# plotting only the right half- for the results plot
def T3_lim_symm(t3):
    """
    Return the x-limits [x_L, x_R] of the isocline where tau_3 = t3,
    restricted to the right half of the triangle (i.e., x ≥ 0).
    """
    x_L = 0.5 * (2 * t3 - 1)  # from T3(t3, x) = 0
    x_R = 0.5 * t3  # from T3(t3, x) = T2(0, x)
    return (max(x_L, 0), x_R)


# Function ternary_coord: Convert Cartesian coordinates to ternary coordinates.
def ternary_coord(x, y):
    """
    Convert Cartesian coordinates (x, y) to ternary coordinates (t1, t2, t3).
    
    Using the inverse transformation:
        Given ternary-to-Cartesian transformation:
            x = t2 - t3
            y = h * t1 = (sqrt(3)/2) * t1
        
        where t1 + t2 + t3 = 1 (normalization constraint).
        
        The inverse system is:
            t1 = y / h = y / (sqrt(3)/2) = 2*y / sqrt(3)
            x = t2 - t3   →   t2 - t3 = x
            t1 + t2 + t3 = 1  →  t2 + t3 = 1 - t1
        
        Solving for t2 and t3:
        A = (0, h), B = (-0.5, 0), C = (+0.5, 0)
    """
    x_arr = np.array(x)
    y_arr = np.array(y)
    
    t1 = 2 * y_arr / sqrt(3)
    t3 = x_arr + 0.5 * (1 - t1)
    t2 = 1 - t1 - t3
    
    return t1, t2, t3


def cartizian(T1, T2, T3):
    """
    Convert ternary (barycentric) coordinates to Cartesian coordinates.
    
    T1, T2, T3 are the barycentric coordinates (should sum to 1).
    
    The transformation uses the standard vertices of an equilateral triangle:
    - A = (0, h)      — top vertex (T1 = 1)
    - B = (-0.5, 0)   — bottom-left vertex (T2 = 1)  
    - C = (0.5, 0)    — bottom-right vertex (T3 = 1)
    
    Cartesian coordinates are computed as:
        x = T1*0 + T2*(-0.5) + T3*(0.5) = T3 - T2
        y = T1*h + T2*0 + T3*0 = T1 * h
    
    Returns:
        tuple: (x, y) Cartesian coordinates
    """
    x = T3 - T2
    y = T1 * h
    
    return x, y


def return_triangle_coord(a1, b1, a2, b2, a3, b3):
    """Calculate triangle coordinates and direction.

    Args:
        a1, b1: First vertex coordinates
        a2, b2: Second vertex coordinates
        a3, b3: Third vertex coordinates

    Returns:
        tuple: (trianglex, triangley, direction) where:
            - trianglex: List of x coordinates (in order: left, spitz, right, left)
            - triangley: List of y coordinates (in order: left, spitz, right, left)
            - direction: "up" or "down"
    """
    direction = None
    triangley = None

    x_l = (a3 - b2) / 2  # left node
    x_r = (b3 - a2) / 2  # right node
    x_s = (b3 - b2) / 2  # spitz node

    # A triangle is classified as an "up-triangle" if and only if the y-component
    # of any of the base nodes intersects with a1*h (if the base nodes are down- so the spitz can be )
    if abs(T2(b2, x_l) - a1 * h) <= 10**-14:  # T2(b2,x_l) == (a1*h)
        direction = "up"

        triangley = [
            a1 * h,
            b1 * h,
            a1 * h,
            a1 * h,
        ]  # y coordinates of [x_l,x_s,x_r,x_l]

    elif abs(T2(b2, x_l) - b1 * h) <= 10**-14:  # T2(b2,x_l) == (b1*h)
        direction = "down"

        triangley = [
            b1 * h,
            a1 * h,
            b1 * h,
            b1 * h,
        ]  # y coordinates of [x_l,x_s,x_r,x_l]

    else:
        raise ValueError(
            "Couldn't classify triangle as 'up' or 'down'. Check coordinates."
        )

    trianglex = [x_l, x_s, x_r, x_l]
    return (trianglex, triangley, direction)


# reflected coordinate of the reflected triangle
def ref(a1, b1, a2, b2, a3, b3):
    return (a1, b1, a3, b3, a2, b2)


def mid_point_triangle(a1, b1, a2, b2, a3, b3):
    """
    Compute the midpoint (in Cartesian coordinates) of a triangle defined by ternary intervals.
    This is used to place p-values or other statistics in the triangle when plotting.

    Args:
        a1, b1, a2, b2, a3, b3 (float): Ternary interval coordinates for T1, T2, T3.

    Returns:
        tuple: (x, y) midpoint for the triangle, with y rounded to 4 decimal places.
    """
    trianglex, triangley, direction = return_triangle_coord(a1, b1, a2, b2, a3, b3)
    apex_x = trianglex[1]
    base_x = trianglex[0]
    base_y = triangley[0]

    if direction == "up":
        mid_y = base_y + (apex_x - base_x) / 2
    else:  # "down"
        mid_y = base_y + 2 * (base_x - apex_x) / 3

    return apex_x, round(mid_y, 4)


def build_conditions(a1, b1, a2, b2, a3, b3, data):
    """
    Construct inclusion conditions for points within a ternary subtriangle.

    This helper function defines half-open intervals of the form (a, b], except when a == 0,
    in which case the interval becomes [0, b] to ensure inclusion of boundary points and
    avoid double-counting along shared edges.

    Args:
        a1, b1, a2, b2, a3, b3 (float): Interval endpoints for T1, T2, T3.
        data (pd.DataFrame): DataFrame with ternary coordinate columns: T1, T2, T3.

    Returns:
        tuple: Boolean Series for inclusion in T1, T2, T3 directions.
    """
    if a1 == 0:
        c1 = a1 <= data.T1
    else:
        c1 = a1 < data.T1
    if a2 == 0:
        c2 = a2 <= data.T2
    else:
        c2 = a2 < data.T2
    if a3 == 0:
        c3 = a3 <= data.T3
    else:
        c3 = a3 < data.T3
    return c1, c2, c3


def n(a1, b1, a2, b2, a3, b3, data):
    """
    Count the number of points inside a ternary subtriangle and its reflection.

    Given the interval-based definition of a triangle in ternary space using coordinates
    (a1, b1), (a2, b2), (a3, b3), this function counts how many data points fall within:
        - the triangle defined by those intervals, and
        - its reflected counterpart (where T2 and T3 are swapped).

    It returns the counts as (n_r, n_l), where the direction is determined by checking
    the x-coordinate of the triangle's top node: if the node lies to the right of the
    vertical axis, the original triangle is considered a "right" triangle.

    Interval logic is defined as half-open: (a, b], except at a = 0, where [0, b] is used
    to include the first bin edge.

    Args:
        a1, b1, a2, b2, a3, b3 (float): Interval endpoints for T1, T2, T3.
        data (pd.DataFrame): DataFrame with ternary coordinate columns: T1, T2, T3.

    Returns:
        tuple[int, int]: (n_r, n_l), the number of data points in the right and left triangles.
    """
    condition_a1, condition_a2, condition_a3 = build_conditions(
        a1, b1, a2, b2, a3, b3, data
    )
    n_1 = (
        condition_a1
        & (data.T1 <= b1)
        & condition_a2
        & (data.T2 <= b2)
        & condition_a3
        & (data.T3 <= b3)
    ).sum()

    a1, b1, a2, b2, a3, b3 = ref(a1, b1, a2, b2, a3, b3)
    condition_a1, condition_a2, condition_a3 = build_conditions(
        a1, b1, a2, b2, a3, b3, data
    )
    n_2 = (
        condition_a1
        & (data.T1 <= b1)
        & condition_a2
        & (data.T2 <= b2)
        & condition_a3
        & (data.T3 <= b3)
    ).sum()

    # checking which triangle is left in which is right
    # flipping the triangle to the original coordinates given
    a1, b1, a2, b2, a3, b3 = ref(a1, b1, a2, b2, a3, b3)
    # The check involves determining whether the x-axis component of the top node of a subtriangle is positive or negative
    # in Cartesian coordinates. The top node is positive along the x-axis if and only if it belongs to the right-hand side triangle.
    trianglex, triangley, direction = return_triangle_coord(a1, b1, a2, b2, a3, b3)
    top_node = trianglex[
        1
    ]  # this is the x_axis coordinate of the top node in the given triangle

    if top_node > 0:  # the coordinates we were given were of a right-side triangle
        n_r = n_1
        n_l = n_2
    else:  # the coordinates we were given were of a left-side triangle
        n_r = n_2
        n_l = n_1

    return int(n_r), int(n_l)


def dump_data(file_path, logger=None):
    """
    Load and validate ternary data from CSV file.
    
    This function reads a CSV file containing topology weights or ternary coordinates
    and performs basic validation to ensure the data is suitable for analysis.
    
    Args:
        file_path (str): Path to the CSV file to load
        logger: Optional logger to log messages (default: None)
        
    Returns:
        pandas.DataFrame: Loaded and validated data with columns T1, T2, T3
        
    Raises:
        FileNotFoundError: If the file doesn't exist
        ValueError: If the data format is invalid or missing required columns
        
    Expected CSV format:
        - Must contain columns named 'T1', 'T2', 'T3' (case-sensitive)
        - Each row represents one data point with ternary coordinates
        - Values should be numeric and ideally sum to ~1.0 (barycentric constraint)
        
    Example:
        T1,T2,T3
        0.333,0.333,0.334
        0.5,0.25,0.25
        ...
    """
    # Check if file exists
    file_path = Path(file_path)
    if not file_path.exists():
        error_msg = f"File not found: {file_path}"
        if logger:
            logger.error(error_msg)
        raise FileNotFoundError(error_msg)
    
    if logger:
        logger.info(f"Loading data from: {file_path}")
    
    # Try to detect delimiter
    try:
        with open(file_path, 'r') as f:
            first_line = f.readline().strip()
            
        # Check common delimiters
        if ',' in first_line:
            delimiter = ','
        elif '\t' in first_line:
            delimiter = '\t'
        elif ';' in first_line:
            delimiter = ';'
        elif '|' in first_line:
            delimiter = '|'
        else:
            if logger:
                logger.warning("Could not detect delimiter, assuming comma")
            delimiter = ","  # default to comma if no delimiter found
            
        if logger:
            logger.debug(f"Detected delimiter: '{delimiter}'")
            
    except Exception as e:
        if logger:
            logger.warning(f"Could not detect delimiter: {e}, assuming comma")
        delimiter = ","
    
    # Load the data
    try:
        data = pd.read_csv(file_path, delimiter=delimiter)
        if logger:
            logger.info(f"Successfully loaded CSV with shape: {data.shape}")
            
    except Exception as e:
        error_msg = f"Failed to load CSV file: {e}"
        if logger:
            logger.error(error_msg)
        raise ValueError(error_msg)
    
    # Validate required columns
    required_columns = ['T1', 'T2', 'T3']
    missing_columns = [col for col in required_columns if col not in data.columns]
    
    if missing_columns:
        error_msg = f"Missing required columns: {missing_columns}. Found columns: {list(data.columns)}"
        if logger:
            logger.error(error_msg)
        raise ValueError(error_msg)
    
    if logger:
        logger.info(f"Found required columns: {required_columns}")
    
    # Convert to numeric if needed
    for col in required_columns:
        if not pd.api.types.is_numeric_dtype(data[col]):
            try:
                data[col] = pd.to_numeric(data[col], errors='coerce')
                if logger:
                    logger.info(f"Converted column '{col}' to numeric")
            except Exception as e:
                error_msg = f"Could not convert column '{col}' to numeric: {e}"
                if logger:
                    logger.error(error_msg)
                raise ValueError(error_msg)
    
    # Check for missing values
    missing_count = data[required_columns].isnull().sum().sum()
    if missing_count > 0:
        if logger:
            logger.warning(f"Found {missing_count} missing values in ternary columns")
        data = data.dropna(subset=required_columns)
        if logger:
            logger.info(f"Dropped rows with missing values. New shape: {data.shape}")
    
    # Basic validation: check if coordinates sum approximately to 1
    coordinate_sums = data[required_columns].sum(axis=1)
    deviation_from_1 = np.abs(coordinate_sums - 1.0)
    max_deviation = deviation_from_1.max()
    
    if max_deviation > 0.01:  # Allow for small floating point errors
        if logger:
            logger.warning(f"Some coordinates don't sum to 1.0 (max deviation: {max_deviation:.6f})")
            logger.warning("This might indicate non-normalized barycentric coordinates")
    else:
        if logger:
            logger.info("All coordinates sum approximately to 1.0 ✓")
    
    if logger:
        logger.info(f"Data validation completed. Final shape: {data.shape}")
    
    return data


def D_LR(n_r, n_l):
    """
    Calculate the directional asymmetry statistic D_LR = (n_r - n_l) / (n_r + n_l),
    which captures the imbalance between right and left triangle counts.
    """
    total = n_r + n_l
    return np.nan if total == 0 else (n_r - n_l) / total


def number_triangles(alpha):
    """
    Calculate the number of triangle regions in a ternary plot grid at granularity alpha.

    Args:
        alpha (float): Granularity of the grid. Must divide 1 evenly into an even number of intervals.

    Returns:
        int: Total number of subtriangles to consider in the grid.
    """
    inv_alpha = 1 / alpha

    # Ensure 1/alpha is even
    if int(inv_alpha) % 2 != 0:
        print("Error: 1/alpha must be even to avoid alignment with the vertical axis.")
        import sys
        sys.exit()

    half_inv = int(inv_alpha // 2)  # number of up-triangles in bottom row
    total = half_inv  # start with bottom row

    for i in range(1, half_inv):
        total += 4 * (
            half_inv - i
        )  # each diamond row contributes 4 triangles per up-down pair

    return total


def right_triangle_coordinates_list(alpha):
    """
    Generate a list of coordinate tuples for all admissible right-hand-side subtriangles
    in a ternary triangle grid (divided along T1 axis) with granularity alpha.

    A triangle is included only if its leftmost x-coordinate is non-negative
    (i.e., fully lies on the right side of the ternary diagram).

    Parameters:
        alpha (float): Granularity step (must evenly divide 1).

    Returns:
        list of list of tuples: Each triangle is represented by 3 (a,b) coordinate intervals
        corresponding to the T1, T2, and T3 axes.
    """
    coords_list = []

    # Iterate over horizontal rows along the T1 axis
    for row_index in range(int(1 / alpha)):
        a1 = row_index * alpha
        b1 = (row_index + 1) * alpha

        # Indices for T2 and T3 intervals, initially set based on the current row
        k_T2 = 0
        k_T3 = row_index

        while True:
            # --- UP TRIANGLE ---
            a2_up = k_T2 * alpha
            b2_up = (k_T2 + 1) * alpha

            a3_up = 1 - (k_T3 + 1) * alpha
            b3_up = 1 - k_T3 * alpha

            # Get Cartesian coordinates of triangle corners
            triangle_x, _, _ = return_triangle_coord(a1, b1, a2_up, b2_up, a3_up, b3_up)

            # Stop if triangle extends past the vertical midline (i.e. x < 0)
            if round(triangle_x[0], 4) < 0:
                break

            # Store the triangle's T1-T2-T3 coordinates
            coords_list.append(
                [
                    (round(a1, 4), round(b1, 4)),
                    (round(a2_up, 4), round(b2_up, 4)),
                    (round(a3_up, 4), round(b3_up, 4)),
                ]
            )

            # --- DOWN TRIANGLE ---
            k_T3 += 1  # Shift the base right to the next diamond-pair

            a2_down = k_T2 * alpha
            b2_down = (k_T2 + 1) * alpha

            a3_down = 1 - (k_T3 + 1) * alpha
            b3_down = 1 - k_T3 * alpha

            triangle_x, _, _ = return_triangle_coord(
                a1, b1, a2_down, b2_down, a3_down, b3_down
            )

            # Stop again if this triangle crosses over the y-axis
            if round(triangle_x[0], 4) < 0:
                break

            coords_list.append(
                [
                    (round(a1, 4), round(b1, 4)),
                    (round(a2_down, 4), round(b2_down, 4)),
                    (round(a3_down, 4), round(b3_down, 4)),
                ]
            )

            k_T2 += 1  # Prepare next pair of up/down triangles along the current row

    coords_list = pd.DataFrame(coords_list)
    # Assign descending index for plotting (bottom-up row indexing)
    coords_list["index"] = list(range(number_triangles(alpha), 0, -1))

    return coords_list


def log_likelihood_ratio_test(n_left, n_right):
    """
    Perform G-test (likelihood ratio test) for symmetry.
    
    Tests the null hypothesis that the probability of falling on either side
    is equal (p = 0.5) against the alternative that they differ.
    
    Args:
        n_left (int): Number of observations on the left side
        n_right (int): Number of observations on the right side
        
    Returns:
        tuple: (G_statistic, p_value) where:
               - G_statistic: The G-test statistic (follows chi-square distribution)
               - p_value: Two-tailed p-value from chi-square distribution
               
    The G-test is more appropriate than chi-square test for count data and
    has better properties for small sample sizes.
    
    Mathematical formulation:
        - Null hypothesis: p = 0.5 (symmetric)
        - G = 2 * Σ(observed * ln(observed/expected))
        - Under H0: G ~ χ²(df=1)
    """
    total = n_left + n_right
    
    if total == 0:
        return np.nan, 1.0
        
    # Expected counts under null hypothesis (p = 0.5)
    expected = total / 2.0
    
    # Handle zero counts (add small constant to avoid log(0))
    n_left_adj = max(n_left, 1e-10)
    n_right_adj = max(n_right, 1e-10)
    
    # G-statistic calculation
    G = 2 * (n_left_adj * log(n_left_adj / expected) + 
             n_right_adj * log(n_right_adj / expected))
    
    # Calculate p-value using chi-square distribution with df=1
    p_value = 1.0 - chi2.cdf(G, df=1)
    
    # For two-tailed test
    if p_value > 0.5:
        return np.nan, 0.0
        
    return G, p_value