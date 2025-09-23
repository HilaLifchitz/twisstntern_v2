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
    # Match legacy mapping so points stay within the [-0.5, 0.5] triangle bounds.
    x = (T3 - T2) / 2
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


def dump_data(file, logger=None, axis_order=None, normalize=True, remove_equal_t2_t3=True):
    """Load and process data from CSV file, matching legacy heuristics."""
    warnings.filterwarnings("ignore", category=FutureWarning)

    file_path = Path(file)
    if not file_path.exists():
        msg = f"File not found: {file_path}"
        if logger:
            logger.error(msg)
        raise FileNotFoundError(msg)

    print(f"\nAttempting to read file: {file_path}")
    if logger:
        logger.info(f"Attempting to read file: {file_path}")

    best_delimiter = None
    best_skiprows = 0
    best_data = None

    # Legacy heuristic: probe several delimiters and skipped header rows
    for delimiter in [",", "\t", ";", " "]:
        print(f"\nTrying delimiter: '{delimiter}'")
        if logger:
            logger.debug(f"Trying delimiter: '{delimiter}'")

        for skiprows in range(10):
            try:
                test_data = pd.read_csv(
                    file_path, sep=delimiter, skiprows=skiprows, nrows=10
                )
            except Exception:
                continue

            # Coerce to numeric and drop empty columns
            for col in test_data.columns:
                test_data[col] = pd.to_numeric(test_data[col], errors="coerce")
            test_data = test_data.dropna(axis=1, how="all")

            numeric_cols = [col for col in test_data.columns if test_data[col].notna().any()]
            print(f"  Skip {skiprows} rows: Found {len(numeric_cols)} numeric columns")
            if logger:
                logger.debug(
                    "Skip %d rows: found %d numeric columns", skiprows, len(numeric_cols)
                )

            if len(numeric_cols) == 3:
                best_delimiter = delimiter
                best_skiprows = skiprows
                best_data = test_data[numeric_cols]
                print(
                    f"  Found 3 numeric columns with delimiter '{delimiter}' skipping {skiprows} rows"
                )
                if logger:
                    logger.info(
                        "Using delimiter '%s' skipping %d rows", delimiter, skiprows
                    )
                break

        if best_data is not None:
            break

    if best_data is None:
        raise ValueError(
            "Could not find 3 numeric columns in the first 10 rows with any common delimiter. "
            "Please check your data format."
        )

    print(f"\nUsing delimiter: '{best_delimiter}', skipping {best_skiprows} rows")

    try:
        data = pd.read_csv(file_path, sep=best_delimiter, skiprows=best_skiprows)
        if logger:
            logger.info("Successfully read entire file with detected parameters")
    except Exception as exc:
        raise ValueError(f"Failed to read file with best parameters: {exc}")

    # Reduce to three numeric columns, matching legacy behaviour
    if len(data.columns) > 3:
        t_cols = [col for col in data.columns if any(f"T{i}" in str(col).upper() for i in [1, 2, 3])]
        if len(t_cols) == 3:
            data = data[t_cols]
        else:
            data = data[best_data.columns]

    if axis_order is None:
        axis_order = ["T1", "T2", "T3"]

    data.columns = axis_order

    for col in data.columns:
        data[col] = pd.to_numeric(data[col], errors="coerce")

    data = data.dropna()
    if len(data) == 0:
        raise ValueError("No valid numeric data found after processing")

    if normalize:
        n_rows = data.shape[0]
        for i in range(n_rows):
            row_sum = data.iloc[i, :].sum()
            if row_sum == 0:
                continue
            data.iloc[i, :] = data.iloc[i, :] / row_sum

    if remove_equal_t2_t3 and data.shape[1] >= 3:
        before = len(data)
        data = data.loc[data.iloc[:, 1] != data.iloc[:, 2]]
        removed = before - len(data)
    else:
        removed = 0

    status_parts = ["Data loaded"]
    if normalize:
        status_parts.append("normalized")
    if remove_equal_t2_t3:
        status_parts.append(f"removed {removed} points where T2 == T3")

    status_message = ". ".join(status_parts) + f". Remaining rows: {len(data)}"
    print(f"\n{status_message}")
    if logger:
        logger.info(
            "%s | normalize=%s remove_equal_t2_t3=%s remaining=%d",
            status_message,
            normalize,
            remove_equal_t2_t3,
            len(data),
        )

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


def log_likelihood_ratio_test(n_r, n_l):
    """Compute likelihood ratio (G) statistic and p-value for asymmetry."""
    N = n_r + n_l

    if N == 0:
        return np.nan, np.nan

    p_null = 0.5
    p_alt = n_l / N

    L_null = binom.pmf(n_l, N, p_null)
    L_alt = binom.pmf(n_l, N, p_alt)

    if L_alt == 0 or L_null == 0 or L_null / L_alt < 1e-16:
        return np.nan, 0.0

    test_stat = -2 * log(L_null / L_alt)
    p_value = 1 - chi2.cdf(test_stat, df=1)
    return float(test_stat), float(p_value)
