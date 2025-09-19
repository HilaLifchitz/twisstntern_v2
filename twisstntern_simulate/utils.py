#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import warnings
from math import sqrt, log

# from sympy import Eq, Symbol as sym, solve
from scipy.stats import binom, chi2
import sys

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


# Function ternary_coord: Convert Cartesian coordinates to ternary coordinates.
def ternary_coord(x, y):
    """
    Convert Cartesian coordinates (x, y) to ternary coordinates (t1, t2, t3).

    This function accepts either:
    - Scalars (floats), returning a tuple of floats (t1, t2, t3)
    - NumPy arrays, returning a tuple of arrays (t1, t2, t3)

    Assumes the triangle has vertices at:
        A = (0, h), B = (-0.5, 0), C = (+0.5, 0)
    """
    x_arr = np.asarray(x)
    y_arr = np.asarray(y)

    t1 = y_arr / h
    t3 = x_arr + 0.5 * (1 - t1)
    t2 = 1 - t1 - t3

    if np.isscalar(x) and np.isscalar(y):
        # Return plain floats if scalar input
        return float(t1), float(t2), float(t3)
    else:
        # Return NumPy arrays
        return t1, t2, t3


# Function cartizian: Convert ternary coordinates to Cartesian coordinates.
def cartizian(t1, t2, t3):
    """Convert ternary coordinates to Cartesian coordinates.

    Args:
        t1, t2, t3: Ternary coordinates

    Returns:
        tuple: (x, y) Cartesian coordinates
    """
    return ((t3 - t2) / 2, t1 * h)


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
    Return the x-limits [x_L, x_R] of the isocline where tau_2 = t2.
    """
    x_L = -0.5 * t2
    x_R = 0.5 - t2
    return (x_L, x_R)


def T3_lim(t3):
    """
    Return the x-limits [x_L, x_R] of the isocline where tau_3 = t3.
    """
    x_L = 0.5 * (2 * t3 - 1)
    x_R = 0.5 * t3
    return (x_L, x_R)


# plotting only the right half- for the results plot
def T3_lim_symm(t3):
    """
    Return the x-limits [x_L, x_R] of the isocline where tau_3 = t3,
    restricted to the right half of the triangle (i.e., x ≥ 0).
    """
    x_L = 0.5 * (2 * t3 - 1)  # from T3(t3, x) = 0
    x_R = 0.5 * t3  # from T3(t3, x) = T2(0, x)
    return (max(x_L, 0), x_R)


# Function: return_triangle_coord
#
# Given barycentric-like coordinates for the three vertices of a subtriangle
# (a1, b1), (a2, b2), and (a3, b3), this function computes:
#   1. The Cartesian x- and y-coordinates of the triangle corners.
#   2. The orientation ("up" or "down") of the triangle.
#
# The orientation is determined based on the vertical alignment (y-values)
# of the triangle's base relative to the apex (spitz) vertex.


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


###########################################


def dump_data(file, logger=None, axis_order=None):
    """Load and process data from CSV file.

    Args:
        file: Path to the input CSV file
        logger: Optional logger instance for logging messages
        axis_order: Optional list of column names to reorder axes (e.g., ["T2", "T1", "T3"])
                   Default is ["T1", "T2", "T3"]

    Returns:
        pandas.DataFrame: Processed data with normalized coordinates

    Raises:
        ValueError: If no valid numeric data is found
    """
    warnings.filterwarnings("ignore", category=FutureWarning)

    print(f"\nAttempting to read file: {file}")

    # First, read the first 10 non-comment lines to analyze
    with open(file, "r", encoding="utf-8") as f:
        lines = []
        for line in f:
            if not line.strip().startswith("#"):
                lines.append(line)
                if len(lines) >= 10:
                    break

    # Try each delimiter and find the one that gives us 3 numeric columns
    best_delimiter = None
    best_skiprows = 0
    best_data = None

    for delimiter in [",", "\t", ";", " "]:
        print(f"\nTrying delimiter: '{delimiter}'")

        # Try different numbers of rows to skip
        for skiprows in range(10):  # Try skipping 0-9 rows
            try:
                # Read just the first 10 rows to check
                test_data = pd.read_csv(
                    file, sep=delimiter, skiprows=skiprows, nrows=10
                )

                # Convert all columns to numeric, coercing errors to NaN
                for col in test_data.columns:
                    test_data[col] = pd.to_numeric(test_data[col], errors="coerce")

                # Drop completely empty columns
                test_data = test_data.dropna(axis=1, how="all")

                # Count how many columns have numeric data
                numeric_cols = []
                for col in test_data.columns:
                    if test_data[col].notna().any():
                        numeric_cols.append(col)

                print(
                    f"  Skip {skiprows} rows: Found {len(numeric_cols)} numeric columns"
                )

                # If we found exactly 3 numeric columns, this is our best match so far
                if len(numeric_cols) == 3:
                    best_delimiter = delimiter
                    best_skiprows = skiprows
                    best_data = test_data[numeric_cols]
                    print(
                        f"  Found 3 numeric columns with delimiter '{delimiter}' skipping {skiprows} rows"
                    )
                    break

            except Exception as e:
                continue

        if best_data is not None:
            break

    if best_data is None:
        raise ValueError(
            "Could not find 3 numeric columns in the first 10 rows with any common delimiter. "
            "Please check your data format."
        )

    print(f"\nUsing delimiter: '{best_delimiter}', skipping {best_skiprows} rows")

    # Now read the entire file with the best parameters
    try:
        data = pd.read_csv(file, sep=best_delimiter, skiprows=best_skiprows)
        print(f"Successfully read entire file")
    except Exception as e:
        raise ValueError(f"Failed to read file with best parameters: {str(e)}")

    # If we have more than 3 columns, try to identify the relevant columns
    if len(data.columns) > 3:
        # First try to find columns with names containing T1, T2, T3
        t_cols = []
        for col in data.columns:
            if any(f"T{i}" in str(col).upper() for i in [1, 2, 3]):
                t_cols.append(col)

        if len(t_cols) == 3:
            data = data[t_cols]
        else:
            # Use the same columns we found in our test
            data = data[best_data.columns]

    # Rename columns to standard format or use provided axis order
    if axis_order is None:
        axis_order = ["T1", "T2", "T3"]
    data.columns = axis_order

    # Log axis order if logger is provided
    if logger:
        logger.info(f"Using axis order: {axis_order}")

    # Convert to numeric and handle any non-numeric values
    for col in data.columns:
        data[col] = pd.to_numeric(data[col], errors="coerce")

    data = data.dropna()

    if len(data) == 0:
        raise ValueError("No valid numeric data found after processing")

    # Normalize the data
    n_rows = data.shape[0]
    for i in range(n_rows):
        s = sum(data.iloc[i, :])
        if s == 0:
            continue
        data.iloc[i, :] = (data.iloc[i, :]) / s

    data = data.loc[data.iloc[:, 1] != data.iloc[:, 2]]
    print(
        f"\nData loaded and normalized. Points with T2 == T3 removed. Remaining rows: {len(data)}"
    )

    return data


###################################################


def dump_data1(file):
    """Load and process data from CSV file.

    Args:
        file: Path to the input CSV file

    Returns:
        pandas.DataFrame: Processed data with normalized coordinates

    Raises:
        ValueError: If no valid numeric data is found
    """
    warnings.filterwarnings("ignore", category=FutureWarning)

    print(f"\nAttempting to read file: {file}")

    # First, try to detect the delimiter and number of header rows
    with open(file, "r", encoding="utf-8") as f:
        lines = []
        for line in f:
            if not line.strip().startswith("#"):
                lines.append(line)

        # Try to detect delimiter from first non-comment line
        if lines:
            first_line = lines[0]
            print(f"First non-comment line: {first_line.strip()}")
            # Check common delimiters
            for delim in [",", "\t", ";", " "]:
                if delim in first_line:
                    delimiter = delim
                    break
            else:
                delimiter = ","  # default to comma if no delimiter found

            # Try to detect header rows by checking consecutive lines
            skiprows = 0
            for i, line in enumerate(lines[:10]):  # Check first 10 lines at most
                try:
                    # Try to convert first value to float
                    float(line.split(delimiter)[0])
                    skiprows = i  # Found numeric data at line i
                    break
                except (ValueError, IndexError):
                    continue

    print(f"\nDetected delimiter: '{delimiter}', skipping {skiprows} header rows")

    # Try to read with pandas
    try:
        data = pd.read_csv(file, sep=delimiter, skiprows=skiprows)
        print(f"Successfully read with detected delimiter")
        print(f"Columns found: {list(data.columns)}")
    except Exception as e:
        print(f"Failed to read with detected delimiter: {str(e)}")
        # Try alternative delimiters as fallback
        for alt_delim in [",", "\t", ";", " "]:
            if alt_delim != delimiter:
                try:
                    data = pd.read_csv(file, sep=alt_delim, skiprows=skiprows)
                    print(
                        f"Successfully read with alternative delimiter: '{alt_delim}'"
                    )
                    print(f"Columns found: {list(data.columns)}")
                    break
                except Exception:
                    continue
        else:
            raise ValueError(f"Could not read file with any common delimiter")

    # If we have more than 3 columns, try to identify the relevant columns
    if len(data.columns) > 3:
        print(
            f"\nFound {len(data.columns)} columns. Attempting to identify relevant columns..."
        )

        # First try to find columns with names containing T1, T2, T3
        t_cols = []
        for col in data.columns:
            if any(f"T{i}" in str(col).upper() for i in [1, 2, 3]):
                t_cols.append(col)
                print(f"Found T-column: {col}")

        if len(t_cols) == 3:
            print("Using T1, T2, T3 columns")
            data = data[t_cols]
        else:
            print("\nTrying to find numeric columns...")
            # If we can't find named columns, try to find numeric columns
            # First convert all columns to numeric, coercing errors to NaN
            for col in data.columns:
                data[col] = pd.to_numeric(data[col], errors="coerce")

            # Drop any completely empty columns
            data = data.dropna(axis=1, how="all")
            print(f"Columns after dropping empty: {list(data.columns)}")

            # Get the first 3 columns that have at least one non-NaN value
            numeric_cols = []
            for col in data.columns:
                if data[col].notna().any():
                    numeric_cols.append(col)
                    print(f"Found numeric column: {col}")
                    if len(numeric_cols) == 3:
                        break

            if len(numeric_cols) == 3:
                print("Using first three numeric columns")
                data = data[numeric_cols]
            else:
                print("\nTrying to find numeric data in first 10 rows...")
                # If we still don't have 3 columns, try to find any numeric data
                # by examining the first 10 rows of each column
                numeric_cols = []
                for col in data.columns:
                    # Check if any of the first 10 values are numeric
                    if any(pd.to_numeric(data[col].iloc[:10], errors="coerce").notna()):
                        numeric_cols.append(col)
                        print(f"Found column with numeric data in first 10 rows: {col}")
                        if len(numeric_cols) == 3:
                            break

                if len(numeric_cols) == 3:
                    print("Using columns with numeric data in first 10 rows")
                    data = data[numeric_cols]
                else:
                    print("\nColumn detection failed. Available columns:")
                    for col in data.columns:
                        print(f"Column: {col}")
                        print(f"First few values: {data[col].iloc[:5].tolist()}")
                    raise ValueError(
                        "Could not identify three numeric columns for ternary coordinates. "
                        "Please check your data format."
                    )

    # Rename columns to standard format
    data.columns = ["T1", "T2", "T3"]

    # Convert to numeric and handle any non-numeric values
    for col in data.columns:
        data[col] = pd.to_numeric(data[col], errors="coerce")

    data = data.dropna()

    if len(data) == 0:
        raise ValueError("No valid numeric data found after processing")

    # Normalize the data
    n_rows = data.shape[0]
    for i in range(n_rows):
        s = sum(data.iloc[i, :])
        if s == 0:
            continue
        data.iloc[i, :] = (data.iloc[i, :]) / s

    data = data.loc[data.iloc[:, 1] != data.iloc[:, 2]]
    print(f"\nSuccessfully read {len(data)} rows of data")

    return data


#  Data point counting in ternary space in the subtriangle corresponding to the given coordinates a1,b1,a2,b2,a3,b3
#
# These functions compute how many points from a ternary-encoded dataset fall within:
#   - a given subtriangle defined by axis-aligned intervals on (T1, T2, T3), and
#   - its reflected counterpart (swapping T2 and T3).
#
# The data space is divided into half-open intervals (a, b], with a special case of [0, b]
# to ensure boundary points are included without double-counting.
#
# The function `n(...)` returns (n_r, n_l), the number of points in the right and left
# subtriangles, determined by the x-position of the triangle's apex in Cartesian space.
#


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


# Calculate the number of subtriangles in a ternary plot grid given a granularity parameter alpha.
#
# The triangle grid is defined by horizontal rows of triangle bases aligned with T1 coordinates,
# starting at T1 = 0 and increasing in steps of size alpha. Each row contains alternating
# up- and down-oriented triangles ("diamonds") formed by splitting each rhombus along its horizontal base.
#
# To ensure that the vertical y-axis (x=0) does not intersect the midpoint of any triangle,
# alpha is chosen such that 1/alpha is an even number.
#
# In the bottom row (T1 = 0), there are exactly (1 / (2 * alpha)) upward-facing triangles.
# Each subsequent row contains two fewer upward triangles. For each such triangle, its "diamond"
# (comprising one up and one down triangle) contributes four triangles in total.
#
# The function sums all triangles across rows until only one up-triangle remains, and returns
# the total number of subtriangles that must be considered.


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
        sys.exit()

    half_inv = int(inv_alpha // 2)  # number of up-triangles in bottom row
    total = half_inv  # start with bottom row

    for i in range(1, half_inv):
        total += 4 * (
            half_inv - i
        )  # each diamond row contributes 4 triangles per up-down pair

    return total


# Generate a list of coordinate tuples for all admissible right-hand-side subtriangles to be used in the analysis.
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


######################### Symmetry tests (move to analysis.py?) ################################


def D_LR(n_r, n_l):
    """
    Calculate the directional asymmetry statistic D_LR = (n_r - n_l) / (n_r + n_l),
    which captures the imbalance between right and left triangle counts.
    """
    total = n_r + n_l
    return np.nan if total == 0 else (n_r - n_l) / total


def log_likelihood_ratio_test(n_r, n_l):
    """
    Compute the likelihood ratio test statistic for asymmetry between right and left triangle counts.

    The test compares:
        - Null hypothesis: p = 0.5 (symmetric)
        - Alternative hypothesis: p = n_l / (n_r + n_l) (empirical)

    The test statistic is:
        -2 * log(likelihood_null / likelihood_alt)
    which follows a chi-squared distribution with 1 degree of freedom under the null.

    Args:
        n_r (int): Count in the right triangle
        n_l (int): Count in the left triangle

    Returns:
        tuple[float, float]: (G-statistic, p-value). Returns (np.nan, np.nan) if no data.
    """
    N = n_r + n_l

    if N == 0:
        return np.nan, np.nan

    p_null = 0.5
    p_alt = n_l / N

    L_null = binom.pmf(n_l, N, p_null)
    L_alt = binom.pmf(n_l, N, p_alt)

    if (
        L_alt == 0 or L_null == 0 or L_null / L_alt < 1e-16
    ):  # this means sth extremely unlikely has occured
        return np.nan, 0.0

    test_stat = -2 * log(L_null / L_alt)
    # degress of freedom = 1, according to https://stats.libretexts.org/Bookshelves/Applied_Statistics/Biological_Statistics_(McDonald)/02%3A_Tests_for_Nominal_Variables/2.04%3A_GTest_of_Goodness-of-Fit
    p_value = 1 - chi2.cdf(test_stat, df=1)
    return float(test_stat), float(p_value)


def download_config_template(destination_path=None):
    """
    Download the config_template.yaml file from GitHub.

    This function downloads the latest config template directly from the
    GitHub repository, ensuring users always get the most up-to-date version.

    Args:
        destination_path (str, optional): Where to save the config file.
                                        Defaults to current directory.

    Returns:
        str: Path where the config file was saved

    Examples:
        >>> import twisstntern_simulate.utils as utils
        >>> config_path = utils.download_config_template()
        >>> print(f"Config template downloaded to: {config_path}")
    """
    import requests
    import os

    if destination_path is None:
        destination_path = os.path.join(os.getcwd(), "config_template.yaml")

    # GitHub raw URL for the config template
    url = "https://raw.githubusercontent.com/HilaLifchitz/twisstntern_v2/main/config_template.yaml"

    try:
        print("📥 Downloading config template from GitHub...")
        response = requests.get(url, timeout=10)
        response.raise_for_status()

        with open(destination_path, "w", encoding="utf-8") as f:
            f.write(response.text)

        print(f"✅ Config template downloaded to: {destination_path}")
        return destination_path

    except requests.exceptions.RequestException as e:
        print(f"❌ Failed to download config template: {e}")
        print("🔗 Please download manually from:")
        print(
            "   https://github.com/HilaLifchitz/twisstntern_v2/raw/main/config_template.yaml"
        )
        return None
    except Exception as e:
        print(f"❌ Error saving config template: {e}")
        return None
