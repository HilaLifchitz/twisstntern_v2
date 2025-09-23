#!/usr/bin/env python
# coding: utf-8

import pandas as pd
from .utils import (
    cartizian,
    return_triangle_coord,
    n,
    D_LR,
    log_likelihood_ratio_test,
    number_triangles,
)


def fundamental_asymmetry(data):
    """
    Perform basic symmetry analysis between the two main subtriangles (left vs. right of the y-axis).

    Assumes:
        - Data has columns 'T1', 'T2', 'T3' (ternary coordinates).
        - Points on the y-axis (x == 0) have been removed during preprocessing.

    Returns:
        tuple: (
            main_n_r       - number of points in the right subtriangle (x > 0),
            main_n_l       - number of points in the left subtriangle (x < 0),
            main_d_lr      - directional asymmetry statistic D_LR,
            main_g_test    - G-statistic from likelihood ratio test,
            main_p_value   - corresponding p-value
        )
    """
    # Convert ternary coordinates to Cartesian x-coordinates
    data["x-axis"] = cartizian(data["T1"], data["T2"], data["T3"])[0]

    # Boolean masks for each side
    right_mask = data["x-axis"] > 0
    left_mask = data["x-axis"] < 0

    # Count points on each side
    main_n_r = right_mask.sum()
    main_n_l = left_mask.sum()

    # Compute asymmetry statistics
    main_d_lr = D_LR(main_n_r, main_n_l)
    main_g_test, main_p_value = log_likelihood_ratio_test(main_n_r, main_n_l)
    """
    Reminder: Under the Binomial model, the test compares:
    - Null hypothesis: p = 0.5 (symmetric)
    - Alternative hypothesis: p = n_l / (n_r + n_l) (empirical)
    --->  likelihood_null = binom.pmf(n_l, N, p_null)
          likelihood_alt = binom.pmf(n_l, N, p_alt)
    --->  log_likelihood_ratio_test(n_r, n_l) = -2 * log(likelihood_null / likelihood_alt)
    """
    return (
        int(main_n_r),
        int(main_n_l),
        float(main_d_lr),
        float(main_g_test),
        float(main_p_value),
    )


def triangles_analysis(data, granularity):
    """
    Analyze symmetry across subtriangles of a ternary diagram by computing, for each triangle:
      - the number of data points in left and right reflected triangle pairs,
      - the D_LR asymmetry score,
      - the log-likelihood ratio test statistic and p-value.

    Args:
        data (DataFrame): Ternary-coordinated data with columns T1, T2, T3.
        granularity (str or float): Grid resolution (e.g., 'fine', '0.1').
        file_name (str): Prefix used for saving any outputs (currently unused).

    Returns:
        DataFrame: Table of triangles and their statistics.
    """

    # Map granularity names to alpha values, or parse user-provided float
    if granularity == "superfine":
        alpha = 0.05
    elif granularity == "fine":
        alpha = 0.1
    elif granularity == "coarse":
        alpha = 0.25
    else:
        alpha = float(granularity)  # if the user provides a float value for granularity

    all_results = []

    # Loop over rows of triangles along T1 axis
    for row_index in range(int(1 / alpha)):
        a1 = row_index * alpha
        b1 = (row_index + 1) * alpha

        k_T2 = 0
        k_T3 = row_index

        while True:
            # Process up-triangle
            a2_up = k_T2 * alpha
            b2_up = (k_T2 + 1) * alpha
            a3_up = 1 - (k_T3 + 1) * alpha
            b3_up = 1 - k_T3 * alpha

            triangle_x, _, _ = return_triangle_coord(a1, b1, a2_up, b2_up, a3_up, b3_up)
            if round(triangle_x[0], 4) < 0:
                break  # Exit if triangle crosses outside the valid domain

            n_r, n_l = n(a1, b1, a2_up, b2_up, a3_up, b3_up, data)
            d_lr = D_LR(n_r, n_l)
            g_stat, p_val = log_likelihood_ratio_test(n_r, n_l)

            coords = [
                (round(a1, 4), round(b1, 4)),
                (round(a2_up, 4), round(b2_up, 4)),
                (round(a3_up, 4), round(b3_up, 4)),
            ]

            all_results.append([coords, n_r, n_l, d_lr, g_stat, p_val])

            # Process down-triangle (sharing base with the up-triangle)
            k_T3 += 1

            a2_down = k_T2 * alpha
            b2_down = (k_T2 + 1) * alpha
            a3_down = 1 - (k_T3 + 1) * alpha
            b3_down = 1 - k_T3 * alpha

            triangle_x, _, _ = return_triangle_coord(
                a1, b1, a2_down, b2_down, a3_down, b3_down
            )
            if round(triangle_x[0], 4) < 0:
                break

            n_r, n_l = n(a1, b1, a2_down, b2_down, a3_down, b3_down, data)
            d_lr = D_LR(n_r, n_l)
            g_stat, p_val = log_likelihood_ratio_test(n_r, n_l)

            coords = [
                (round(a1, 4), round(b1, 4)),
                (round(a2_down, 4), round(b2_down, 4)),
                (round(a3_down, 4), round(b3_down, 4)),
            ]

            all_results.append([coords, n_r, n_l, d_lr, g_stat, p_val])

            k_T2 += 1  # Advance to the next pair of up/down-triangles

    # Create DataFrame of results
    triangles = pd.DataFrame(
        all_results,
        columns=[
            "coord. (T1, T2, T3)",
            "n_right",
            "n_left",
            "D-LR",
            "g-test",
            "p-value(g-test)",
        ],
    )

    # Assign descending index for plotting (bottom-up row indexing)
    triangles["index"] = list(range(number_triangles(alpha), 0, -1))

    return triangles
