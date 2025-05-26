#!/usr/bin/env python
# coding: utf-8

import numpy as np
from scipy.stats import binom
from scipy.stats.distributions import chi2
from math import log
import pandas as pd
from .core import cartizian, return_triangle_coord, number_triangles

def fundemental_asymmetry(data):
    """Implement basic symmetry analysis between the two primary subtriangles"""
    data["x-axis"] = cartizian(data["T1"], data["T2"], data["T3"])[0]
    main_n_r = len(data[data["x-axis"] > 0])
    main_n_l = len(data[data["x-axis"] < 0])

    main_d_lr = D_LR(main_n_r, main_n_l)
    main_g_test, main_p_value = log_likelihood_ratio_test(main_n_r, main_n_l)
    return (main_n_r, main_n_l, main_d_lr, main_g_test, main_p_value)

def n(a1, b1, a2, b2, a3, b3, data):
    """Count points in triangle and its reflection"""
    if a1 == 0:
        condition_a1 = a1 <= data.T1
    else:
        condition_a1 = a1 < data.T1
    
    if a2 == 0:
        condition_a2 = a2 <= data.T2
    else:
        condition_a2 = a2 < data.T2
        
    if a3 == 0:
        condition_a3 = a3 <= data.T3
    else:
        condition_a3 = a3 < data.T3
        
    n_1 = len(data[((a1 <= data.T1) & (data.T1 <= b1)) &
                   ((a2 <= data.T2) & (data.T2 <= b2)) &
                   ((a3 <= data.T3) & (data.T3 <= b3))])
    
    n_2 = len(data[((a1 <= data.T1) & (data.T1 <= b1)) &
                   ((a3 <= data.T2) & (data.T2 <= b3)) &
                   ((a2 <= data.T3) & (data.T3 <= b2))])
    
    trianglex, triangley, direction = return_triangle_coord(a1, b1, a2, b2, a3, b3)
    top_node = trianglex[1]
    
    if top_node > 0:
        n_r = n_1
        n_l = n_2
    else:
        n_r = n_2
        n_l = n_1
    
    return (n_r, n_l)

def D_LR(n_r, n_l):
    """Calculate d_lr value"""
    if n_l + n_r != 0:
        d_lr = (n_r - 0.5*(n_l + n_r))/(0.5*(n_l + n_r))
    else:
        d_lr = np.nan
        
    return d_lr

def log_likelihood_ratio_test(n_r, n_l):
    """Calculate log likelihood ratio test"""
    N = n_r + n_l
    if N != 0:
        L_null = binom.pmf(n_l, N, 0.5)
        L_alt = binom.pmf(n_l, N, n_l/N)

        if L_null/L_alt < 10**-16:
            test_res = np.nan
            p_value = 0
        else:
            test_res = float(-2 * log(L_null/L_alt))
            p_value = 1 - chi2.cdf(test_res, 1)
    else:
        test_res = np.nan
        p_value = np.nan

    return (test_res, p_value)

def triangles_analysis(data, granularity, file_name):
    """Perform triangle analysis"""
    if isinstance(granularity, str):
        if granularity == "superfine":
            alpha = 0.05
        elif granularity == "fine":
            alpha = 0.1
        elif granularity == "coarse":
            alpha = 0.25
        else:
            raise ValueError("Invalid granularity string. Must be one of: 'coarse', 'fine', 'superfine'")
    else:
        alpha = float(granularity)

    tri = []

    for i in range(int(1/alpha)):
        a1 = i * alpha
        b1 = (i+1) * alpha 

        k2 = 0
        k3 = 0 + i

        trianglex = [1,1,1,1]
        while trianglex[0] >= 0:
            a2 = k2 * alpha
            b2 = (k2+1) * alpha

            a3 = 1-(k3+1) * alpha
            b3 = 1 - k3 * alpha

            trianglex, triangley, direction = return_triangle_coord(a1, b1, a2, b2, a3, b3)
            if round(trianglex[0],4) < 0:
                continue

            n_r, n_l = n(a1, b1, a2, b2, a3, b3, data)
            d_lr = D_LR(n_r, n_l)
            g_test, p_value = log_likelihood_ratio_test(n_r, n_l)

            coord = [(round(a1,4), round(b1,4)),
                    (round(a2,4), round(b2,4)),
                    (round(a3,4), round(b3,4))]
            tri.append([coord, n_r, n_l, d_lr, g_test, p_value])

            k3 = k3 + 1

            a2 = k2 * alpha
            b2 = (k2+1) * alpha

            a3 = 1-(k3+1) * alpha
            b3 = 1 - k3 * alpha

            trianglex, triangley, direction = return_triangle_coord(a1, b1, a2, b2, a3, b3)
            if round(trianglex[0],4) < 0:
                continue

            n_r, n_l = n(a1, b1, a2, b2, a3, b3, data)
            d_lr = D_LR(n_r, n_l)
            g_test, p_value = log_likelihood_ratio_test(n_r, n_l)

            coord = [(round(a1,4), round(b1,4)),
                    (round(a2,4), round(b2,4)),
                    (round(a3,4), round(b3,4))]
            tri.append([coord, n_r, n_l, d_lr, g_test, p_value])

            k2 = k2 + 1

    triangles = pd.DataFrame(tri)
    triangles.columns = ["coord. (T1, T2, T3)", "n_right", "n_left", "D-LR", "g-test", "p-value(g-test)"]
    index = list(range(number_triangles(alpha), 0, -1))
    triangles["index"] = index
    
    return triangles 