#!/usr/bin/env python
# coding: utf-8

import os
from pathlib import Path
import numpy as np
from scipy.stats import binom
from scipy.stats.distributions import chi2
from math import log
import pandas as pd
import matplotlib.pyplot as plt
from twisstntern.core import (cartizian, return_triangle_coord, dump_data,
                            n, D_LR, log_likelihood_ratio_test, fundemental_asymmetry)
from twisstntern.visualization import (plot, plot_results, plotting_triangle_index, plot_fundemental_asymmetry)

def number_triangles(alpha):
    """Calculate number of triangles for given granularity"""
    if int(1/alpha) % 2:
        raise ValueError("1/alpha must be even")
    
    a = int((1-2*alpha)/(2*alpha))
    n = int(1/(2*alpha))
    for i in range(1, a+1):
        n = n + 4*int((1/(2*alpha))-i)
    
    return n

def triangles_analysis(data, granularity, file_name):
    """Perform triangle analysis"""
    alpha = granularity
    if granularity == "superfine":
        alpha = 0.05
    elif granularity == "fine":
        alpha = 0.1
    elif granularity == "coarse":
        alpha = 0.25
    
    n_triangles = number_triangles(alpha)
    
    results = {
        "D-LR": np.zeros(n_triangles),
        "p-value(g-test)": np.zeros(n_triangles),
        "coord. (T1, T2, T3)": [],
        "index": np.arange(n_triangles)
    }
    
    triangle_index = 0
    
    # Analyze each triangle
    for i in range(int(1/(2*alpha))):
        for j in range(int(1/(2*alpha))-i):
            # Up triangle
            a1 = i*alpha
            b1 = (i+1)*alpha
            a2 = j*alpha
            b2 = (j+1)*alpha
            a3 = (j+1)*alpha
            b3 = (j+2)*alpha
            
            n_r = n(a1, b1, a2, b2, a3, b3, data)
            n_l = n(a1, b1, -b3, -a3, -b2, -a2, data)
            
            d_lr = D_LR(n_r, n_l)
            g_test, p_value = log_likelihood_ratio_test(n_r, n_l)
            
            results["D-LR"][triangle_index] = d_lr
            results["p-value(g-test)"][triangle_index] = p_value
            results["coord. (T1, T2, T3)"].append([(a1, b1), (a2, b2), (a3, b3)])
            
            triangle_index += 1
            
            # Down triangle
            if i < int(1/(2*alpha))-1:
                a1 = (i+1)*alpha
                b1 = (i+2)*alpha
                a2 = j*alpha
                b2 = (j+1)*alpha
                a3 = (j+1)*alpha
                b3 = (j+2)*alpha
                
                n_r = n(a1, b1, a2, b2, a3, b3, data)
                n_l = n(a1, b1, -b3, -a3, -b2, -a2, data)
                
                d_lr = D_LR(n_r, n_l)
                g_test, p_value = log_likelihood_ratio_test(n_r, n_l)
                
                results["D-LR"][triangle_index] = d_lr
                results["p-value(g-test)"][triangle_index] = p_value
                results["coord. (T1, T2, T3)"].append([(a1, b1), (a2, b2), (a3, b3)])
                
                triangle_index += 1
    
    return results

def run_analysis(file, granularity):
    """Run full analysis pipeline"""
    # Create Results directory
    results_dir = Path("Results")
    results_dir.mkdir(exist_ok=True)
    
    # Load and process data
    data = dump_data(file)
    
    # Run analyses
    results = triangles_analysis(data, granularity, file)
    fundamental_results = fundemental_asymmetry(data)
    
    # Generate all plots
    plot_fundemental_asymmetry(data, file)
    plot(data, granularity, file)
    plot_results(results, file)
    plotting_triangle_index(results, granularity)
    
    return results, fundamental_results 