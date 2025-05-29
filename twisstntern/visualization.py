#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import pandas as pd
from math import sqrt
from twisstntern.core import (cartizian, return_triangle_coord, dump_data,
                            T1, T2, T3, T1_lim, T2_lim, T3_lim, T3_lim_symm, h,
                            fundemental_asymmetry)

def plot_fundemental_asymmetry(data, file_name):
    """Plot fundamental asymmetry analysis"""
    fig = plt.figure(figsize=(6, 4))
    ax = plt.axes()

    x_side_T2 = np.linspace(0, 0.5, 100)
    x_side_T3 = np.linspace(-0.5, 0, 100)

    ax.plot(x_side_T2, T2(0,x_side_T2), "k", linewidth=1)
    ax.plot(x_side_T3, T3(0,x_side_T3), "k", linewidth=1)
    plt.hlines(y=0, xmin=-0.5, xmax=0.5, color="k", linewidth=1)
    
    ax.set_xticks([])
    ax.set_yticks([])

    plt.vlines(x=0, ymin=0, ymax=h, colors="black")

    trianglex_R = [0, 0, 0.5, 0]
    triangley_R = [0, h, 0, 0]

    trianglex_L = [-0.5, 0, 0, -0.5]
    triangley_L = [0, h, 0, 0]

    main_n_r, main_n_l, main_d_lr, main_g_test, main_p_value = fundemental_asymmetry(data)

    main_d_lr_left = (main_n_l-0.5*(main_n_l+main_n_r))/(0.5*(main_n_l+main_n_r))

    d_lr_color_score_R = (main_d_lr + 1)/2
    d_lr_color_score_L = (main_d_lr_left + 1)/2

    plt.fill(trianglex_R, triangley_R, color=(d_lr_color_score_R, 1-d_lr_color_score_R, 1-d_lr_color_score_R))
    plt.fill(trianglex_L, triangley_L, color=(d_lr_color_score_L, 1-d_lr_color_score_L, 1-d_lr_color_score_L))

    x = 0.15
    y = 0.4*h
    p = main_p_value
    if p < 0.05 and p >= 0.001:
        ax.scatter(x, y, color='yellow', marker="*", alpha=0.4, s=9)
    if p < 0.001 and p >= 10**(-5):
        ax.scatter(x, y, color="darkslateblue", marker="*", alpha=0.9, s=22)
    if p <= 10**(-5):
        ax.scatter(x, y, color="black", marker="*", alpha=1, s=25)

    d_lr_color_score = 0
    pp1 = plt.Rectangle((0.2, 0.85), 0.05, 0.05, color=(d_lr_color_score, 1-d_lr_color_score, 1-d_lr_color_score))
    d_lr_color_score = 0.5
    pp2 = plt.Rectangle((0.25, 0.85), 0.05, 0.05, color=(d_lr_color_score, 1-d_lr_color_score, 1-d_lr_color_score))
    d_lr_color_score = 1
    pp3 = plt.Rectangle((0.3, 0.85), 0.05, 0.05, color=(d_lr_color_score, 1-d_lr_color_score, 1-d_lr_color_score))
    ax.add_patch(pp1)
    ax.add_patch(pp2)
    ax.add_patch(pp3)
    plt.text(0.135, 0.865, "$D_{lr} =-1$", size=7)
    plt.text(0.26, 0.865, '0', size=7)
    plt.text(0.31, 0.865, '1', size=7)
    
    ax.scatter(0.2, 0.8, color="black", marker="*", alpha=1, s=25)
    plt.text(0.214, 0.8, "$p < 10^{-5}$", size=8)
    ax.scatter(0.2, 0.77, color="darkslateblue", marker="*", alpha=0.9, s=22)
    plt.text(0.214, 0.77, "$p < 0.001$", size=8)
    ax.scatter(0.2, 0.74, color="darkgoldenrod", marker="*", alpha=0.4, s=9)
    plt.text(0.214, 0.74, "$p < 0.05$", size=8)

    title_left = str(main_n_l)
    title_right = str(main_n_r)
    plt.text(-0.5, -0.1, " n =", size=12)
    plt.text(-0.3, -0.1, title_left, size=12, color="grey")
    plt.text(0.2, -0.1, title_right, size=12, color="grey")

    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['top'].set_color('none')

    title = file_name + "_fundamental_asymmetry.png"  
    plt.savefig(title)
    return (main_d_lr, main_g_test, main_p_value)

def plot(data, alpha, file_name):
    """Plot ternary coordinates"""
    fig = plt.figure(figsize=(6, 4))
    ax = plt.axes()

    x_side_T2 = np.linspace(0, 0.5, 100)
    x_side_T3 = np.linspace(-0.5, 0, 100)

    ax.plot(x_side_T2, T2(0,x_side_T2), "k", linewidth=1)
    ax.plot(x_side_T3, T3(0,x_side_T3), "k", linewidth=1)
    plt.hlines(y=0, xmin=-0.5, xmax=0.5, color="k", linewidth=1)
    
    ax.set_xticks([])
    ax.set_yticks([])

    plt.vlines(x=0, ymin=0, ymax=h, colors="black")

    for i in range(data.shape[0]):
        x = cartizian(data.iloc[i,0], data.iloc[i,1], data.iloc[i,2])[0]
        y = cartizian(data.iloc[i,0], data.iloc[i,1], data.iloc[i,2])[1]
        ax.scatter(x, y, color="black", alpha=0.1, s=1)

    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['top'].set_color('none')

    title = file_name + "_granuality_" + str(alpha) + ".png"  
    plt.savefig(title)

def plot_results(results, file_name):
    """Plot results summary"""
    fig = plt.figure(figsize=(6, 4))
    ax = plt.axes()

    x_side_T2 = np.linspace(0, 0.5, 100)
    x_side_T3 = np.linspace(-0.5, 0, 100)

    ax.plot(x_side_T2, T2(0,x_side_T2), "k", linewidth=1)
    ax.plot(x_side_T3, T3(0,x_side_T3), "k", linewidth=1)
    plt.hlines(y=0, xmin=-0.5, xmax=0.5, color="k", linewidth=1)
    
    ax.set_xticks([])
    ax.set_yticks([])

    plt.vlines(x=0, ymin=0, ymax=h, colors="black")

    for i in range(len(results["coord. (T1, T2, T3)"])):
        a1 = results["coord. (T1, T2, T3)"][i][0][0]
        b1 = results["coord. (T1, T2, T3)"][i][0][1]
        a2 = results["coord. (T1, T2, T3)"][i][1][0]
        b2 = results["coord. (T1, T2, T3)"][i][1][1]
        a3 = results["coord. (T1, T2, T3)"][i][2][0]
        b3 = results["coord. (T1, T2, T3)"][i][2][1]
        
        trianglex, triangley, direction = return_triangle_coord(a1, b1, a2, b2, a3, b3)
        
        d_lr = results["D-LR"][i]
        d_lr_color_score = (d_lr + 1)/2
        
        plt.fill(trianglex, triangley, color=(d_lr_color_score, 1-d_lr_color_score, 1-d_lr_color_score))
        
        p = results["p-value(g-test)"][i]
        if p < 0.05 and p >= 0.001:
            ax.scatter(trianglex[0], triangley[0], color='yellow', marker="*", alpha=0.4, s=9)
        if p < 0.001 and p >= 10**(-5):
            ax.scatter(trianglex[0], triangley[0], color="darkslateblue", marker="*", alpha=0.9, s=22)
        if p <= 10**(-5):
            ax.scatter(trianglex[0], triangley[0], color="black", marker="*", alpha=1, s=25)

    d_lr_color_score = 0
    pp1 = plt.Rectangle((0.2, 0.85), 0.05, 0.05, color=(d_lr_color_score, 1-d_lr_color_score, 1-d_lr_color_score))
    d_lr_color_score = 0.5
    pp2 = plt.Rectangle((0.25, 0.85), 0.05, 0.05, color=(d_lr_color_score, 1-d_lr_color_score, 1-d_lr_color_score))
    d_lr_color_score = 1
    pp3 = plt.Rectangle((0.3, 0.85), 0.05, 0.05, color=(d_lr_color_score, 1-d_lr_color_score, 1-d_lr_color_score))
    ax.add_patch(pp1)
    ax.add_patch(pp2)
    ax.add_patch(pp3)
    plt.text(0.135, 0.865, "$D_{lr} =-1$", size=7)
    plt.text(0.26, 0.865, '0', size=7)
    plt.text(0.31, 0.865, '1', size=7)
    
    ax.scatter(0.2, 0.8, color="black", marker="*", alpha=1, s=25)
    plt.text(0.214, 0.8, "$p < 10^{-5}$", size=8)
    ax.scatter(0.2, 0.77, color="darkslateblue", marker="*", alpha=0.9, s=22)
    plt.text(0.214, 0.77, "$p < 0.001$", size=8)
    ax.scatter(0.2, 0.74, color="darkgoldenrod", marker="*", alpha=0.4, s=9)
    plt.text(0.214, 0.74, "$p < 0.05$", size=8)

    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['top'].set_color('none')

    title = file_name + "_triangle_analysis.png"  
    plt.savefig(title)

def plotting_triangle_index(results, granularity):
    """Plot triangle index"""
    fig = plt.figure(figsize=(6, 4))
    ax = plt.axes()

    x_side_T2 = np.linspace(0, 0.5, 100)
    x_side_T3 = np.linspace(-0.5, 0, 100)

    ax.plot(x_side_T2, T2(0,x_side_T2), "k", linewidth=1)
    ax.plot(x_side_T3, T3(0,x_side_T3), "k", linewidth=1)
    plt.hlines(y=0, xmin=-0.5, xmax=0.5, color="k", linewidth=1)
    
    ax.set_xticks([])
    ax.set_yticks([])

    plt.vlines(x=0, ymin=0, ymax=h, colors="black")

    for i in range(len(results["coord. (T1, T2, T3)"])):
        a1 = results["coord. (T1, T2, T3)"][i][0][0]
        b1 = results["coord. (T1, T2, T3)"][i][0][1]
        a2 = results["coord. (T1, T2, T3)"][i][1][0]
        b2 = results["coord. (T1, T2, T3)"][i][1][1]
        a3 = results["coord. (T1, T2, T3)"][i][2][0]
        b3 = results["coord. (T1, T2, T3)"][i][2][1]
        
        trianglex, triangley, direction = return_triangle_coord(a1, b1, a2, b2, a3, b3)
        
        plt.fill(trianglex, triangley, color="white")
        ax.text(trianglex[0], triangley[0], str(i), size=5)

    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['top'].set_color('none')

    title = "index_granulality_" + str(granularity) + ".png"  
    plt.savefig(title) 