#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from .core import (T1, T2, T3, T1_lim, T2_lim, T3_lim, T3_lim_symm,
                  cartizian, return_triangle_coord, mid_point_triangle, h)

# Set up matplotlib
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
import matplotlib
matplotlib.rcParams.update(matplotlib.rcParamsDefault)

def plot(data, alpha, file_name):
    """Initial plotting of data-points in ternary coordinates"""
    fig = plt.figure(figsize=(8, 6))
    ax = plt.axes()

    x_side_T2 = np.linspace(0, 0.5, 100)
    x_side_T3 = np.linspace(-0.5, 0, 100)

    ax.plot(x_side_T2, T2(0,x_side_T2), "k", linewidth=1)
    ax.plot(x_side_T3, T3(0,x_side_T3), "k", linewidth=1)
    plt.hlines(y=0, xmin=-0.5, xmax=0.5, color="k", linewidth=1)
    
    ax.set_xticks([])
    ax.set_yticks([])

    for i in range(1, int(1/alpha)):
        y = i*alpha
        plt.hlines(y=y*h, xmin=T1_lim(y)[0], xmax=T1_lim(y)[1], color="crimson", linewidth=1)

        x2 = np.linspace(T2_lim(y)[0], T2_lim(y)[1], 100)
        ax.plot(x2, T2(y,x2), "dodgerblue", linewidth=1)

        x3 = np.linspace(T3_lim(y)[0], T3_lim(y)[1], 100)
        ax.plot(x3, T3(y,x3), "gold", linewidth=1)

    plt.vlines(x=0, ymin=0, ymax=h, colors="grey", ls=':')

    x_data = cartizian(data["T1"], data["T2"], data["T3"])[0]
    y_data = cartizian(data["T1"], data["T2"], data["T3"])[1]
    
    plt.scatter(x_data, y_data, color="lightsteelblue", alpha=0.5, s=9)
    
    plt.text(-0.02, 0.88, 'T1', size=12, color="crimson")
    plt.text(0.54, -0.01, 'T3', size=12, color="darkgoldenrod")
    plt.text(-0.58, -0.01, 'T2', size=12, color="dodgerblue")
    
    coord = np.arange(0, 1+alpha, alpha)
    T_1 = np.arange(0, 0.5+alpha/2, alpha/2)
    T_2 = np.arange(-0.5, 0+alpha/2, alpha/2)
    T_3 = np.arange(-0.5, 0.5+alpha, alpha)
        
    for i in range(len(T_1)):
        label = str(round(1-coord[i],2))
        x = T_1[i]
        y = T2(0,x)
        plt.text(x+0.01, y, label, size=7, color="crimson")
    
    for i in range(len(T_2)):
        label = str(round(1-coord[i],2))
        x = T_2[i]
        y = T3(0,x)
        plt.text(x-0.04, y, label, size=7, color="dodgerblue")  
    
    for i in range(len(coord)):
        label = str(round(coord[i],2))
        x = T_3[i]
        plt.text(x, -0.03, label, size=7, color="darkgoldenrod")   
        
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['top'].set_color('none')
     
    title = file_name + "_granuality_" + str(alpha) + ".png"  
    plt.savefig(title)
    return fig

def plot_results(res, granuality, file_name):
    """Plot analysis results"""
    alpha = granuality
    if granuality == "superfine":
        alpha = 0.05
    if granuality == "fine":
        alpha = 0.1
    if granuality == "coarse":
        alpha = 0.25

    fig = plt.figure(figsize=(5, 8))
    ax = plt.axes()

    x_side_T2 = np.linspace(0, 0.5, 100)
    ax.plot(x_side_T2, T2(0,x_side_T2), "k", linewidth=1)
    plt.hlines(y=0, xmin=0, xmax=0.5, color="k", linewidth=1)
    ax.set_xticks([])
    ax.set_yticks([])

    for i in range(1, int(1/(2*alpha))):
        y = i*alpha
        x2 = np.linspace(0, T2_lim(y)[1], 100)
        ax.plot(x2, T2(y,x2), "dodgerblue", linewidth=1)

    for i in range(1, int(1/alpha)):
        y = i*alpha
        plt.hlines(y=y*h, xmin=0, xmax=T1_lim(y)[1], color="crimson", linewidth=1)
        x3 = np.linspace(T3_lim_symm(y)[0], T3_lim_symm(y)[1], 100)
        ax.plot(x3, T3(y,x3), "gold", linewidth=1)

    plt.vlines(x=0, ymin=0, ymax=h, colors="grey", ls=':')

    for i in range(res["D-LR"].size):
        a1 = res["coord. (T1, T2, T3)"][i][0][0]
        b1 = res["coord. (T1, T2, T3)"][i][0][1]
        a2 = res["coord. (T1, T2, T3)"][i][1][0]
        b2 = res["coord. (T1, T2, T3)"][i][1][1]
        a3 = res["coord. (T1, T2, T3)"][i][2][0]
        b3 = res["coord. (T1, T2, T3)"][i][2][1]

        trianglex, triangley, direction = return_triangle_coord(a1, b1, a2, b2, a3, b3)

        if np.isnan(res["D-LR"][i]):
            plt.fill(trianglex, triangley, color="black")
        else:
            d_lr_color_score = (res["D-LR"][i] + 1)/2
            plt.fill(trianglex, triangley, color=(d_lr_color_score, 1-d_lr_color_score, 1-d_lr_color_score))

            x, y = mid_point_triangle(a1, b1, a2, b2, a3, b3)
            p = res["p-value(g-test)"][i]
            if p < 0.05 and p >= 0.001:
                ax.scatter(x, y, color='yellow', marker="*", alpha=0.4, s=9)
            if p < 0.001 and p >= 10**(-5):
                ax.scatter(x, y, color="darkslateblue", marker="*", alpha=0.9, s=22)
            if p <= 10**(-5):
                ax.scatter(x, y, color="black", marker="*", alpha=1, s=25)

    # Legends
    d_lr_color_score = 0
    pp1 = plt.Rectangle((0.2, 0.85), 0.05, 0.05, color=(d_lr_color_score, 1-d_lr_color_score, 1-d_lr_color_score))
    d_lr_color_score = 0.5
    pp2 = plt.Rectangle((0.25, 0.85), 0.05, 0.05, color=(d_lr_color_score, 1-d_lr_color_score, 1-d_lr_color_score))
    d_lr_color_score = 1
    pp3 = plt.Rectangle((0.3, 0.85), 0.05, 0.05, color=(d_lr_color_score, 1-d_lr_color_score, 1-d_lr_color_score))
    pp4 = plt.Rectangle((0.518, 0.85), 0.05, 0.05, color="black")

    ax.add_patch(pp1)
    ax.add_patch(pp2)
    ax.add_patch(pp3)
    ax.add_patch(pp4)
    plt.text(0.16, 0.865, "$D_{lr} =      -1$", size=8)
    plt.text(0.26, 0.865, '0', size=8)
    plt.text(0.31, 0.865, '1', size=8)
    plt.text(0.38, 0.865, 'empty triangle', size=8)

    ax.scatter(0.2, 0.8, color="black", marker="*", alpha=1, s=25)
    plt.text(0.214, 0.8, "$p < 10^{-5}$", size=8)
    ax.scatter(0.2, 0.77, color="darkslateblue", marker="*", alpha=0.9, s=22)
    plt.text(0.214, 0.77, "$p < 0.001$", size=8)
    ax.scatter(0.2, 0.74, color="darkgoldenrod", marker="*", alpha=0.4, s=9)
    plt.text(0.214, 0.74, "$p < 0.05$", size=8)

    plt.text(-0.02, 0.88, 'T1', size=12, color="crimson")
    plt.text(-0.03, -0.01, 'T2', size=12, color="dodgerblue")
    plt.text(0.54, -0.01, 'T3', size=12, color="darkgoldenrod")

    coord = np.arange(0, 1+alpha, alpha)
    T_1 = np.arange(0, 0.5+alpha/2, alpha/2)
    T_2_3 = np.arange(0, 0.5+alpha, alpha)

    for i in range(len(T_1)):
        label = str(round(1-coord[i],2))
        x = T_1[i]
        y = T2(0,x)
        plt.text(x+0.01, y, label, size=7, color="firebrick")

    for i in range(len(T_2_3)-2):
        label = str(round(coord[i]+alpha,2))
        x = 0
        y = T2(coord[i]+alpha,0)
        plt.text(x-0.033, y, label, size=7, color="dodgerblue")  

    for i in range(len(T_2_3)):
        label = str(round(coord[i]+0.5,2))
        x = T_2_3[i]
        plt.text(x, -0.03, label, size=7, color="darkgoldenrod")   

    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['top'].set_color('none')

    title = file_name + "_analysis_granuality_" + str(alpha) + ".png"  
    plt.savefig(title)

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

def plotting_triangle_index(res, granuality):
    """Plot triangle index"""
    alpha = granuality
    fig_size = (5, 4)
    font_size = 8

    if granuality == "superfine":
        alpha = 0.05
        fig_size = (7, 6)
        font_size = 7
    elif granuality == "fine":
        alpha = 0.1
        fig_size = (5, 4)
        font_size = 8
    elif granuality == "coarse":
        alpha = 0.25
        fig_size = (4, 3)
        font_size = 9
    
    fig = plt.figure(figsize=fig_size)
    ax = plt.axes()

    x_side_T2 = np.linspace(0, 0.5, 100)
    ax.plot(x_side_T2, T2(0,x_side_T2), "k", linewidth=1)
    plt.hlines(y=0, xmin=0, xmax=0.5, color="k", linewidth=1)
    ax.set_xticks([])
    ax.set_yticks([])

    for i in range(1, int(1/(2*alpha))):
        y = i*alpha
        x2 = np.linspace(0, T2_lim(y)[1], 100)
        ax.plot(x2, T2(y,x2), "dodgerblue", linewidth=1)

    for i in range(1, int(1/alpha)):
        y = i*alpha
        plt.hlines(y=y*h, xmin=0, xmax=T1_lim(y)[1], color="crimson", linewidth=1)
        x3 = np.linspace(T3_lim_symm(y)[0], T3_lim_symm(y)[1], 100)
        ax.plot(x3, T3(y,x3), "gold", linewidth=1)

    plt.vlines(x=0, ymin=0, ymax=h, colors="grey", ls=':')

    for i in range(res["D-LR"].size):
        a1 = res["coord. (T1, T2, T3)"][i][0][0]
        b1 = res["coord. (T1, T2, T3)"][i][0][1]
        a2 = res["coord. (T1, T2, T3)"][i][1][0]
        b2 = res["coord. (T1, T2, T3)"][i][1][1]
        a3 = res["coord. (T1, T2, T3)"][i][2][0]
        b3 = res["coord. (T1, T2, T3)"][i][2][1]

        x, y = mid_point_triangle(a1, b1, a2, b2, a3, b3)
        index = res["index"][i]
        plt.text(x-0.01, y, str(index), size=font_size)

    plt.text(-0.02, 0.88, 'T1', size=12, color="crimson")
    plt.text(-0.03, -0.01, 'T2', size=12, color="dodgerblue")
    plt.text(0.54, -0.01, 'T3', size=12, color="darkgoldenrod")

    coord = np.arange(0, 1+alpha, alpha)
    T_1 = np.arange(0, 0.5+alpha/2, alpha/2)
    T_2_3 = np.arange(0, 0.5+alpha, alpha)

    for i in range(len(T_1)):
        label = str(round(1-coord[i],2))
        x = T_1[i]
        y = T2(0,x)
        plt.text(x+0.01, y, label, size=7, color="firebrick")

    for i in range(len(T_2_3)-2):
        label = str(round(coord[i]+alpha,2))
        x = 0
        y = T2(coord[i]+alpha,0)
        plt.text(x-0.033, y, label, size=7, color="dodgerblue")  

    for i in range(len(T_2_3)):
        label = str(round(coord[i]+0.5,2))
        x = T_2_3[i]
        plt.text(x, -0.03, label, size=7, color="darkgoldenrod")   

    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['top'].set_color('none')

    title = "index_granulality_" + str(alpha) + ".png"  
    plt.savefig(title)
    return fig 