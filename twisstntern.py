#!/usr/bin/env python
# coding: utf-8

# In[1]:


# libraries
import os
from pathlib import Path

import numpy as np
import random
from math import *
import scipy

from scipy.stats import binom
from scipy.stats.distributions import chi2
import matplotlib.pyplot as plt
import pandas as pd

from sympy import Eq, Symbol as sym, solve
import sys

# presentation in the dataframe
pd.set_option('display.float_format', '{:.4e}'.format)
# activating latex printing in matplotlib
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})
import matplotlib
matplotlib.rcParams.update(matplotlib.rcParamsDefault)

import sys 


# # Preliminaries - Ternary setup

# ### Coordinates

# In[2]:


h = sqrt(3)/2 # = the hight of an equilateral triangle, where each side = 1


# Formulate the appropriate isoclines within the ternary triangle
# for a given y-coordinate value in the range [0, 1] of Ti

def T1(y,x): # -
    return y*h

def T2(y,x): # \
    a = 2*h
    b = sqrt(3)*(0.5-y) 
    return -a*x + b

def T3(y,x): # /
    a = 2*h
    b = sqrt(3)*(0.5-y) 
    return a*x + b


# Convert a Cartesian coordinate in R^2 to ternary coordinates.
# 1st coordinate is T_1 == -
# 2nd coordinate is T_2 == \
# 3rd coordinate is T_3 == /

def ternary_coord(a,b): # converting Cartisian -> ternary coordinates
    t1=b/h
    
    y = sym('y')
    eqa_2 = Eq(T2(y,a), b)
    t2 = solve(eqa_2,y)[0]
    
    z = sym('z')
    eqa_3 = Eq(T3(z,a), b)
    t3 = solve(eqa_3,z)[0]
    
    return (t1,t2,t3)


def cartizian(x,y,z): # converting ternary -> Cartisian coordinates
    return ((z-y)/2, x*h)


# In[3]:


# verification- it works :)

# x=random.random()
# y=random.random()

# print(cartizian(ternary_coord(x,y)[0],ternary_coord(x,y)[1],ternary_coord(x,y)[2])[0]-x)
# print(cartizian(ternary_coord(x,y)[0],ternary_coord(x,y)[1],ternary_coord(x,y)[2])[1]-y)


# ## Plotting utility functions

# In[4]:


# Utility functions for plotting.

# Return the x-axis limits for accurately sketching Ti(y,.) within our triangle.


def T1_lim(y): # - given y coordinate of T1, returns the delimeters [x_l,x_r] to plot it in 
    
    # The general logic: y =T3(0,x_left)= T2(1, x_right)
    
    y=y*h # the real height of coordinate y in T1 (as it is all scaled by h)
    x_l = sym('x_l') # x_left
    eqa_l = Eq(T3(0,x_l),y)
    x_L = float(solve(eqa_l,x_l)[0])

    x_r = sym('x_r') # x_right
    eqa_r = Eq(T2(0,x_r),y)
    x_R = float(solve(eqa_r,x_r)[0])
    return (x_L,x_R) 



#full plot
def T2_lim(y): # Given the y-coordinate of T2, calculate the delimiters [x_l, x_r] for plotting it within the triangle.

    
    # Logic:  T3(0,x_left)= T2(y, x_left),  0= T2(y, x_right) should be =1-y-0.5=0.5-y
    
    x_l = sym('x_l') # x_left
    eqa_l = Eq(T3(0,x_l),T2(y, x_l))
    x_L = float(solve(eqa_l,x_l)[0])

    x_r = sym('x_r') # x_right
    eqa_r = Eq(T2(y,x_r),0)
    x_R = float(solve(eqa_r,x_r)[0])
    return (x_L,x_R) 

    
#full plot
def T3_lim(y): # Given the y-coordinate of T3, calculate the delimiters [x_l, x_r] for plotting it within the triangle.
    
    # Logic:   0= T3(y, x_left) should be = y-0.5 ,  T3(y,x_right)= T2(0, x_right)  
   
    x_l = sym('x_l') # x_left
    eqa_l = Eq(T3(y,x_l),0)
    x_L = float(solve(eqa_l,x_l)[0])

    x_r = sym('x_r') # x_right
    eqa_r = Eq(T3(y,x_r),T2(0, x_r) )
    x_R = float(solve(eqa_r,x_r)[0])
    return (x_L,x_R)

# plotting only the right half- for the results plot
def T3_lim_symm(y): # Given the y-coordinate of T3, calculate the delimiters [x_l, x_r] for plotting it within 
                    # right half of the triangle.

    
    # Logic:  T3(0,x_left)= T2(y, x_left),  0= T2(y, x_right) should be =1-y-0.5=0.5-y
   
    x_l = sym('x_l') # x_left
    eqa_l = Eq(T3(y,x_l),0)
    x_L = float(solve(eqa_l,x_l)[0])

    x_r = sym('x_r') # x_right
    eqa_r = Eq(T3(y,x_r),T2(0, x_r) )
    x_R = float(solve(eqa_r,x_r)[0])
    
    
    return (max(x_L,0),x_R) 


# In[5]:


# Given the coordinates of a subtriangle, returns two arrays in Cartesian coordinates,
# along with the triangle orientation as a string ("up" or "down").

def return_triangle_coord(a1,b1,a2,b2,a3,b3): 
    
    x_s = sym('x_s') # s- spitz node
    x_l = sym('x_l') # l- left node
    x_r = sym('x_r') # r- right node

    eqa_r = Eq(T2(a2, x_r) ,T3(b3,x_r)) # true for both up & down triangles
    x_r = float(solve(eqa_r,x_r)[0])

    eqa_l = Eq(T2(b2, x_l) ,T3(a3,x_l))   # true for both up & down triangles
    x_l = float(solve(eqa_l,x_l)[0])
    
 
    # A triangle is classified as an "up-triangle" if and only if the y-component, such as T2(b2, x_l),
    # of any of the base nodes intersects with b1/h, which represents the higher T1 coordinate.

    if abs(T2(b2,x_l)-a1*h) <= 10**-14: # T2(b2,x_l) == (a1*h)  
        direction = "up"  # think later if we want to store this as strings really, maybe not for some condition later-on
        
        eqa_s = Eq(T2(a2, x_s) ,T3(a3,x_s))
        x_s = float(solve(eqa_s,x_s)[0])

        triangley = [ a1*h, b1*h, a1*h, a1*h ]  # y coordinates of [ x_l,x_s ,x_r , x_l ] 
    
    
    if  abs(T2(b2,x_l)-b1*h) <= 10**-14 :# T2(b2,x_l) == (b1*h)  :
        direction = "down"
    
        eqa_s = Eq(T2(b2, x_s) ,T3(b3,x_s))
        x_s = float(solve(eqa_s,x_s)[0])
    
        triangley = [ b1*h, a1*h, b1*h, b1*h ]   # y coordinates of [ x_l,x_s ,x_r , x_l ] 

    
    trianglex = [ x_l, x_s ,x_r , x_l ] 
    return (trianglex, triangley, direction) # this format is ready to be printed with: plt.fill(trianglex, triangley)


# In[6]:


# reflected coordinate of the reflected triangle
def ref(a1,b1,a2,b2,a3,b3):
    return (a1,b1,a3,b3,a2,b2)


# In[7]:


# Given the coordinates of a (right) subtriangle, determine the midpoint's location
# to be utilized for plotting the p-value results within the sub-triangles.

def mid_point_triangle(a1,b1,a2,b2,a3,b3): # in cartesian coordinates
    trianglex, triangley, direction=return_triangle_coord(a1,b1,a2,b2,a3,b3)
    mid_x= trianglex[1]# x coordinate of the top
    
    if direction == "up":          
            mid_y= triangley[0] +(trianglex[1]-trianglex[0])/2 # simple geometry
    else: #    direction == "down"     
        mid_y= triangley[0] +2*(trianglex[0]-trianglex[1])/3 # simple geometry
        
    return (mid_x,round(mid_y, 4))    


# ### Plotting the Data

# In[8]:


def plot(data, alpha,file_name): # Initial plotting of data-points in ternary coordinates

    fig = plt.figure(figsize=(8, 6))
    ax = plt.axes()

    x_side_T2 = np.linspace(0, 0.5, 100)
    x_side_T3 = np.linspace(-0.5, 0, 100)

    ax.plot(x_side_T2, T2(0,x_side_T2),"k",linewidth=1);
    ax.plot(x_side_T3, T3(0,x_side_T3),"k",linewidth=1);
    plt.hlines(y = 0, xmin =-0.5, xmax = 0.5,color="k",linewidth=1)
    
    # Hide X and Y axes tick marks
    ax.set_xticks([])
    ax.set_yticks([])


    # coordinates!
    for i in range(1, int(1/alpha)):
        y=i*alpha
        # ploting T1 coordinates
        plt.hlines(y = y*h, xmin = T1_lim(y)[0], xmax = T1_lim(y)[1],color="crimson",linewidth=1)

        # ploting T2 coordinates
        x2 = np.linspace(T2_lim(y)[0], T2_lim(y)[1], 100)
        ax.plot(x2, T2(y,x2), "dodgerblue", linewidth=1)

        # ploting T3 coordinates
        x3 = np.linspace(T3_lim(y)[0], T3_lim(y)[1], 100)
        ax.plot(x3, T3(y,x3), "gold", linewidth=1)

    # single vline with full ymin and ymax
    plt.vlines(x=0, ymin=0, ymax=h, colors="grey", ls=':')

    x_data = cartizian(data["T1"],data["T2"],data["T3"])[0] # x coordinates of the data points to be plotted
    y_data = cartizian(data["T1"],data["T2"],data["T3"])[1] # y coordinates of the data points to be plotted
    
    plt.scatter(x_data,y_data, color = "lightsteelblue",alpha=0.5,s=9)
#     plt.scatter(x_data,y_data, color = "blue",alpha=1,s=15)
    
    plt.text(-0.02, 0.88, 'T1', size=12, color="crimson")
    plt.text(0.54, -0.01, 'T3', size=12, color= "darkgoldenrod")
    plt.text(-0.58, -0.01, 'T2', size=12, color = "dodgerblue")
    
    
    # printing the coordinates
    
    coord= np.arange(0,1+alpha,alpha)
    T_1 = np.arange(0,0.5+alpha/2,alpha/2)
    T_2 = np.arange(-0.5,0+alpha/2,alpha/2)
    T_3 = np.arange(-0.5,0.5+alpha,alpha)
        
    # coordintes T1
    for i in range(len(T_1)):
        label= str(round(1-coord[i],2))
        x=T_1[i]
        y=T2(0,x)
        plt.text(x+0.01, y, label , size=7, color = "crimson")
    
    # coordintes T2
    for i in range(len(T_2)):
        label= str(round(1-coord[i],2))
        x=T_2[i]
        y=T3(0,x)
        plt.text(x-0.04, y, label , size=7, color = "dodgerblue")  
    
    # coordintes T3    
    for i in range(len(coord)):
        label= str(round(coord[i],2))
        x=T_3[i]
        plt.text(x, -0.03, label , size=7, color = "darkgoldenrod")   
        
    # removing the box lines around the plot
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['top'].set_color('none')
     
    #saving the plot
    title = file_name + "_granuality_" + str(alpha)  +".png"  
    plt.savefig(title)
    return fig


# ### Plotting Results

# In[9]:


# Generates a plot of the analysis results using the 'res' DataFrame, which is the output of 
# 'triangles_analysis(data, granularity)'. 
# The granularity should be provided as one of the strings; "course", "fine" or "superfine".
# If a non-standard choice was used, specify the different granularity as a float.

def plot_results(res, granuality, file_name):

    # If the user didn't specify a granularity level, it implies they used one of the default settings.
    alpha = granuality 
    if granuality == "superfine":
        alpha = 0.05 # granulity "super fine" was used
    if granuality == "fine":
        alpha = 0.1 # granulity "fine" was used
    if granuality == "coarse":
        alpha = 0.25 # granulity "coarse" was used


    fig = plt.figure(figsize=(5, 8))
    ax = plt.axes()

    x_side_T2 = np.linspace(0, 0.5, 100)

    ax.plot(x_side_T2, T2(0,x_side_T2),"k",linewidth=1);

    plt.hlines(y = 0, xmin =0, xmax = 0.5,color="k",linewidth=1)
    # Hide X and Y axes tick marks
    ax.set_xticks([])
    ax.set_yticks([])


    # coordinates!

    # ploting T2 coordinates
    for i in range(1,int(1/(2*alpha))):

        y=i*alpha

        x2 = np.linspace(0, T2_lim(y)[1], 100)
        ax.plot(x2, T2(y,x2), "dodgerblue", linewidth=1)

    # ploting T1 & T3 coordinates
    for i in range(1, int(1/alpha)):

        y=i*alpha

        # ploting T1 coordinates
        plt.hlines(y = y*h, xmin = 0, xmax = T1_lim(y)[1],color="crimson",linewidth=1)


        # ploting T3 coordinates
        x3 = np.linspace(T3_lim_symm(y)[0], T3_lim_symm(y)[1], 100)
        ax.plot(x3, T3(y,x3),"gold", linewidth=1)

    # single vline with full ymin and ymax
    plt.vlines(x=0, ymin=0, ymax=h, colors="grey", ls=':')

    #p_plt=[] # a list of the mid -points of the subtraingles, to plot the p-values via scatter function

    for i in range(res["D-LR"].size): # going through all the right sub-triangles

        a1=res["coord. (T1, T2, T3)"][i][0][0]
        b1=res["coord. (T1, T2, T3)"][i][0][1]

        a2=res["coord. (T1, T2, T3)"][i][1][0]
        b2=res["coord. (T1, T2, T3)"][i][1][1]

        a3=res["coord. (T1, T2, T3)"][i][2][0]
        b3=res["coord. (T1, T2, T3)"][i][2][1]

        # sketching the subtriangle
        trianglex, triangley, direction=return_triangle_coord(a1,b1,a2,b2,a3,b3)

        if np.isnan(res["D-LR"][i]): # this means this triangle is empty and so is his reflection
            plt.fill(trianglex, triangley,color="black")

        else:  
            # Normalizing to obtain a score within the [0, 1] range for color-coding based on the d_lr result.
            d_lr_color_score = (res["D-LR"][i] + 1 )/2 
            plt.fill(trianglex, triangley,color=(d_lr_color_score, 1-d_lr_color_score, 1-d_lr_color_score))

            # Plotting points that represent statistically significant p-values within the relevant subtriangles.
            x,y = mid_point_triangle(a1,b1,a2,b2,a3,b3)
            p=res["p-value(g-test)"][i]
            if p < 0.05 and p >= 0.001 :
                ax.scatter(x,y,color='yellow', marker="*", alpha = 0.4, s =9)
            if p < 0.001 and p>= 10**(-5):
                ax.scatter(x,y,color="darkslateblue", marker="*", alpha =0.9, s =22)
            if p<= 10**(-5):
                ax.scatter(x,y,color="black", marker="*", alpha = 1, s =25)



    # creating the legends
    
    # colors scheme of the sub-triangles
    d_lr_color_score = 0
    pp1 = plt.Rectangle((0.2, 0.85), 0.05, 0.05, color=(d_lr_color_score, 1-d_lr_color_score, 1-d_lr_color_score)) # d_lr_color_score = 0
    d_lr_color_score = 0.5
    pp2 = plt.Rectangle((0.25, 0.85), 0.05, 0.05, color=(d_lr_color_score, 1-d_lr_color_score, 1-d_lr_color_score)) # d_lr_color_score = 0.5
    d_lr_color_score = 1
    pp3 = plt.Rectangle((0.3, 0.85), 0.05, 0.05, color=(d_lr_color_score, 1-d_lr_color_score, 1-d_lr_color_score))# d_lr_color_score = 1
    pp4 = plt.Rectangle((0.518, 0.85), 0.05, 0.05, color="black")# empty triangle

    ax.add_patch(pp1)
    ax.add_patch(pp2)
    ax.add_patch(pp3)
    ax.add_patch(pp4)
    plt.text(0.16, 0.865, "$D_{lr} =      -1$", size=8)
    plt.text(0.26, 0.865, '0', size=8)
    plt.text(0.31, 0.865, '1', size=8)
    plt.text(0.38, 0.865, 'empty triangle', size=8)

    # the p-values legend
    #p <= 10**-5
    ax.scatter(0.2, 0.8,color="black", marker="*", alpha = 1, s =25)
    plt.text(0.214, 0.8, "$p < 10^{-5}$", size=8)
    # p <= 0.001
    ax.scatter(0.2, 0.77,color="darkslateblue", marker="*", alpha =0.9, s =22)
    plt.text(0.214, 0.77, "$p < 0.001$", size=8)
    # p_ex <= 0.05
    ax.scatter(0.2, 0.74,color="darkgoldenrod", marker="*", alpha = 0.4, s =9)
    plt.text(0.214, 0.74, "$p < 0.05$", size=8)


    # dubbing the nodes
    plt.text(-0.02, 0.88, 'T1', size=12, color="crimson")
    plt.text(-0.03, -0.01, 'T2', size=12, color = "dodgerblue")
    plt.text(0.54, -0.01, 'T3', size=12, color= "darkgoldenrod")


    # printing the coordinates

    coord= np.arange(0,1+alpha,alpha)
    T_1 = np.arange(0,0.5+alpha/2,alpha/2)
    T_2_3 = np.arange(0,0.5+alpha,alpha)

    # coordintes T1
    for i in range(len(T_1)):
        label= str(round(1-coord[i],2))
        x=T_1[i]
        y=T2(0,x)
        plt.text(x+0.01, y, label , size=7, color = "firebrick")

    # coordintes T2
    for i in range(len(T_2_3)-2):
        label= str(round(coord[i]+alpha,2))
        x=0
        y=T2(coord[i]+alpha,0)
        plt.text(x-0.033, y, label , size=7, color = "dodgerblue")  

    # coordintes T3    
    for i in range(len(T_2_3)):
        label= str(round(coord[i]+0.5,2))
        x=T_2_3[i]
        plt.text(x, -0.03, label , size=7, color = "darkgoldenrod")   

    # removing the box lines around the plot
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['top'].set_color('none')

    title = file_name+"_analysis_granuality_" + str(alpha)  +".png"  
    plt.savefig(title)
#     return fig


# In[10]:


# Plotting the fundamental two main subtriangles for symmetry analysis,
# saving the analysis plot, and returning the results (d_lr, g_test, _p_value) of the comparison.

def plot_fundemental_asymmetry(data, file_name):
    fig = plt.figure(figsize=(6, 4))
    ax = plt.axes()

    x_side_T2 = np.linspace(0, 0.5, 100)
    x_side_T3 = np.linspace(-0.5, 0, 100)

    ax.plot(x_side_T2, T2(0,x_side_T2),"k",linewidth=1);
    ax.plot(x_side_T3, T3(0,x_side_T3),"k",linewidth=1);
    plt.hlines(y = 0, xmin =-0.5, xmax = 0.5,color="k",linewidth=1)
    
    # Hide X and Y axes tick marks
    ax.set_xticks([])
    ax.set_yticks([])

    # single vline with full ymin and ymax
    plt.vlines(x=0, ymin=0, ymax=h, colors="black")


    # sketching the subtriangles
    # coordinates of right-main subtriangle
    trianglex_R=[0,0,0.5,0] 
    triangley_R = [0,h,0,0]

    # coordinates of right-main subtriangle
    trianglex_L=[-0.5,0,0,-0.5] 
    triangley_L = [0,h,0,0]

    main_n_r, main_n_l, main_d_lr,main_g_test, main_p_value = fundemental_asymmetry(data)

    main_d_lr_left = (main_n_l-0.5*(main_n_l+main_n_r) )/(0.5*(main_n_l+main_n_r)) # for the plotting

    d_lr_color_score_R = (main_d_lr + 1 )/2 # normalzing to get a score between [0,1]
    d_lr_color_score_L = (main_d_lr_left + 1 )/2 # normalzing to get a score between [0,1]

    # Filling the subtriangles color according to the d_lr result
    plt.fill(trianglex_R, triangley_R,color=(d_lr_color_score_R, 1-d_lr_color_score_R, 1-d_lr_color_score_R))
    plt.fill(trianglex_L, triangley_L,color=(d_lr_color_score_L, 1-d_lr_color_score_L, 1-d_lr_color_score_L))

    # sketching points that indicate significant p-values in the right subtraingle
    x= 0.15
    y= 0.4*h
    p=main_p_value
    if p < 0.05 and p >= 0.001 :
        ax.scatter(x,y,color='yellow', marker="*", alpha = 0.4, s =9)
    if p < 0.001 and p>= 10**(-5):
        ax.scatter(x,y,color="darkslateblue", marker="*", alpha =0.9, s =22)
    if p<= 10**(-5):
        ax.scatter(x,y,color="black", marker="*", alpha = 1, s =25)


    # creating the legends
    # colors of the sub-triangles
    d_lr_color_score = 0
    pp1 = plt.Rectangle((0.2, 0.85), 0.05, 0.05, color=(d_lr_color_score, 1-d_lr_color_score, 1-d_lr_color_score)) # d_lr_color_score = 0
    d_lr_color_score = 0.5
    pp2 = plt.Rectangle((0.25, 0.85), 0.05, 0.05, color=(d_lr_color_score, 1-d_lr_color_score, 1-d_lr_color_score)) # d_lr_color_score = 0.5
    d_lr_color_score = 1
    pp3 = plt.Rectangle((0.3, 0.85), 0.05, 0.05, color=(d_lr_color_score, 1-d_lr_color_score, 1-d_lr_color_score))# d_lr_color_score = 1
    ax.add_patch(pp1)
    ax.add_patch(pp2)
    ax.add_patch(pp3)
    plt.text(0.135, 0.865, "$D_{lr} =-1$", size=7)
    plt.text(0.26, 0.865, '0', size=7)
    plt.text(0.31, 0.865, '1', size=7)
    
    # the p-values legend
    #p <= 10**-5
    ax.scatter(0.2, 0.8,color="black", marker="*", alpha = 1, s =25)
    plt.text(0.214, 0.8, "$p < 10^{-5}$", size=8)
    # p <= 0.001
    ax.scatter(0.2, 0.77,color="darkslateblue", marker="*", alpha =0.9, s =22)
    plt.text(0.214, 0.77, "$p < 0.001$", size=8)
    # p_ex <= 0.05
    ax.scatter(0.2, 0.74,color="darkgoldenrod", marker="*", alpha = 0.4, s =9)
    plt.text(0.214, 0.74, "$p < 0.05$", size=8)


    # adding text
    title_left = str(main_n_l)
    title_right = str(main_n_r)
    plt.text(-0.5, -0.1," n =" , size=12)
    plt.text(-0.3, -0.1, title_left, size=12, color = "grey")
    plt.text(0.2, -0.1, title_right, size=12, color= "grey")


    # removing the box lines around the plot
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['top'].set_color('none')

    #saving figure
    title = file_name + "_fundamental_asymmetry.png"  
    plt.savefig(title)
    return (main_d_lr,main_g_test, main_p_value)


# In[11]:


# Generating an index plot with distinct numbers assigned to each subtriangle for clear visualization.
# These numbers will correspond to the indices in the result file.

def plotting_triangle_index(res, granuality):
    # Initialize default values
    alpha = granuality
    fig_size = (5, 4)  # default size
    font_size = 8      # default font size

    # Set specific values based on granularity
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

    ax.plot(x_side_T2, T2(0,x_side_T2),"k",linewidth=1);

    plt.hlines(y = 0, xmin =0, xmax = 0.5,color="k",linewidth=1)
    # Hide X and Y axes tick marks
    ax.set_xticks([])
    ax.set_yticks([])

    # coordinates!

    # ploting T2 coordinates
    for i in range(1,int(1/(2*alpha))):
        y=i*alpha
        x2 = np.linspace(0, T2_lim(y)[1], 100)
        ax.plot(x2, T2(y,x2), "dodgerblue", linewidth=1)

    # ploting T1 & T3 coordinates
    for i in range(1, int(1/alpha)):
        y=i*alpha

        # ploting T1 coordinates
        plt.hlines(y = y*h, xmin = 0, xmax = T1_lim(y)[1],color="crimson",linewidth=1)

        # ploting T3 coordinates
        x3 = np.linspace(T3_lim_symm(y)[0], T3_lim_symm(y)[1], 100)
        ax.plot(x3, T3(y,x3),"gold", linewidth=1)

    # single vline with full ymin and ymax
    plt.vlines(x=0, ymin=0, ymax=h, colors="grey", ls=':')

    N=number_triangles(alpha)
    for i in range(res["D-LR"].size): # going through all right sub-triangles
        a1=res["coord. (T1, T2, T3)"][i][0][0]
        b1=res["coord. (T1, T2, T3)"][i][0][1]

        a2=res["coord. (T1, T2, T3)"][i][1][0]
        b2=res["coord. (T1, T2, T3)"][i][1][1]

        a3=res["coord. (T1, T2, T3)"][i][2][0]
        b3=res["coord. (T1, T2, T3)"][i][2][1]

        # sketching points that indicate significant p-values in the relevent subtriangles
        x,y = mid_point_triangle(a1,b1,a2,b2,a3,b3)
        index= res["index"][i]
        
        plt.text(x-0.01, y, str(index), size=font_size) # different font sizes for the different granulaties

    # dubbing the nodes
    plt.text(-0.02, 0.88, 'T1', size=12, color="crimson")
    plt.text(-0.03, -0.01, 'T2', size=12, color = "dodgerblue")
    plt.text(0.54, -0.01, 'T3', size=12, color= "darkgoldenrod")

    # printing the coordinates
    coord= np.arange(0,1+alpha,alpha)
    T_1 = np.arange(0,0.5+alpha/2,alpha/2)
    T_2_3 = np.arange(0,0.5+alpha,alpha)

    # coordintes T1
    for i in range(len(T_1)):
        label= str(round(1-coord[i],2))
        x=T_1[i]
        y=T2(0,x)
        plt.text(x+0.01, y, label , size=7, color = "firebrick")

    # coordintes T2
    for i in range(len(T_2_3)-2):
        label= str(round(coord[i]+alpha,2))
        x=0
        y=T2(coord[i]+alpha,0)
        plt.text(x-0.033, y, label , size=7, color = "dodgerblue")  

    # coordintes T3    
    for i in range(len(T_2_3)):
        label= str(round(coord[i]+0.5,2))
        x=T_2_3[i]
        plt.text(x, -0.03, label , size=7, color = "darkgoldenrod")   

    # removing the box lines around the plot
    ax.spines['right'].set_color('none')
    ax.spines['left'].set_color('none')
    ax.spines['bottom'].set_color('none')
    ax.spines['top'].set_color('none')

    # saving figure
    title = "index_granulality_" + str(alpha)  +".png"  
    plt.savefig(title)
    return fig


# # Analysis- help functions

# In[12]:


# Implementing basic symmetry analysis between the two primary subtriangles, i.e., those divided by the y-axis.
# A new column containing the x-axis values of each point is created.
# A point is considered to be in the right subtriangle if and only if its x-coordinate is strictly greater than 0.
# (Note that points on the y-axis were eliminated during preprocessing.)

def fundemental_asymmetry(data):
    data["x-axis"]= cartizian(data["T1"],data["T2"],data["T3"])[0] 
    main_n_r = len(data[data["x-axis"]>0])
    main_n_l = len(data[data["x-axis"]<0])

    main_d_lr=D_LR(main_n_r,main_n_l)
    main_g_test, main_p_value= log_likelihood_ratio_test(main_n_r,main_n_l)
    return (main_n_r, main_n_l, main_d_lr,main_g_test, main_p_value)


# In[13]:


# Given coordinates for a valid subtriangle (either right or left), determine the number of data points within the
# triangle and its corresponding reflected triangle. Regardless of whether the coordinates were provided for the left
# or right subtriangle, it returns (n_r, n_l) for the respective directions.

#The input coordinates ai < bi signifies the interval point in the coordinates of T_i, i=1,2,3

def n(a1,b1,a2,b2,a3,b3, data): 
    
    # This additional logic is implemented to avoid double-counting of data. The coordinates divide the range [0, 1]
    # into a series of half-open intervals of the form (i * alpha, (i + 1) * alpha], 
    # but the very first interval is closed, covering the range [0, 1 * alpha].
    if a1 == 0:
        condition_a1 = a1<=data.T1
    else   : 
        condition_a1 = a1<data.T1 
    
    if a2 == 0: 
        condition_a2 = a2<=data.T2
    else    :
        condition_a2 = a2<data.T2 
        
    if a3 == 0: 
        condition_a3 = a3<=data.T3
    else    :
        condition_a3 = a3<data.T3     
        
    n_1=    len(data[ ((a1<=data.T1) & (data.T1<=b1)) &    # the count in the given subtriangle 
                ((a2<=data.T2) & (data.T2<=b2))&
                ((a3<=data.T3) & (data.T3<=b3))])
    

    
    n_2=    len(data[ ((a1<=data.T1) & (data.T1<=b1)) &  # the count in the reflected subtriangle 
            ((a3<=data.T2) & (data.T2<=b3))&
            ((a2<=data.T3) & (data.T3<=b2))])
    
    # checking which triangle is left in which is right
    
    # The check involves determining whether the x-axis component of the top node of a subtriangle is positive or negative
    # in Cartesian coordinates. The top node is positive along the x-axis if and only if it belongs to the right-hand side triangle.

    trianglex, triangley, direction = return_triangle_coord(a1,b1,a2,b2,a3,b3)
    top_node = trianglex[1] # this is the x_axis coordinate of the top node in the given triangle 
    
    if top_node > 0: # the coordinates we were given were of a right-side triangle
        n_r =  n_1
        n_l =  n_2
    else:            # the coordinates we were given were of a left-side triangle
        n_r =  n_2
        n_l =  n_1
    
    
    return (n_r,n_l) 



# the D_LR comparison 
# Calculate the d_lr value, where n_l represents the number of points in a left triangle, 
# and n_r represents the number of points in the corresponding right triangle.

def D_LR(n_r,n_l):
    
    if n_l+n_r != 0 : # finally the test
        d_lr = (n_r-0.5*(n_l+n_r))/(0.5*(n_l+n_r)) 
    else: # the case we would have divided by 0
        d_lr = np.nan # in this case we return a NaN value
        
    return d_lr


# In[14]:


# Determine the required number of triangles to compare for a given granularity.
# We define a "row" as a coordinate line of T1, where the bottom line is represented by T1(0),
# and each subsequent row is incremented by T1(1 * alpha) until we reach T1(1 / alpha - 1) * alpha.

# Our choice of alpha is designed to ensure that the y-axis, passing through 0, does not intersect any triangle.
# Specifically, we select an alpha value such that 1/alpha is an even number.
# Consequently, in the bottom row, the count of "up-triangles" (where the base is parallel to the x-axis and facing downward)
# is equal to 0.5 * (1/alpha).
# We then systematically traverse row by row, tallying both "up triangles" and "down triangles."


# Given a granularity level, we determine the number of triangles that require examination.
def number_triangles(alpha): 
    
    if  int(1/alpha) % 2:
        print(" 1/alpha is odd")
        sys.exit() 
    
    # Explanation:
    # In the first row, there are 1/(2 * alpha) up triangles, and then for each subsequent row, we have two fewer up triangles,
    # where the count decreases according to the expression (1/(2 * alpha) - i), for i = 1 until (1/(2 * alpha) - i) == 1.
    # To find the value of i where (1/(2 * alpha) - i) == 1, we solve for i using the equation i == (1 - 2 * alpha) / (2 * alpha),
    # which can be represented as 'a'.
  
    a = int((1-2*alpha)/(2*alpha))
    n = int(1/(2*alpha)) # n=# triangles. initiation is of  up-triangles in row 0 - bottom of the triangle
    for i in range(1,a+1):
        # Up to the top of the triangles, there are double rows of "diamonds," which we define as a combination of a bottom triangle
        # and an up triangle sharing the same base. Since we have double rows in total, we multiply this count by 4.
        n = n + 4*int((1/(2*alpha))-i) 

    return n


# In[15]:


# returns -2* log(likelihood(n_l))/log(likelihood(symmetric)), under the binomial
def log_likelihood_ratio_test(n_r,n_l): 
   N= n_r +n_l
   if N != 0: # the non-pathalogical case        
      
       # Calculating the likelihood of landing 'n_l' times in the left triangle out of 'N' trials, each with a probability of 0.5
       # to land in either 'n_l' or 'n_r'.
       L_null = binom.pmf(n_l,N , 0.5) 
       
       # The alternative hypothesis posits that the actual proportion of 'n_l' occurrences is precisely
       # equal to the observed proportion in the experiment.
       L_alt = binom.pmf(n_l,N , n_l/N) 

       if L_null/L_alt < 10**-16: # this means sth extremely unlikely has occured
           test_res= np.nan
           p_value = 0  
       else:    
           test_res = float(-2 * log(L_null/L_alt))
           p_value = 1 - chi2.cdf(test_res, 1) # degress of freedom = 1, according to https://stats.libretexts.org/Bookshelves/Applied_Statistics/Biological_Statistics_(McDonald)/02%3A_Tests_for_Nominal_Variables/2.04%3A_GTest_of_Goodness-of-Fit
   else: # the case where n_l=n_r =0
           test_res= np.nan
           p_value = np.nan  

   return (test_res,p_value)


# In[16]:


# To perform the comparison, we iterate row by row, following the T1(alpha,.) coordinates,
# and tally the count of up-triangles followed by down-triangles in each row.

def triangles_analysis(data, granularity, file_name):
    # Convert granularity to alpha value
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
        # If granularity is already a number, use it directly
        alpha = float(granularity)

    plot(data, alpha, file_name)  
    tri = [] # temporary list storing the values

    for i in range(int(1/alpha)): # we go row by row, T1 coordinates
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

            trianglex, triangley, direction = return_triangle_coord(a1,b1,a2,b2,a3,b3)
            if round(trianglex[0],4) < 0:
                continue

            n_r, n_l = n(a1,b1,a2,b2,a3,b3, data)
            d_lr = D_LR(n_r,n_l)
            g_test, p_value = log_likelihood_ratio_test(n_r,n_l)

            coord = [(round(a1,4),round(b1,4)),(round(a2,4),round(b2,4)),(round(a3,4),round(b3,4))]
            tri.append([coord,n_r,n_l, d_lr, g_test, p_value])

            k3 = k3 + 1 # going from an 'up' triangle to a 'down' triangle, from right to left 

            a2 = k2 * alpha
            b2 = (k2+1) * alpha

            a3 = 1-(k3+1) * alpha
            b3 = 1 - k3 * alpha

            trianglex, triangley, direction = return_triangle_coord(a1,b1,a2,b2,a3,b3)
            if round(trianglex[0],4) < 0:
                continue

            n_r, n_l = n(a1,b1,a2,b2,a3,b3, data)
            d_lr = D_LR(n_r,n_l)
            g_test, p_value = log_likelihood_ratio_test(n_r,n_l)

            coord = [(round(a1,4),round(b1,4)),(round(a2,4),round(b2,4)),(round(a3,4),round(b3,4))]
            tri.append([coord,n_r,n_l, d_lr, g_test, p_value])

            k2 = k2 + 1

    # saving the list as a dataframe
    triangles = pd.DataFrame(tri)
    triangles.columns = ["coord. (T1, T2, T3)","n_right","n_left","D-LR", "g-test", "p-value(g-test)"]
    index = list(range(number_triangles(alpha), 0, -1)) # indexing the triangles
    triangles["index"] = index # the index is fitting the index-plot
    
    return triangles


# # Data- processing

# In[17]:


# Accept a CSV data file as input, assuming that the first column represents T1 coordinates, 
# the second column represents T2 coordinates, and the third column represents T3 coordinates. 
# We filter out points located on the y-axis, as they may contribute to both the left and right subtriangles.

def dump_data(file):
    import csv
    import warnings
    
    # Suppress the FutureWarning about delim_whitespace
    warnings.filterwarnings('ignore', category=FutureWarning)
    
    print(f"\nAttempting to read file: {file}")
    
    # First try to read the file and detect its format
    with open(file, 'r', encoding='utf-8') as f:
        # Read a sample to detect the format
        sample = f.read(1024)
        f.seek(0)  # Go back to start
        
        # Try to detect the dialect
        try:
            dialect = csv.Sniffer().sniff(sample)
            print(f"Detected format: delimiter='{dialect.delimiter}', quoting={dialect.quoting}")
        except:
            print("Could not detect format automatically, will try multiple formats")
            dialect = None
        
        # Read all lines
        lines = f.readlines()
    
    print("\nFirst few lines of the file:")
    for i, line in enumerate(lines[:5]):
        print(f"Line {i}: {line.strip()}")
    
    # Try multiple approaches to read the data
    data = None
    delimiters = [',', '\t', ';', ' ', '|']  # Common delimiters
    
    # First try with detected dialect if available
    if dialect:
        try:
            data = pd.read_csv(file, dialect=dialect, header=None, names=['T1', 'T2', 'T3'])
            print(f"Successfully read with detected dialect: {dialect.delimiter}")
        except Exception as e:
            print(f"Failed to read with detected dialect: {str(e)}")
    
    # If that didn't work, try each delimiter
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
    
    # If still no success, try reading line by line
    if data is None or data.shape[1] != 3:
        print("\nTrying to read line by line...")
        rows = []
        for line in lines:
            line = line.strip()
            if not line:
                continue
                
            # Try different delimiters
            for sep in delimiters:
                cols = [col.strip() for col in line.split(sep)]
                if len(cols) == 3:
                    try:
                        # Try to convert to numbers
                        row = [float(col) for col in cols]
                        rows.append(row)
                        break
                    except ValueError:
                        continue
        
        if rows:
            data = pd.DataFrame(rows, columns=['T1', 'T2', 'T3'])
            print("Successfully read data line by line")
    
    # Verify we have data
    if data is None or data.shape[1] != 3:
        raise ValueError(f"Could not read file with 3 columns. Please check the file format.")
    
    print(f"\nSuccessfully read {len(data)} rows of data")
    
    # Convert all columns to numeric, coercing errors to NaN
    for col in data.columns:
        data[col] = pd.to_numeric(data[col], errors='coerce')
    
    # Drop any rows with NaN values
    data = data.dropna()
    
    if len(data) == 0:
        raise ValueError("No valid numeric data found after processing")

    # We normalize the rows, as a precaution, in case they were not already normalized.
    n_rows = data.shape[0]
    for i in range(n_rows):
        s = sum(data.iloc[i,:])
        if s == 0:
            continue  # Skip rows that sum to zero
        data.iloc[i,:] = (data.iloc[i,:])/s
    
    data = data.loc[data.iloc[:,1] != data.iloc[:,2]] 
   
    return data


# # User Interface

# In[19]:


# The main function- runs the anlysis
# 
# Input
# 1. file: Ensure the input file is a CSV file containing precisely three columns. Assume the first column represents T1, the second represents T2, and the third represents T3.
# 2. granuality: Specify the desired granularity by choosing one of the following options: "coarse", "fine", or "superfine". Keep in mind that "coarse" corresponds to a granularity of 0.25, "fine" to 0.1, and "superfine" to 0.05.

# Output : The program will generate a Results file comprising the following components:
# 1. CSV Results Table: This file will contain the results in CSV format.
# 2. Subtriangles Index Plot: An index plot with distinct numbers assigned to each subtriangle for clear visualization. These numbers correspond to the indices in the result file.
# 3. Ternary Coordinates Plot: A plot of the data in ternary coordinates.
# 4. Results Summary Plot: A plot displaying the results with D_lr results together with significant G-test scoring.
# 5. Basic Asymmetry Plot: A plot showcasing the basic asymmetry between the left and right halves of the triangle.


def run_analysis(file,granuality): 
    try:
        # creating and redirecting to a Results folder
        original_path=os.getcwd() #current path
        desired_directory=Path(original_path+"/Results") #creating a desired directory, if it doesnt yet exist
        if not desired_directory.exists(): #so we avoid an error of trying to create it more than once
            # this is the directory that will be created
            path = os.getcwd() +"/Results" # Gets your current directory and make a path to a folder for the results
            os.mkdir(path) # creating the directory (realizing the path)
        os.chdir(desired_directory) # we move to the new directory

        # data upload/initial processing
        data= dump_data(original_path+"/"+file) # data processing, data is at precious directory
        
        # to deal with pathological cases of an empty data-set
        if data.size == 0:
            print("Dataset is incompatible")
            print("Perform a validation check to determine whether the uploaded data is either empty or "+ 
                   "exclusively consists of points located on the reflection axis")
            sys.exit()
            
        # Get just the filename without path and extension
        file_name = os.path.splitext(os.path.basename(file))[0]
        
        #plotting data
        basic_plot=plot_fundemental_asymmetry(data,file_name)
        
        # analysis
        result = triangles_analysis(data, granuality,file_name) 
        
        # plotting
        fig_index=plotting_triangle_index(result, granuality);# plotting the index and saving the figure
        fig_result=plot_results(result,granuality,file_name) # plotting results and saving the figure
        
        # turning the coordinate trinagle plotting index to the actual index of the dataframe
        result.set_index("index", inplace=True)
        
        # adding a row with the information for the fundemental asymmetry
        (main_n_r, main_n_l, main_d_lr,main_g_test, main_p_value)=fundemental_asymmetry(data)
        result.loc['full dataset'] = ['NA',main_n_r,main_n_l,main_d_lr,main_g_test,main_p_value]
        
        result =result.iloc[::-1] # reversing order so the first row is index 1...
        
        # saving as csv
        title_csv = file_name +"_results_granuality_" + str(granuality)  +".csv"  
        result.to_csv(title_csv, index=True) # saving results dataFrame as a csv-file
        
        os.chdir(original_path) # return to the original directory
        return result
        
    except Exception as e:
        print(f"\nError processing data: {str(e)}")
        print("Please make sure the file is a valid CSV with three columns (T1, T2, T3)")
        if os.getcwd() != original_path:
            os.chdir(original_path)  # Make sure we return to the original directory even if there's an error
        raise  # Re-raise the exception after handling it


# In[20]:


# Conducts an initial assessment of basic asymmetry between the two primary subtriangles, divided by the y-axis.
# Returns the counts of data points in each triangle, the d_lr result, the G-test statistic, and its associated p-value.
# Additionally, it saves a figure illustrating the results.

# Input
# file: Ensure the input file is a CSV file containing precisely three columns. Assume the first column represents T1, the second represents T2, and the third represents T3.


def run_basic_analysis(file):    
    data= dump_data(file) # data processing
    if data.size == 0:
        print("Dataset is incompatible")
        print("Perform a validation check to determine whether the uploaded data is either empty or "+ 
               "exclusively consists of points located on the reflection axis")
        sys.exit()
        
    (main_n_r, main_n_l, main_d_lr,main_g_test, main_p_value)=fundemental_asymmetry(data)
    
    file_name = file[:len(file) - 4] # getting rid of ".csv" in the file's name
    basic_plot=plot_fundemental_asymmetry(data, file_name)
    print("d_lr =", main_d_lr)
    print("\ng-test =", main_g_test)
    print("\ng-test- p value=", main_p_value)
    return (main_n_r, main_n_l, main_d_lr,main_g_test, main_p_value)










