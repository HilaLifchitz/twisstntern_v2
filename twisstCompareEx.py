#!/usr/bin/env python
# coding: utf-8

# # COMPARING METRICS FOR COMPARING TERNARY DATA

# In[23]:


# importing libraries
from pathlib import Path
import ternary

import numpy as np
import random
import math
import scipy
import seaborn as sns
import math


import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import ListedColormap, Normalize
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#from scipy.stats import chisquare
from scipy.stats import chi2


from scipy.stats import wasserstein_distance

# presentation in the dataframe
pd.set_option('display.float_format', '{:.4e}'.format)
# activating latex printing in matplotlib
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "Helvetica"
})

# importing the twisstntern library
import twisstntern

# importing the optimal transport library
import ot
import os


# ## Data Processing
# ### **NOTICE**  im sampling "only" AmountSamplingPoints (~ 10,000) points from the data set to make computations faster

# In[24]:


################################################################## DETERMINE AMOOUNT OF SAMPLING POINTS##########################
AmountSamplingPoints= 13000;
########################################################################################################################


# In[49]:


# Principally we expect the data that we compare to already have been normalized (output of twisstentern)
# making sure the data file fits- a specific function for this data set  
# OUTPUT: an array



# def data_process(file_name):
#     original_dir = os.getcwd()  # Save the current working directory
#     try:
#         os.chdir('gene_flow_sims') # going to where Sean's data is stored
#         data = pd.read_csv(file_name,sep= ';|\t|,', header=None,engine='python')
#         data=data.iloc[4:]
#         data=data.iloc[:, [0,1,2]]
#         data = data.reset_index(drop=True)

#         if file_name == "G_m0.04.csv":
#             data = data.drop(data.index[np.where(data.values == 'topo1')[0][0]])

#         data=data.astype(int) #for some reason it reads the data.values as strings, this command fixes this

#         data = data.sample(n=AmountSamplingPoints) ########################## HERE IS THE SAMPLING!!!
#         data_array = np.divide(data.values,data.values.sum(1, keepdims=True)) #normalizing

#         data_array=data_array.astype(np.float32)

#     finally:
#         os.chdir(original_dir)  # Always go back, even if error occur
        
#     return data_array    

# #data Array to data in dataframe
# def dumpData(dataArray):
#     data = pd.DataFrame(dataArray, columns=['T1', 'T2', 'T3'])
#     n_rows = data.shape[0]

#     # Normalize rows
#     for i in range(n_rows):
#         s = sum(data.iloc[i, :])
#         data.iloc[i, :] = (data.iloc[i, :]) / s

#     data = data.loc[data.iloc[:, 1] != data.iloc[:, 2]]
#     return data


# # Metrics

# ## Wassersteins (Regular Euclidean  and the KL tweaked version)

# In[26]:


def wasserstein_distance(data1, data2):
    # In the ususal wasserstein the cost matrix: Euclidean distance between all points, even as they are in the simplex
    
    # Uniform weights for each point
    a = np.ones((data1.shape[0],)) / data1.shape[0]
    b = np.ones((data2.shape[0],)) / data2.shape[0]
    
    #Cost matrix: Euclidean distance between all point
    M = ot.dist(data1, data2)
    
    # Compute the (regular- euclidean) Wasserstein distance
    dist = ot.emd2(a, b, M, numItermax=1000000)
    
    return dist 


# ### The cost function for the KL-Wasserstein:
# Let $x=(x_1,x_2,x_3) \in data_1, y=(y_1,y_2,y_3) \in data_2$ be ternary points in each respective data set. Both $x,y$ can be regarded as distributions(over $\{1,2,3\}$).<br>We define $Cost(x,y):= KL(x||y) + KL(y||x):= \sum_{k=1}^3 x_k log(\frac{x_k}{y_k}) + y_k log(\frac{y_k}{x_k})$

# In[27]:


# data1, data2 are pandas datasets, containing lists of ternary points whocse two distributions we are comparing 
# each row in data_i is a point on the simplex (the three columns of the columns are alwats T1,T2 and T3 according to the
# convention T1 on top, T2 on the right and T3 on the left of the triangle.) 
# Hence data1.shape=(n1,3),data2.shape=(n2,3), its allowed that n1 != n2. 

def simplex_wasserstein_distance(data1,data2):
    
    # adding an epsilon value to the data set to avoid division by zero in the KL calculation
    data1_plus_eps = data1 + 1e-8
    data1_plus_eps = data1_plus_eps / data1_plus_eps.sum(1, keepdims=True) #normalizing

    data2_plus_eps = data2 + 1e-8
    data2_plus_eps = data2_plus_eps / data2_plus_eps.sum(1, keepdims=True) #normalizing


    # calculating the symmetric KL distance between each data point in data1 p=(p1,p2,p3) to data points in data2 q=(q1,q2,q3)
    # according to the kl divergence formula:sum_{i=1 to 3}p_i log(p_i\q_i)
    #
    kl_pq = (data1_plus_eps[None, :, :] * (np.log(data1_plus_eps[None, :, :]) - np.log(data2_plus_eps[:, None, :]))).sum(-1)
    kl_qp = (data2_plus_eps[None, :, :] * (np.log(data2_plus_eps[None, :, :]) - np.log(data1_plus_eps[:, None, :]))).sum(-1)
    kl_symm = (kl_pq + kl_qp.T)

    # now we call the optimal transport package to solve the cost matrix
    dist = ot.emd2(np.ones(data2_plus_eps.shape[0])/data2_plus_eps.shape[0], 
            np.ones(data1_plus_eps.shape[0])/data1_plus_eps.shape[0],
            M=kl_symm, numItermax=int(3e6))
    
    #kl_symm=kl_symm.T # now the dimensions kl_symm.shape = (n1,n2)

    # now we call the optimal transport package to solve the cost matrix
 #   dist = ot.emd2(np.ones(data1_plus_eps.shape[0])/data1_plus_eps.shape[0], 
  #          np.ones(data2_plus_eps.shape[0])/data2_plus_eps.shape[0],
   #         M=kl_symm, numItermax=int(3e6))
    
    return dist


# ## $L^2$ and $\chi^2$ comparisons
# ### Auxilary functions for counting points in subtraingles

# In[28]:


#auxilary function to ocunt the number of data point in a sub triangle of *data* with the coordinates:
# T1:(a1,b1), T2:(a2,b2) T3:(a3,b3)
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
        
    n=    len(data[ (condition_a1 & (data.T1<=b1)) &    # the count in the given subtriangle 
                 (condition_a2 & (data.T2<=b2))&
                (condition_a3 & (data.T3<=b3))])
    return n



# taking a dataframe of the Data, and its granuality level and returning 
#a list of the subtriangles and how many data points are there in each such subtriangle
def SubtrianglesDataCount(data,alpha):
    triangles=[]
    steps = int(1 / alpha)
    for k in range(steps):
        #print("k = ", k)
        a1 = round(k * alpha, 10) # we use the round because of the quirk of representing alpha intervals in binary in a computer
        b1 = round((k + 1) * alpha,10)
        
        # T2 goes from 0 to (1 - (k+1)*alpha) in steps of alpha
        T2_upper_limit = round(1-k * alpha, 10)
        T2_steps = round(T2_upper_limit / alpha) 

        
        # First triangle (T3 from [1 - (k+1)*α, 1 - k*α]), progressing from left to right
        a3_1 = round(1 - (k + 1)* alpha,10)
        b3_1 = round(1 - k * alpha,10)
        
        for T2_step in range(T2_steps):
         #   print("T2 step= ", T2_step)
            a2 = round(T2_step*alpha,10)
            b2 = round((T2_step+1)*alpha, 10)

        
            if a3_1 >= 0:
                triangles.append({
                    'T1': (a1, b1),
                    'T2': (a2, b2),
                    'T3': (a3_1, b3_1),
                    'n' : n(a1,b1,a2,b2,a3_1,b3_1,data)
                })
                #print(triangles[-1])
           
        
            # Second triangle- T3 is -alpha from the last coords. some k values we wont have this second triagle
            a3_2 = round(a3_1- alpha,10)
            b3_2 = round(b3_1- alpha,10)
        
            if a3_2 >= 0:
                triangles.append({
                    'T1': (a1, b1),
                    'T2': (a2, b2),
                    'T3': (a3_2, b3_2),
                    'n' : n(a1,b1,a2,b2,a3_2,b3_2,data)
        
                })
            #print(triangles[-1])
            
             # updating the T3 coordinates for the next T2 step, the new first T3 triangle has the first coordiantes of the last T3 triangle
            a3_1 = a3_2;
            b3_1 = b3_2;
    triangles = pd.DataFrame(triangles)
    return triangles    


# Reminder about the chi square test: we are comparing datasets to a "true" distribution (given by dataTruth) over the triangles.
# chi2 square is the measures the discrepancy of the observed frequencies (data2) from the expected
# frequencies (dataTruth counts) under the null hypothesis (that the data2 is a realization of the trute distribution==the dataTruth).
# 
# Let d:= number of histogram slots= number of triangles, then
# $X^2 = \sum_{k=1}^{d}\frac{(\Delta_{k}^{1}-\Delta_{k}^{2})^2}{\Delta_{k}^{1}}=\sum_{k=1}^{d}\frac{(\Delta_{k}^{2}-\mathbb{E}_{H_0}(\Delta_{k}))^2}{\mathbb{E}_{H_0}(\Delta_{k})} \sim \chi^2(d-1)$

# In[29]:


#INPUT: two data frames of the data in our usual form-- we assume theu have the same # points!!!! (at the moment only?)
#OUTPUT: returns a data fram with with the (L^2)^2 distance per subtriangle, calculates chi^2 test, L^2 distance between the datasets
#and a figure that visualizes both of them

def compare_triangle_counts(dataTruth, data2, alpha):
    # creating two dataframes with coordinates of the subtriangles and their respective data point counts
    triangles1= SubtrianglesDataCount(dataTruth,alpha) # we make the differentiation of what we consider truth as required for the chi^2 test
    triangles2= SubtrianglesDataCount(data2,alpha)
    
    result = triangles1.copy()
    # Rename 'n' to 'Truth'
    result = result.rename(columns={"n": "nTruth"})
    
    # Add column from triangles2
    result["n2"] = triangles2["n"]

    # Element-wise differnce dataTruth-data2
    result["data1-data2"] = (triangles1["n"] - triangles2["n"])
    
    # Element-wise difference of squares
    result["L2"] = (triangles1["n"] - triangles2["n"]) ** 2
    
    # Replace 0 with NaN to skip those rows
    safe_n = triangles1["n"].replace(0, np.nan)
    # Normalize by triangles1["n"] per row
    #result["L2_normalized"] = result["L2"] / triangles1["n"] -- but without dividing by 0
    result["L2_normalized"] = result["L2"] / safe_n
    
    # NOTICE: We cannot do Normalize by triangles1["n"] per row
    #(result["L2_normalized"] = result["L2"] / triangles1["n"])
    # because some of the triangles have 0 points in it, and we need to avoid dividing by 0.
    # this also turns the chi square to a bit less meaningfull btw
    
    # To get the this is d for the chi^2(d-1)-- we need to correct for the right amount of subtriangles investiated
    # so instead of the tottal amount of subtraingles:  (degreesFreedome= triangles1.shape[0] - 1)
    # we take the actuak amount without the 0th subtraingle of the in the dataTruth:
    # We cannot do Normalize by triangles1["n"] per row
    degreesFreedome= safe_n.notna().sum() -1
    
    
    #L2 distance = sqrt(\sum_{triangles} L^2 * (arean of triangle=1/(#triangles)))
    L2=math.sqrt(result["L2"].sum()/result.shape[0])
    chi_Statistic = result["L2_normalized"].sum()
    p_value = chi2.sf(chi_Statistic,degreesFreedome )
   
    return result, L2, chi_Statistic, p_value


# In[30]:


# os.chdir('gene_flow_sims')
# os.getcwd() #current path
# original_path=os.getcwd() #current path
# os.chdir('..')
# print(os.getcwd())


# # Plotting L2 per subtriangle

# In[31]:


# functions
# AUXILARY FUNCTION FOR THE GRADIENT AND THE HEAT MAP LEGEND
def gradient_ignite_log(x):
    x = min(1, max(1e-9, x))  # Clamp and avoid log(0)

    # Aggressive log transform to emphasize small values
    log_x = np.log10(x * 999 + 1) / np.log10(1000)  # maps [0,1] → [0,1]

    # Fire-inspired gradient: white → fiery orange → blood red → deep purple → black
    key_colors = {
        0.00: (1.00, 1.00, 1.00),  # White
        0.02: (1.00, 0.60, 0.10),  # Bright Orange
        0.08: (0.95, 0.20, 0.10),  # Blood Red
        0.25: (0.50, 0.00, 0.30),  # Wine Purple
        0.50: (0.20, 0.00, 0.30),  # Deep Purple
        1.00: (0.00, 0.00, 0.00),  # Black
    }

    keys = sorted(key_colors.keys())
    for i in range(len(keys) - 1):
        x0, x1 = keys[i], keys[i + 1]
        if x0 <= log_x <= x1:
            t = (log_x - x0) / (x1 - x0)
            c0 = key_colors[x0]
            c1 = key_colors[x1]
            r = c0[0] * (1 - t) + c1[0] * t
            g = c0[1] * (1 - t) + c1[1] * t
            b = c0[2] * (1 - t) + c1[2] * t
            return (r, g, b)

    return key_colors[1.0]  # fallback


# Define key color steps (for dramatic fiery visual)
key_colors = {
    0.00: (1.00, 1.00, 1.00),  # White
    0.02: (1.00, 0.60, 0.10),  # Bright Orange
    0.08: (0.95, 0.20, 0.10),  # Blood Red
    0.25: (0.50, 0.00, 0.30),  # Wine Purple
    0.50: (0.20, 0.00, 0.30),  # Deep Purple
    1.00: (0.00, 0.00, 0.00),  # Black
}

def interpolate_gradient(x):
    x = min(1, max(0, x))
    keys = sorted(key_colors.keys())
    for i in range(len(keys) - 1):
        x0, x1 = keys[i], keys[i + 1]
        if x0 <= x <= x1:
            t = (x - x0) / (x1 - x0)
            c0, c1 = key_colors[x0], key_colors[x1]
            r = c0[0] * (1 - t) + c1[0] * t
            g = c0[1] * (1 - t) + c1[1] * t
            b = c0[2] * (1 - t) + c1[2] * t
            return (r, g, b)
    return key_colors[1.0]

def make_plain_colormap(resolution=256):
    colors = [interpolate_gradient(x / (resolution - 1)) for x in range(resolution)]
    return ListedColormap(colors)

def add_plain_gradient_legend(ax, label="Proportion of total points"):
    cmap = make_plain_colormap()
    norm = plt.Normalize(vmin=0, vmax=1)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    # Position the inset colorbar
    cax = inset_axes(ax, width="3%", height="80%", loc='center right',
                     bbox_to_anchor=(0.05, 0, 1, 1), bbox_transform=ax.transAxes, borderpad=1)

    cbar = plt.colorbar(sm, cax=cax)
    cbar.set_label(label, fontsize=10)
    cbar.set_ticks([0.0, 0.02, 0.08, 0.25, 0.5, 1.0])
    cbar.set_ticklabels(["0.00", "0.02", "0.08", "0.25", "0.5", "1.0"])

###################################################################################################################################################
#### FUNCTION FOR GENERAL COUNT###########################################################################################
def PlotSubtrianglePositiveCount(data, alpha, ax=None):
    # Get triangle data
    triangles = SubtrianglesDataCount(data, alpha)

    # Create new axis if none is provided
    if ax is None:
        fig = plt.figure(figsize=(8, 6))
        ax = plt.axes()

    # coordinate skeleton    
    h = math.sqrt(3) / 2

    x_side_T2 = np.linspace(0, 0.5, 100)
    x_side_T3 = np.linspace(-0.5, 0, 100)

    ax.plot(x_side_T2, twisstntern.T2(0, x_side_T2), "k", linewidth=1)
    ax.plot(x_side_T3, twisstntern.T3(0, x_side_T3), "k", linewidth=1)
    ax.hlines(y=0, xmin=-0.5, xmax=0.5, color="k", linewidth=1)

    ax.set_xticks([])
    ax.set_yticks([])

    for i in range(1, int(1 / alpha)):
        y = i * alpha

        ax.hlines(y=y * h, xmin=twisstntern.T1_lim(y)[0], xmax=twisstntern.T1_lim(y)[1], color="#C74375", linewidth=1)

        x2 = np.linspace(twisstntern.T2_lim(y)[0], twisstntern.T2_lim(y)[1], 100)
        ax.plot(x2, twisstntern.T2(y, x2), color="#3182BD", linewidth=1)

        x3 = np.linspace(twisstntern.T3_lim(y)[0], twisstntern.T3_lim(y)[1], 100)
        ax.plot(x3, twisstntern.T3(y, x3), color="#D4A017", linewidth=1)

    ax.vlines(x=0, ymin=0, ymax=h, colors="#888888", ls=':')

    ax.text(-0.02, 0.88, 'T1', size=12, color="#C74375")
    ax.text(0.54, -0.01, 'T3', size=12, color="#D4A017")
    ax.text(-0.58, -0.01, 'T2', size=12, color="#3182BD")

    coord = np.arange(0, 1 + alpha, alpha)
    T_1 = np.arange(0, 0.5 + alpha / 2, alpha / 2)
    T_2 = np.arange(-0.5, 0 + alpha / 2, alpha / 2)
    T_3 = np.arange(-0.5, 0.5 + alpha, alpha)

    for i in range(len(T_1)):
        label = str(round(1 - coord[i], 2))
        x = T_1[i]
        y = twisstntern.T2(0, x)
        ax.text(x + 0.01, y, label, size=7, color="#C74375")

    for i in range(len(T_2)):
        label = str(round(1 - coord[i], 2))
        x = T_2[i]
        y = twisstntern.T3(0, x)
        ax.text(x - 0.04, y, label, size=7, color="#3182BD")

    for i in range(len(coord)):
        label = str(round(coord[i], 2))
        x = T_3[i]
        ax.text(x, -0.03, label, size=7, color="#D4A017")

         
    # filling the subtriangles with the counts
    numberPoints = triangles["n"].sum()
    numerTriangles = triangles.shape[0]
    expectedNormal = numerTriangles / numberPoints
    
    for index, row in triangles.iterrows():
        trianglex, triangley, direction = twisstntern.return_triangle_coord(
            row["T1"][0], row["T1"][1], row["T2"][0], row["T2"][1], row["T3"][0], row["T3"][1]
        )
        x = row["n"] / numberPoints
        color = gradient_ignite_log(x)
        ax.fill(trianglex, triangley, color=color)
    
    
    add_plain_gradient_legend(ax)
    return ax
####################################################################################################################################################################
# BASED ON THIS, WE ARE DOING A COMPARE TWO DATASETS PLOT FUNCTIONS
#############################################################################################################################
#INPUT SHOULD  the title of the plot, the table output of compare_triangle_counts(data1,data2,alpha), alpha

def PlotL2(title,results, alpha,ax=None):
     # Create new axis if none is provided
    if ax is None:
        fig = plt.figure(figsize=(8, 6))
        ax = plt.axes()

    # coordinate skeleton    
    h = math.sqrt(3) / 2

    x_side_T2 = np.linspace(0, 0.5, 100)
    x_side_T3 = np.linspace(-0.5, 0, 100)

    ax.plot(x_side_T2, twisstntern.T2(0, x_side_T2), "k", linewidth=1)
    ax.plot(x_side_T3, twisstntern.T3(0, x_side_T3), "k", linewidth=1)
    ax.hlines(y=0, xmin=-0.5, xmax=0.5, color="k", linewidth=1)

    ax.set_xticks([])
    ax.set_yticks([])

    for i in range(1, int(1 / alpha)):
        y = i * alpha

        ax.hlines(y=y * h, xmin=twisstntern.T1_lim(y)[0], xmax=twisstntern.T1_lim(y)[1], color="#C74375", linewidth=1)

        x2 = np.linspace(twisstntern.T2_lim(y)[0], twisstntern.T2_lim(y)[1], 100)
        ax.plot(x2, twisstntern.T2(y, x2), color="#3182BD", linewidth=1)

        x3 = np.linspace(twisstntern.T3_lim(y)[0], twisstntern.T3_lim(y)[1], 100)
        ax.plot(x3, twisstntern.T3(y, x3), color="#D4A017", linewidth=1)

    ax.vlines(x=0, ymin=0, ymax=h, colors="#888888", ls=':')

    ax.text(-0.02, 0.88, 'T1', size=12, color="#C74375")
    ax.text(0.54, -0.01, 'T3', size=12, color="#D4A017")
    ax.text(-0.58, -0.01, 'T2', size=12, color="#3182BD")

    coord = np.arange(0, 1 + alpha, alpha)
    T_1 = np.arange(0, 0.5 + alpha / 2, alpha / 2)
    T_2 = np.arange(-0.5, 0 + alpha / 2, alpha / 2)
    T_3 = np.arange(-0.5, 0.5 + alpha, alpha)

    for i in range(len(T_1)):
        label = str(round(1 - coord[i], 2))
        x = T_1[i]
        y = twisstntern.T2(0, x)
        ax.text(x + 0.01, y, label, size=7, color="#C74375")

    for i in range(len(T_2)):
        label = str(round(1 - coord[i], 2))
        x = T_2[i]
        y = twisstntern.T3(0, x)
        ax.text(x - 0.04, y, label, size=7, color="#3182BD")

    for i in range(len(coord)):
        label = str(round(coord[i], 2))
        x = T_3[i]
        ax.text(x, -0.03, label, size=7, color="#D4A017")

         
    # filling the subtriangles with the counts-
    numberPoints = results["nTruth"].sum()

    
    for index, row in results.iterrows():
        trianglex, triangley, direction = twisstntern.return_triangle_coord(
            row["T1"][0], row["T1"][1], row["T2"][0], row["T2"][1], row["T3"][0], row["T3"][1]
        )
        x = row["L2"] / ((numberPoints**2)/5)
        color = gradient_ignite_log(x)
        ax.fill(trianglex, triangley, color=color)

    ax.set_title(title, fontsize=14)
    
    add_plain_gradient_legend(ax)
    return ax


# # The examples from Sean:

# ## Investigating m

# In[32]:


# alpha=0.1;
# # Investigating m
# dataTruthm = data_process("TRUTH_m.csv")
# # List of file name identifiers
# name_list_num = ["A_m0","B_m0.001","C_m0.01_C","D_m0.02","E_m0.03",
#                  "F_m0.035","G_m0.04","H_m0.05","I_m0.06"]


# # Create an empty list to store results
# results = []

# # Loop through the list and calculate distances
# for i in name_list_num:
#     file_name = i + ".csv"
#     print(f"Processing: {file_name}")

#     # Process data and compute distance
#     data2 = data_process(file_name)

#     # COMPARING
    
#     # Wasserstein- KL cost matrix
#     distKL = simplex_wasserstein_distance(dataTruthm, data2)
#     # Wasserstein- euclidean cost matrix
#     distEU = wasserstein_distance(dataTruthm, data2)

#     # L2, chi^2 -statistic, p-value for chi^2
#     table, L2, chi_Statistic, p_value = compare_triangle_counts(dumpData(dataTruthm), dumpData(data2), alpha);
#     # Store results
#     results.append((file_name, distEU, distKL, L2, chi_Statistic, p_value ));

# # Create DataFrame
# df_m = pd.DataFrame(results, columns=["File Name", "Wasserstein", "Wasserstein-KL", "L2", "chi Statistic", "p-value"])


# # Extract m values
# df_m['m_value'] = df_m['File Name'].str.extract(r'm([\d\.]+)')[0].str.rstrip('.').astype(float)

# # Sort
# df_m = df_m.sort_values(by='m_value')

# #PLOTTING
# sns.set(style="whitegrid")

# # Melt the DataFrame
# df_melted = df_m.melt(
#     id_vars=["m_value"],
#     value_vars=["Wasserstein", "Wasserstein-KL", "L2", "chi Statistic", "p-value"],
#     var_name="Metric",
#     value_name="Value"
# )

# # Create subplots
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), sharex=True)

# # Plot 1: Wasserstein metrics
# for metric in ["Wasserstein", "Wasserstein-KL"]:
#     sns.lineplot(
#         data=df_melted[df_melted["Metric"] == metric],
#         x="m_value", y="Value", ax=ax1, marker="o", label=metric
#     )
# ax1.set_title("Wasserstein Distances")
# ax1.set_ylabel("Distance Value")
# ax1.set_xlabel("m-value")
# ax1.legend()

# # Plot L2 and Chi² Statistic on ax2
# for metric in ["L2", "chi Statistic"]:
#     sns.lineplot(
#         data=df_melted[df_melted["Metric"] == metric],
#         x="m_value", y="Value", ax=ax2, marker="s", label=metric
#     )

# ax2.set_yscale("log")
# ax2.set_ylabel("Value (log scale)")
# ax2.set_xlabel("m-value")
# ax2.set_title("$L^2$, $\chi^2$ Statistic and p-value")

# # Plot p-value on secondary y-axis (ax2b)
# ax2b = ax2.twinx()
# pval_line = sns.lineplot(
#     data=df_melted[df_melted["Metric"] == "p-value"],
#     x="m_value", y="Value", ax=ax2b, marker="^", color="black"
# )

# # Set label manually
# pval_line.legend_ = None  # Prevent Seaborn from messing with it
# pval_line_label = ax2b.plot([], [], color="black", marker="^", linestyle="none", label="p-value")[0]

# ax2b.set_yscale("log")
# ax2b.set_ylabel("p-value (log)", color="black")
# ax2b.tick_params(axis='y', labelcolor='black')

# # Combine legends from both axes and place them automatically
# lines1, labels1 = ax2.get_legend_handles_labels()
# lines2, labels2 = ax2b.get_legend_handles_labels()
# ax2.legend(lines1 + lines2, labels1 + labels2, loc='best')  


# In[33]:


# note to self, plot things and make sure this matches intution for the distnace


# ### Lets plot L2 distances

# In[34]:


# alpha=0.1
# dataTruthm = data_process("TRUTH_m.csv")
# # List of file name identifiers
# name_list_num = ["A_m0","B_m0.001","C_m0.01_C","D_m0.02","E_m0.03",
#                  "F_m0.035","G_m0.04","H_m0.05","I_m0.06"]

# fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(17, 12))  # 3x3 grid

# k=0;
# for i in name_list_num:
#     file_name = i + ".csv"
#     row, col = divmod(k, 3) 
#     ax = axes[row, col]

#     # Process data and compute distance
#     datam2 = data_process(file_name)
#     # calculating L2 per subtriangle, stored inside --table["L2"]
#     tablem, L2, chi_Statistic, p_value = compare_triangle_counts(dumpData(dataTruthm), dumpData(datam2), alpha);

#     title = "L2(Truth," + "m" +i+")" 
#     PlotL2(title,tablem, alpha,ax)
#     k += 1


# In[ ]:





# ## Investigating Ne

# In[35]:


# alpha=0.1;
# # Investigating Ne
# name_list_num = ["0.05","0.15","0.20","0.23","0.25","0.27","0.30","0.35","0.45"]
# dataTruth= data_process("TRUTH_ne0.25_1.csv")

# # Create an empty list to store results
# resultsN = []

# # Loop through the list and calculate distances
# for i in name_list_num:
#     file_name = "ne" +i + ".csv"
#     print(f"Processing: {file_name}")

#     # Process data and compute distance
#     dataNe2 = data_process(file_name)

#     # Wasserstein KL
#     distKLN = simplex_wasserstein_distance(dataTruth, dataNe2)
#     #Wasserstein Euclidean
#     distN = wasserstein_distance(dataTruth, dataNe2)
#     # L2, chi^2 -statistic, p-value for chi^2
#     table_N, L2_N, chi_Statistic_N, p_value_N = compare_triangle_counts(dumpData(dataTruth), dumpData(dataNe2), alpha);
    
#     # Store results
#     resultsN.append((file_name, distN, distKLN, L2_N, chi_Statistic_N, p_value_N ));



# # Create DataFrame
# df_N = pd.DataFrame(resultsN, columns=["File Name", "Wasserstein", "Wasserstein-KL", "L2", "chi Statistic", "p-value"])


# # Extract m values
# df_N['Ne_value'] = df_N['File Name'].str.extract(r'ne([\d\.]+)')[0].str.rstrip('.').astype(float)

# # Sort
# df_N = df_N.sort_values(by='Ne_value')

# #PLOTTING
# sns.set(style="whitegrid")

# # Melt the DataFrame
# df_melted = df_N.melt(
#     id_vars=["Ne_value"],
#     value_vars=["Wasserstein", "Wasserstein-KL", "L2", "chi Statistic", "p-value"],
#     var_name="Metric",
#     value_name="Value"
# )

# # Create subplots
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), sharex=True)

# # Plot 1: Wasserstein metrics
# for metric in ["Wasserstein", "Wasserstein-KL"]:
#     sns.lineplot(
#         data=df_melted[df_melted["Metric"] == metric],
#         x="Ne_value", y="Value", ax=ax1, marker="o", label=metric
#     )
# ax1.set_title("Wasserstein Distances")
# ax1.set_ylabel("Distance")
# ax1.set_xlabel("Ne-value")
# ax1.legend()

# # Plot L2 and Chi² Statistic on ax2
# for metric in ["L2", "chi Statistic"]:
#     sns.lineplot(
#         data=df_melted[df_melted["Metric"] == metric],
#         x="Ne_value", y="Value", ax=ax2, marker="s", label=metric
#     )

# ax2.set_yscale("log")
# ax2.set_ylabel("Value (log scale)")
# ax2.set_xlabel("Ne-value")
# ax2.set_title("$L^2$, $\chi^2$ Statistic and p-value")

# # Plot p-value on secondary y-axis (ax2b)
# ax2b = ax2.twinx()
# pval_line = sns.lineplot(
#     data=df_melted[df_melted["Metric"] == "p-value"],
#     x="Ne_value", y="Value", ax=ax2b, marker="^", color="black"
# )

# # Set label manually
# pval_line.legend_ = None  # Prevent Seaborn from messing with it
# pval_line_label = ax2b.plot([], [], color="black", marker="^", linestyle="none", label="p-value")[0]

# ax2b.set_yscale("log")
# ax2b.set_ylabel("p-value (log)", color="black")
# ax2b.tick_params(axis='y', labelcolor='black')

# # Combine legends from both axes and place them automatically
# lines1, labels1 = ax2.get_legend_handles_labels()
# lines2, labels2 = ax2b.get_legend_handles_labels()
# ax2.legend(lines1 + lines2, labels1 + labels2, loc='best')  



# In[36]:


# dataTruth= data_process("TRUTH_ne0.25_1.csv")
# data1=dumpData(dataTruth)
# data2=dumpData(data_process("ne0.30.csv"))
# data3=dumpData(data_process("ne0.45.csv"))
# table, L2, chi_Statistic, p_value = compare_triangle_counts(data1,data3,0.05)
# PlotL2("My first Compare",table, 0.1)


# In[37]:


# # OK LETS RUN OVER ALL DATA DIFFERENCES

# name_list_num = ["0.05","0.15","0.20","0.23","0.25","0.27","0.30","0.35","0.45"]
# dataTruth= data_process("TRUTH_ne0.25_1.csv")
# fig, axes = plt.subplots(nrows=3, ncols=3, figsize=(17, 12))  # 3x3 grid

# k=0;
# for i in name_list_num:
#     file_name = "ne" +i + ".csv"
#     row, col = divmod(k, 3) 
#     ax = axes[row, col]

#     # Process data and compute distance
#     dataNe2 = data_process(file_name)
#     # calculating L2 per subtriangle, stored inside --table["L2"]
#     table, L2, chi_Statistic, p_value = compare_triangle_counts(dumpData(dataTruth), dumpData(dataNe2), alpha);

#     title = "L2(Truth=Ne0.25," + "Ne" +i+")" 
#     PlotL2(title,table, 0.1,ax)
#     k += 1


    


# ## General Plotting

# In[38]:


# A sleek new data look
def plotData(data, alpha, ax=None):
    if ax is None:
        fig = plt.figure(figsize=(8, 6))
        ax = plt.axes()

    h = math.sqrt(3) / 2

    x_side_T2 = np.linspace(0, 0.5, 100)
    x_side_T3 = np.linspace(-0.5, 0, 100)

    ax.plot(x_side_T2, twisstntern.T2(0, x_side_T2), "k", linewidth=1)
    ax.plot(x_side_T3, twisstntern.T3(0, x_side_T3), "k", linewidth=1)
    ax.hlines(y=0, xmin=-0.5, xmax=0.5, color="k", linewidth=1)

    ax.set_xticks([])
    ax.set_yticks([])

    for i in range(1, int(1 / alpha)):
        y = i * alpha

        # T1 lines
        ax.hlines(
            y=y * h,
            xmin=twisstntern.T1_lim(y)[0],
            xmax=twisstntern.T1_lim(y)[1],
            color="#7B1E1E", linewidth=1
        )

        # T2 lines
        x2 = np.linspace(twisstntern.T2_lim(y)[0], twisstntern.T2_lim(y)[1], 100)
        ax.plot(x2, twisstntern.T2(y, x2), color="#277DA1", linewidth=1)

        # T3 lines
        x3 = np.linspace(twisstntern.T3_lim(y)[0], twisstntern.T3_lim(y)[1], 100)
        ax.plot(x3, twisstntern.T3(y, x3), color="#F4A261", linewidth=1)

    ax.vlines(x=0, ymin=0, ymax=h, colors="#3E3E3E", ls=':')

    x_data, y_data = twisstntern.cartizian(data["T1"], data["T2"], data["T3"])
    ax.scatter(
        x_data, y_data,
        facecolors='#A2C5F2',
        edgecolors='gray',
        alpha=0.7,
        linewidths=0.5,
        s=10
    )

    ax.text(-0.02, 0.88, 'T1', size=12, color="#D72638")
    ax.text(0.54, -0.01, 'T3', size=12, color="#F4A261")
    ax.text(-0.58, -0.01, 'T2', size=12, color="#277DA1")

    coord = np.arange(0, 1 + alpha, alpha)
    T_1 = np.arange(0, 0.5 + alpha / 2, alpha / 2)
    T_2 = np.arange(-0.5, 0 + alpha / 2, alpha / 2)
    T_3 = np.arange(-0.5, 0.5 + alpha, alpha)

    for i in range(len(T_1)):
        label = str(round(1 - coord[i], 2))
        x = T_1[i]
        y = twisstntern.T2(0, x)
        ax.text(x + 0.01, y, label, size=7, color="#D72638")

    for i in range(len(T_2)):
        label = str(round(1 - coord[i], 2))
        x = T_2[i]
        y = twisstntern.T3(0, x)
        ax.text(x - 0.04, y, label, size=7, color="#1F5F7F")

    for i in range(len(coord)):
        label = str(round(coord[i], 2))
        x = T_3[i]
        ax.text(x, -0.03, label, size=7, color="#E68A45")

    # Clean up axis frame
    for spine in ax.spines.values():
        spine.set_visible(False)

    return ax  # optional: return the axis if further use is needed


# In[39]:


#Plotting subtriangle heat map 
# color pallete is based on the number of triangles and data poitns
def gradient_gray_blue_red_purple(x, exptNrm):
    x = min(1, max(0, x))  # Clamp to [0, 1]

    if x <= exptNrm:
        # Stage 1: Light gray (0.98, 0.98, 0.98) → Gentle Blue (0.5, 0.72, 0.92)
        t = x / exptNrm
        r = 0.98 * (1 - t) + 0.5 * t   # 0.98 → 0.5
        g = 0.98 * (1 - t) + 0.72 * t   # 0.98 → 0.7
        b = 0.98 * (1 - t) + 0.92 * t   # 0.98 → 0.9

    elif x <= 3 * exptNrm:
        # Stage 2: Gentle Blue → Flashy Coral Red (0.87, 0.26, 0.2)
        t = (x - exptNrm) / (2 * exptNrm)
        r = 0.5 * (1 - t) + 0.87 * t   # 0.5 → 0.87
        g = 0.7 * (1 - t) + 0.26 * t   # 0.7 → 0.26
        b = 0.9 * (1 - t) + 0.20 * t   # 0.9 → 0.20

    else:
        # Stage 3: Coral Red (0.87, 0.26, 0.2) → Deep Purple (0.2, 0.0, 0.4)
        t = (x - 3 * exptNrm) / (1 - 3 * exptNrm)
        r = 0.87 * (1 - t) + 0.2 * t    # 0.87 → 0.2
        g = 0.26 * (1 - t) + 0.0 * t    # 0.26 → 0.0
        b = 0.20 * (1 - t) + 0.4 * t    # 0.20 → 0.4

    return (r, g, b)




def make_colormap_from_function(gradient_func, resolution=256):
    gradient = [gradient_func(x / (resolution - 1)) for x in range(resolution)]
    return ListedColormap(gradient)

########################################################################################################################
## for plotting only one plot, to see what xceddes normal subtriangle quantity- maybe this is confusing for comparison-- ask Sean's opinion    

def PlotSubtrianglePositiveCount2COLORS(data, alpha, ax=None):
    # Get triangle data
    triangles = SubtrianglesDataCount(data, alpha)

    # Create new axis if none is provided
    if ax is None:
        fig = plt.figure(figsize=(8, 6))
        ax = plt.axes()

    # coordinate skeleton    
    h = math.sqrt(3) / 2

    x_side_T2 = np.linspace(0, 0.5, 100)
    x_side_T3 = np.linspace(-0.5, 0, 100)

    ax.plot(x_side_T2, twisstntern.T2(0, x_side_T2), "k", linewidth=1)
    ax.plot(x_side_T3, twisstntern.T3(0, x_side_T3), "k", linewidth=1)
    ax.hlines(y=0, xmin=-0.5, xmax=0.5, color="k", linewidth=1)

    ax.set_xticks([])
    ax.set_yticks([])

    for i in range(1, int(1 / alpha)):
        y = i * alpha

        ax.hlines(y=y * h, xmin=twisstntern.T1_lim(y)[0], xmax=twisstntern.T1_lim(y)[1], color="#C74375", linewidth=1)

        x2 = np.linspace(twisstntern.T2_lim(y)[0], twisstntern.T2_lim(y)[1], 100)
        ax.plot(x2, twisstntern.T2(y, x2), color="#3182BD", linewidth=1)

        x3 = np.linspace(twisstntern.T3_lim(y)[0], twisstntern.T3_lim(y)[1], 100)
        ax.plot(x3, twisstntern.T3(y, x3), color="#D4A017", linewidth=1)

    ax.vlines(x=0, ymin=0, ymax=h, colors="#888888", ls=':')

    ax.text(-0.02, 0.88, 'T1', size=12, color="#C74375")
    ax.text(0.54, -0.01, 'T3', size=12, color="#D4A017")
    ax.text(-0.58, -0.01, 'T2', size=12, color="#3182BD")

    coord = np.arange(0, 1 + alpha, alpha)
    T_1 = np.arange(0, 0.5 + alpha / 2, alpha / 2)
    T_2 = np.arange(-0.5, 0 + alpha / 2, alpha / 2)
    T_3 = np.arange(-0.5, 0.5 + alpha, alpha)

    for i in range(len(T_1)):
        label = str(round(1 - coord[i], 2))
        x = T_1[i]
        y = twisstntern.T2(0, x)
        ax.text(x + 0.01, y, label, size=7, color="#C74375")

    for i in range(len(T_2)):
        label = str(round(1 - coord[i], 2))
        x = T_2[i]
        y = twisstntern.T3(0, x)
        ax.text(x - 0.04, y, label, size=7, color="#3182BD")

    for i in range(len(coord)):
        label = str(round(coord[i], 2))
        x = T_3[i]
        ax.text(x, -0.03, label, size=7, color="#D4A017")

        
    # filling the subtriangles with the counts
    numberPoints = triangles["n"].sum()
    numerTriangles = triangles.shape[0]
    expectedNormal = numerTriangles / numberPoints #basically normalized L1 count (by tottal number of points)

    for index, row in triangles.iterrows():
        trianglex, triangley, direction = twisstntern.return_triangle_coord(
            row["T1"][0], row["T1"][1], row["T2"][0], row["T2"][1], row["T3"][0], row["T3"][1]
        )
        x = row["n"] / numberPoints
        color = gradient_gray_blue_red_purple(x, exptNrm=expectedNormal)
        ax.fill(trianglex, triangley, color=color)

    cmap = make_colormap_from_function(lambda x: gradient_gray_blue_red_purple(x, exptNrm=expectedNormal))
    norm = plt.Normalize(vmin=0, vmax=1)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cax = inset_axes(ax, width="3%", height="80%", loc='center right', 
                     bbox_to_anchor=(0.05, 0, 1, 1), bbox_transform=ax.transAxes, borderpad=1)
    cbar = plt.colorbar(sm, cax=cax)
    cbar.set_label("Fraction of total points per subtriangle", fontsize=10)
    cax.axhline(expectedNormal, color='black', linestyle='--', linewidth=1)
    cax.text(1.0, expectedNormal, f"___\n{expectedNormal:.3f}", va='center', ha='left',
             transform=cax.transAxes, fontsize=8, color='black')

    return ax


# In[40]:


#putting both these plots together
#accepting right now only alpha==0.1 or 0.05, but can change it prinicipally
def plot_dual(fileName,data, alpha):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)

    # Left: Raw data
    plotData(data, alpha, ax=ax1);

    # Right: Heatmap
    if alpha == 0.1:
        PlotSubtrianglePositiveCount2COLORS(data, alpha, ax=ax2);
    else:     
        PlotSubtrianglePositiveCount(data, alpha, ax=ax2);

    # Hide spines from ax2
    for spine in ax2.spines.values():
        spine.set_visible(False)

    #Titles
    # One title to rule them all
    #fig.suptitle(rf'{fileName}', fontsize=17)
    fig.suptitle(rf'Dataset: $\mathrm{{{fileName}}}$', fontsize=18)
    ax1.set_title(r'$\mathrm{Raw\ Data}$', fontsize=11)
    ax2.set_title(r'$\mathrm{Subtriangle\ Counts}$', fontsize=11)



    return fig


# In[41]:


# file = "F_m0.035.csv"
# data1 = dumpData(data_process(file))

# file = "TRUTH_ne0.25_1.csv"

# dataTruth= data_process("TRUTH_ne0.25_1.csv")
# dataTruth=dumpData(dataTruth)
# data2= dataTruth


# In[42]:


#fig=plot_dual("Truth, ne=0.25",data1, 0.05)


# In[43]:


#fig=plot_dual("Ne=0.3",data2, 0.05)


# In[44]:


#fig=plot_dual("Truth, Ne=0.25",data1, 0.1)


# In[45]:


#fig=plot_dual("Ne=0.3",data2, 0.1)


# In[ ]:





# In[46]:


# table, L2, chi_Statistic, p_value = compare_triangle_counts(data1,data2,0.05)
# PlotL2("My first Compare",table, 0.05)


# ## Trying visulazing directly data1 - data2 (not even absoloute value)

# In[47]:


# CODE
######### AUXILARY###############################################
#################COLOR GRADIENT FUNCTIONS#################################################################
#THIS for point counts
def gradient_lightgray_to_darkred_to_purple(x):
    x = min(1, max(0, x))  # clamp between 0 and 1

    if x <= 0.2:
        # Stage 1: Light gray (0.988) → Dark Red (0.4, 0, 0)
        t = x / 0.1
        r = 0.988 - (0.58 * t)    # 0.988 → 0.4
        g = 0.988 * (1 - t)       # 0.988 → 0
        b = 0.988 * (1 - t)       # 0.988 → 0
    else:
        # Stage 2: Dark Red → Deep Purple (0.1, 0, 0.2)
        t = (x - 0.1) / 0.9
        r = 0.4 - 0.3 * t        # 0.4 → 0.1
        g = 0.0                  # stays 0
        b = 0.0 + 0.2 * t        # 0.0 → 0.2

    return (r, g, b)
    


# #THIS FOR MINUS- for differnece between point counts of two data sets
def gradient_gray_to_green_to_cyan_peak_at_0_1(x):
    x = min(1, max(0, x))  # Clamp to [0, 1]

    if x <= 0.1:
        # Stage 1: Light gray → dark jazzy green
        t = x / 0.1
        r = 0.988 * (1 - t)             # 0.988 → 0.0
        g = 0.988 - (0.58 * t)          # 0.988 → 0.4
        b = 0.988 - (0.78 * t)          # 0.988 → 0.2
    else:
        # Stage 2: Dark green → cyan
        t = (x - 0.1) / 0.9
        r = 0.0                        # stays 0
        g = 0.4 + 0.6 * t              # 0.4 → 1.0
        b = 0.2 + 0.8 * t              # 0.2 → 1.0

    return (r, g, b)

def diverging_custom_gradient(x):
    x = max(-1, min(1, x))  # Clamp x ∈ [-1, 1]
    
    if x < 0:
        # Use green side for negatives, mirrored
        return gradient_gray_to_green_to_cyan_peak_at_0_1(-x)
    elif x > 0:
        # Use red-purple side for positives
        return gradient_lightgray_to_darkred_to_purple(x)
    else:
        return (1.0, 1.0, 1.0)  # White at center


### Auxilary for heatmap#######################################

def make_diverging_colormap(resolution=256):
    # Map x from -1 to 1 → resolution points
    values = np.linspace(-1, 1, resolution)
    colors = [diverging_custom_gradient(x) for x in values]
    return ListedColormap(colors)
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def add_diverging_legend(ax, label="Difference (Normalized)"):
    cmap = make_diverging_colormap()
    norm = plt.Normalize(vmin=-1, vmax=1)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    cax = inset_axes(ax, width="3%", height="80%", loc='center right',
                     bbox_to_anchor=(0.05, 0, 1, 1), bbox_transform=ax.transAxes, borderpad=1)

    cbar = plt.colorbar(sm, cax=cax)
    cbar.set_label(label, fontsize=10)
    cbar.set_ticks([-1, -0.5, 0, 0.5, 1])
    cbar.set_ticklabels(["-1", "-0.5", "0", "0.5", "1"])
######################################################################################################################
######################################################################################################################

#Plotting the differnece between two data sets: data1-data2 in a heat map over subtriangles
#INPUT: the output of compare_triangle_counts(data1,data2,alpha), alpha


def PlotCountDifference(title,results, alpha,ax=None):
     # Create new axis if none is provided
    if ax is None:
        fig = plt.figure(figsize=(8, 6))
        ax = plt.axes()

    # coordinate skeleton    
    h = math.sqrt(3) / 2

    x_side_T2 = np.linspace(0, 0.5, 100)
    x_side_T3 = np.linspace(-0.5, 0, 100)

    ax.plot(x_side_T2, twisstntern.T2(0, x_side_T2), "k", linewidth=1)
    ax.plot(x_side_T3, twisstntern.T3(0, x_side_T3), "k", linewidth=1)
    ax.hlines(y=0, xmin=-0.5, xmax=0.5, color="k", linewidth=1)

    ax.set_xticks([])
    ax.set_yticks([])

    for i in range(1, int(1 / alpha)):
        y = i * alpha

        ax.hlines(y=y * h, xmin=twisstntern.T1_lim(y)[0], xmax=twisstntern.T1_lim(y)[1], color="#C74375", linewidth=1)

        x2 = np.linspace(twisstntern.T2_lim(y)[0], twisstntern.T2_lim(y)[1], 100)
        ax.plot(x2, twisstntern.T2(y, x2), color="#3182BD", linewidth=1)

        x3 = np.linspace(twisstntern.T3_lim(y)[0], twisstntern.T3_lim(y)[1], 100)
        ax.plot(x3, twisstntern.T3(y, x3), color="#D4A017", linewidth=1)

    ax.vlines(x=0, ymin=0, ymax=h, colors="#888888", ls=':')

    ax.text(-0.02, 0.88, 'T1', size=12, color="#C74375")
    ax.text(0.54, -0.01, 'T3', size=12, color="#D4A017")
    ax.text(-0.58, -0.01, 'T2', size=12, color="#3182BD")

    coord = np.arange(0, 1 + alpha, alpha)
    T_1 = np.arange(0, 0.5 + alpha / 2, alpha / 2)
    T_2 = np.arange(-0.5, 0 + alpha / 2, alpha / 2)
    T_3 = np.arange(-0.5, 0.5 + alpha, alpha)

    for i in range(len(T_1)):
        label = str(round(1 - coord[i], 2))
        x = T_1[i]
        y = twisstntern.T2(0, x)
        ax.text(x + 0.01, y, label, size=7, color="#C74375")

    for i in range(len(T_2)):
        label = str(round(1 - coord[i], 2))
        x = T_2[i]
        y = twisstntern.T3(0, x)
        ax.text(x - 0.04, y, label, size=7, color="#3182BD")

    for i in range(len(coord)):
        label = str(round(coord[i], 2))
        x = T_3[i]
        ax.text(x, -0.03, label, size=7, color="#D4A017")

         
    # filling the subtriangles with the counts-
    numberPoints = results["nTruth"].sum()

    
    for index, row in results.iterrows():
        trianglex, triangley, direction = twisstntern.return_triangle_coord(
            row["T1"][0], row["T1"][1], row["T2"][0], row["T2"][1], row["T3"][0], row["T3"][1]
        )
        x = row["data1-data2"]/ numberPoints
        color = diverging_custom_gradient(x)
        ax.fill(trianglex, triangley, color=color)

    ax.set_title(title, fontsize=14)
    
    add_diverging_legend(ax)
    return ax




# In[48]:


# fig=PlotCountDifference("Comparing Ne0.3 to 0.2",table, 0.05,ax=None)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




