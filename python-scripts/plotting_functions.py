# -*- coding: utf-8 -*-
"""

Functions for plotting 

"""

import pandas as pd
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.cluster import KMeans
from scipy import signal
from scipy.spatial.distance import cdist
from sklearn.metrics import silhouette_samples, silhouette_score
from matplotlib.mlab import detrend_linear

from statsmodels.tsa.tsatools import detrend
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker

import os 
import csv



def reset_cmap():
    colors = [(0.0, 0.0, 0.0, 1.0),(0.0, 0.0, 0.8, 1.0),(0.46875, 0.0, 1.0, 1.0),(1.0, 0.3600000000000002, 0.6399999999999999, 1.0),(1.0, 0.7600000000000001, 0.24, 1.0)]  
    n_bin = 6  # Discretizes the interpolation into bins
    cmap_name = 'gnuplot_nowhite'
    cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)
    return cmap

def save_fig(fig, name):
    fig.savefig('%s.svg'%(name), format='svg', bbox_inches='tight')
    
def plot_clustered_ts_results(T, kmeans, no_k, n_rows, n_cols, n_slices, ts_perslice, cmap_name, colors):
    #prep colormap
    cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=no_k)

    lab = kmeans.labels_
    
    unique_lab = np.unique(lab)
    
    lab = lab.reshape([n_rows, n_cols])
    
    #plot map
    fig, ax = plt.subplots()
    ax1 = ax.imshow(lab, cmap=cmap)
    plt.axis('off')

    #plot time series associated with each cluster
    plt.figure()
    clustered_TS = np.zeros([no_k, n_slices])
    for i in unique_lab:
        indices = np.where(kmeans.labels_ == i)
        mean_ts_i = []
        for j in range(n_slices):
            mean_ts_i.append(np.mean(np.array(ts_perslice[j])[indices[0]]))
        clustered_TS[i,:] = mean_ts_i
        #if i <(no_k-1):
        plt.plot(T[:], mean_ts_i[:], color=cmap(i))
        plt.plot(T[:], mean_ts_i[:], color=cmap(i))
        #elif i >= (no_k-1):
           # plt.plot(Xvals[:96], mean_ts_i[10:106], color='grey')
        #if i == 4:
         #   plt.plot(Xvals[:48], mean_ts_i[10:58], color='red')
        ax.set_ylabel('Normalized fluorescence (A.U.)')
        ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
    plt.xlabel('Time (h)')
    #plt.spines['right'].set_visible(False)
    #plt.spines['top'].set_visible(False)
    return (clustered_TS, lab)


def imshow_specific_cluster(lab, cmap):
    b = int(input("Show specific cluster? Answer 1"))
    while b == True: 
        choice = np.unique(lab)
        print(choice)
        choice_input = int(input("Choose a cluster out of choice shown in console "))
        wanted_cluster = (lab==choice_input)*np.ones(lab.shape)
        plt.figure()
        plt.imshow(wanted_cluster, cmap=cmap)
        plt.axis('off')
        plt.show()
        b = int(input("Answer 1 if you want to see another cluster "))
        
        
def plot_final_map_k6(T, GFP_clustered_TS, cmap, labels_ordered, reclustered_phases_str, GFP_lab_ordered, grand_order, interval, per_str, name):
    # check if T and length of timeseries match 
    if len(T) != len(GFP_clustered_TS[1,:]):
        print('X values and length of timeseries do not match. Make sure you are using corresponding data.')
        new_nslices = len(GFP_clustered_TS[0,:])
        T = np.arange(0,(new_nslices/(1/interval)),interval)
        print(len(T))
        
    #determine the number of cycles to be plotted 
    possible_cycles = (len(GFP_clustered_TS[0])-10)/(24/interval)
    
    lindtr = True
    
    #check if figures should be autosaved 
    sf = int(input('Should figures be automatically saved in current folder? Press 1 if yes. 0 otherwise'))
    
    fig = plt.figure(figsize=(8,5))
    ax = plt.subplot(111)
    
    #plot timeseries in pre-determined order 
    for e,i in enumerate(grand_order):
        if lindtr==True:
            print(i)
            print(len(T))
            ax.plot(T[:], detrend_linear(GFP_clustered_TS[i,:]), color=cmap(e), linewidth=3)
        else:
            ax.plot(T[:], GFP_clustered_TS[i,:], color=cmap(e), linewidth=3)
    
    #ax.set_title('Cluster timeseries')
    ax.set_ylabel('Normalized \n fluorescence (A.U.)', fontsize=40)
    ax.set_xlabel('Time (h)', fontsize=40)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
    ax.tick_params(labelsize=35)


    plt.show()
    if sf == 1:
        save_fig(fig, '%s_TS'%name)
    
    repeat = int(input("Press 1 to replot using different range. 0 otherwise "))
    while repeat == 1: 
        print("Possible number of cycles to plot: %f"%possible_cycles)
        n_cyc = int(input("How many cycles to do you want to plot? Enter integer "))
        skip_at_start = float(input("How many cycles to do you want to skip at the start? Enter integer or float (if interval <1h)  "))
        n_X = int(n_cyc*(24/interval))
        SP = int(skip_at_start*(24/interval))
        EP = SP+n_X
        Iskip = 0
        SP = SP+Iskip
        EP = EP+Iskip
        
        
        fig = plt.figure(figsize=(8,5))
        ax = plt.subplot(111)
        
        for e,i in enumerate(grand_order):
            if lindtr==True:
                ax.plot(T[:n_X], detrend_linear(GFP_clustered_TS[i,(SP):(EP)]), color=cmap(e), linewidth=3)
            else:
                ax.plot(T[:n_X], GFP_clustered_TS[i,(SP):(EP)], color=cmap(e), linewidth=3)
        
        
        ax.set_ylabel('Normalized \n fluorescence (A.U.)', fontsize=40)
        ax.set_xlabel('Time (h)', fontsize=40)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
        ax.tick_params(labelsize=35)
        plt.show()
        if sf == 1:
            save_fig(fig, '%s_TS%i'%(name,n_cyc))
        
        repeat = int(input("Press 1 to replot using different range. 0 otherwise "))
    
    
    fig, ax = plt.subplots()
    ax1 = ax.imshow(GFP_lab_ordered, cmap=cmap)
    num_ticks = 6
    cbar = fig.colorbar(ax1, ticks = np.linspace(0.5, 4.5, num_ticks))
    cbar.set_ticklabels(reclustered_phases_str)
    cbar.ax.invert_yaxis()
    plt.axis('off') 
    if sf == 1:
        save_fig(fig, '%s_imshow'%name)
    
    # plot corresponding periods 
    fig, ax = plt.subplots()
    ax1 = ax.imshow(GFP_lab_ordered, cmap=cmap)
    num_ticks = 6
    cbar = fig.colorbar(ax1, ticks = np.linspace(0.5, 4.5, num_ticks))
    cbar.set_ticklabels(per_str)
    cbar.ax.invert_yaxis()
    cbar.ax.set_title('Period (h)', loc='left')
    plt.axis('off') 
    if sf == 1:
        save_fig(fig, '%s_period'%name)
        



