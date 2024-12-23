# -*- coding: utf-8 -*-
"""

Main script to map spatiotemporal clusters of activity within organotypic SCN slices. 

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

from data_handling_functions import get_files, get_timeseries_for_clustering, save_cluster_TS, import_BioDare_results_finalcluster_ksmall
from plotting_functions import reset_cmap, save_fig, plot_clustered_ts_results, plot_final_map_k6 
from map_functions import Cluster_Timeseries, check_and_format_BioDare_results_finalcluster

#set matplotlib parameters
plt.rcParams.update({'font.size': 20,
                     'axes.labelsize': 22,
                     'axes.linewidth': 1,
                     'xtick.major.size': 10,
                     'ytick.major.size': 10})

#%% Data import, basic parameters and clustering 
cmap = reset_cmap() #reset colormap

#set filenames to use, and destination filenames (can be the same)
TPs = ['C1-BSL', 'C2-BSL']
names = ['C1-BSL', 'C2-BSL']
#invert timeseries if desired
inverts = [False,False]

#additional paramaters
no_k = 6
detrend = False 
interval = 0.5 #acquisition in hours 
exclude_BG = False

#iterate through files and get clusters
for h in range(len(TPs)):
    #get data
    TP = TPs[h] #file prefix for import. note: suffix is _measure.csv
    (df, n_rows, n_cols, n_rois, n_slices) = get_files(TP)
    name = '%s_k%i'%(names[h],no_k)
    invert = inverts[h]
    #set time series X values
    T = np.arange(0,(n_slices/(1/interval)),interval)

    #remove dips in fluorscence or bioluminescence signal if required
    #(df, T, n_slices) = remove_dips_on_mean(df, T, n_slices)
    #get timeseries to be used for clustering
    (TS_list, ts_perslice, X, only_SCN) = get_timeseries_for_clustering(T, df, n_rois, n_slices, no_k, cluster_norm=True, exclude_BG=exclude_BG)

    #clustering
    kmeans = Cluster_Timeseries(no_k, TS_list)
    #plot cluster time series and save
    (clustered_TS, lab) = plot_clustered_ts_results(T, kmeans, no_k, n_rows, n_cols, n_slices, ts_perslice, cmap_name, colors, only_SCN, exclude_BG = exclude_BG)
    save_cluster_TS(name, X, T, lab, ts_perslice, clustered_TS, no_k, n_slices, n_rows, n_cols, invert, detrend)


#%% Post-BioDare analysis 
""" Data need to be analysed in biodare before running this part 
and type of analysis needs to be chosen in parameters """
cmap = reset_cmap() #reset colormap 

if no_k < 10:
    new_k = 6
elif no_k > 10:
    new_k = 5

phasetype = 'Abs Phase To Zero'

names = ['C1_BSL','C2-BSL']
inverts = [False,False]
#decide whether to adjust to a certain mean CT, e.g. set RCaMP phase in SCN to CT6
adjusts = [1,0]
desired_ct=6

for h in range(0,len(names)):
    cmap = reset_cmap() #reset colormap
    name = '%s_k6'%names[h]
    print(name)
    invert=inverts[h]
    adjust=adjusts[h]
    if adjust == 1: 
        ctdiff=0
    

    if new_k == 6: #no reclustering. based on GFP_lab and normal clusters deter2mined in part 1
        (biod_res, GFP_lab, GFP_clustered_TS) = import_BioDare_results_finalcluster_ksmall(name, new_k)
        (cmap, labels_ordered, reclustered_phases_str, GFP_lab_ordered, grand_order, per_str) = check_and_format_BioDare_results_finalcluster(biod_res, cmap, invert, phasetype)
        plot_final_map_k6(T, GFP_clustered_TS, cmap, labels_ordered, reclustered_phases_str, GFP_lab_ordered, grand_order, interval, per_str, name)
        ctdiff = replot_CT(adjust, desired_ct, ctdiff, reclustered_phases_str, phasetype, per_str, GFP_lab_ordered, name)
            
