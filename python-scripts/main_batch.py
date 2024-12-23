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

#%% Getting timeseries of pre-determined clusters 
replotting = False
if replotting==True:
    #name for new data 
    name1 ='C2inC1BSL'
    # Files for reporter to be mapped 
    TM = 'RBM3_32_dtr' #file prefix for import. note: suffix is _measure.csv
    (df, n_rows, n_cols, n_rois, n_slices) = get_files(TM)
    
    no_k = 6 
    invert = False # cluster 
    invert1 = False # to be remapped
    detrend = False
    interval = 0.5 #acquisition in hours  
    exclude_BG=False
    
    #Import clusters 
    name = 'C2-BSL_k6'
    (biod_res, lab) = import_BioDare_results_finalcluster_ksmall(name, no_k)[:2]
    
    T = np.arange(0,(n_slices/(1/interval)),interval)
    
    #(df, T, n_slices) = remove_dips(df, T, n_slices) 
    (TS_list, ts_perslice, X, only_SCN) = get_timeseries_for_clustering(T, df, n_rois, n_slices, no_k, cluster_norm=True, exclude_BG=exclude_BG)
    clustered_TS = plot_on_existing_Clusters(T, lab, no_k, n_rows, n_cols, n_slices, ts_perslice, cmap_name, colors, only_SCN, exclude_BG = False)
    save_cluster_TS(name1, X, T, lab, ts_perslice, clustered_TS, no_k, n_slices, n_rows, n_cols, invert1, detrend)

#%% Show specific clusters if desired
#imshow_specific_cluster(lab, cmap)

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
    
    if replotting == False:
        if new_k == 6: #no reclustering. based on GFP_lab and normal clusters deter2mined in part 1
            (biod_res, GFP_lab, GFP_clustered_TS) = import_BioDare_results_finalcluster_ksmall(name, new_k)
            (cmap, labels_ordered, reclustered_phases_str, GFP_lab_ordered, grand_order, per_str) = check_and_format_BioDare_results_finalcluster(biod_res, cmap, invert, phasetype)
            plot_final_map_k6(T, GFP_clustered_TS, cmap, labels_ordered, reclustered_phases_str, GFP_lab_ordered, grand_order, interval, per_str, name)
            ctdiff = replot_CT(adjust, desired_ct, ctdiff, reclustered_phases_str, phasetype, per_str, GFP_lab_ordered, name)
            
    if replotting==True:
        (biod_res, GFP_lab, GFP_clustered_TS) = import_BioDare_results_finalcluster_ksmall(name, new_k)
        #import phases and timeseries of replotted reporter
        (biod_res1, _, GFP_clustered_TS1) = import_BioDare_results_finalcluster_ksmall(name1, new_k)
        
        (cmap, labels_ordered, labels_ordered1, reclustered_phases_str, reclustered_phases_str1, GFP_lab_ordered, grand_order, per_str, per_str1) = REPLOT_check_and_format_BioDare_results_finalcluster(biod_res, biod_res1, cmap, invert, invert1, phasetype)
        
        replot_final_map_k6(T, GFP_clustered_TS, GFP_clustered_TS1, cmap, labels_ordered, labels_ordered1, reclustered_phases_str, reclustered_phases_str1, GFP_lab_ordered, grand_order, interval, per_str, per_str1)
        plot_difference(GFP_lab_ordered, reclustered_phases_str,reclustered_phases_str1)
        ctdiff = replot_CT(adjust, desired_ct, ctdiff, reclustered_phases_str1, phasetype, per_str1, GFP_lab_ordered, name1)
