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


TPs = ['C2-BSL_dtr','C1-BSL', 'C3-BSL_dtr']
names = ['RFP_BSL_dtr','PER2_BSL','GFP_BSL_dtr']
inverts = [False,False,True]
#Navigate to the correct folder first 

no_k = 6
detrend = False
interval = 0.5 #acquisition in hours 
find_optimal_k = False #apply elbow method to find optimal k 
exclude_BG = False


for h in range(len(TPs)):
    TP = TPs[h] #file prefix for import. note: suffix is _measure.csv
    #TP = 'PER2_BSL_Right_dtr'
    (df, n_rows, n_cols, n_rois, n_slices) = get_files(TP)
    name = '%s_k%i'%(names[h],no_k)
    invert = inverts[h]
    T = np.arange(0,(n_slices/(1/interval)),interval)
    
    #(df, T, n_slices) = remove_dips(df, T, n_slices) 
    (df, T, n_slices) = remove_dips_on_mean(df, T, n_slices)
    (TS_list, ts_perslice, X, only_SCN) = get_timeseries_for_clustering(T, df, n_rois, n_slices, no_k, cluster_norm=True, exclude_BG=exclude_BG)
    
    if find_optimal_k == True:
        #optimal_k = find_k_elbow(TS_list)
        #no_k = optimal_k
        exclude_BG = True
        name = '%s_silh'%name
        silhouette_save = get_silhouette_score(TS_list)
        save_silhouette(silhouette_save, name)
        
    else:  
        kmeans = Cluster_Timeseries(no_k, TS_list)
        (clustered_TS, lab) = plot_clustered_ts_results(T, kmeans, no_k, n_rows, n_cols, n_slices, ts_perslice, cmap_name, colors, only_SCN, exclude_BG = exclude_BG)
        save_cluster_TS(name, X, T, lab, ts_perslice, clustered_TS, no_k, n_slices, n_rows, n_cols, invert, detrend)

#%% Getting timeseries of pre-determined clusters 
#name for new data 
name1 ='RBM3inRFP_32_dtr_k6'
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
name = 'RFP_32_dtr_k6'
(biod_res, lab) = import_BioDare_results_finalcluster_ksmall(name, no_k)[:2]

T = np.arange(0,(n_slices/(1/interval)),interval)

#(df, T, n_slices) = remove_dips(df, T, n_slices) 
(TS_list, ts_perslice, X, only_SCN) = get_timeseries_for_clustering(T, df, n_rois, n_slices, no_k, cluster_norm=True, exclude_BG=exclude_BG)
clustered_TS = plot_on_existing_Clusters(T, lab, no_k, n_rows, n_cols, n_slices, ts_perslice, cmap_name, colors, only_SCN, exclude_BG = False)
save_cluster_TS(name1, X, T, lab, ts_perslice, clustered_TS, no_k, n_slices, n_rows, n_cols, invert1, detrend)

#%% Show specific clusters if desired
imshow_specific_cluster(lab, cmap)

#%% change name only 
name1 ='RBM3inRFP_BSL_k6'
TM = 'C2_BSL_dtr'
name = 'RFP_BSL_dtr_k6'
invert1=False
invert=False
#%% change name only no replotting 
name = 'RBM3_BSL_k6'
invert = False
#%% Post-BioDare analysis 
""" Data need to be analysed in biodare before running this part 
and type of analysis needs to be chosen in parameters """
cmap = reset_cmap() #reset colormap 

if no_k < 10:
    new_k = 6
elif no_k > 10:
    new_k = 5
replotting = False
phasetype = 'Abs Phase To Zero'

names = ['RFP_BSL_dtr','PER2_BSL','GFP_BSL_dtr']
inverts = [False,False,True]
adjusts = [1,0,0]

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


#%% Replot with desired CT

#set parameters for replotting with adjusted CT 
adjust = 0#1 to set phases to desired ct, 0 to adjust them by a previously given ctdiff 
desired_ct = 6
if adjust == 1: 
    ctdiff=0

if replotting == False:
    ctdiff = replot_CT(adjust, desired_ct, ctdiff, reclustered_phases_str, phasetype, per_str, GFP_lab_ordered, name)

elif replotting==True:
    ctdiff = replot_CT(adjust, desired_ct, ctdiff, reclustered_phases_str1, phasetype, per_str1, GFP_lab_ordered, name1)
#%% save all silhouette scores in one document

cwd = os.getcwd()
print(cwd)

silmaster = 'C:/Users/natal/Dropbox (UK Dementia Research Institute)/Brancaccio Lab/Natalie/Microglia-Astrocyte crosstalk project/Spatiotemporal mapping/SilhouetteScoresMaster.csv'

with open(silmaster, 'a', newline="") as sm:
    writer = csv.writer(sm)
    #writer.writerow([cwd, name])
    for i in range(silhouette_save.shape[0]):
        writer.writerow([cwd, name, silhouette_save[i,0], silhouette_save[i,1]])