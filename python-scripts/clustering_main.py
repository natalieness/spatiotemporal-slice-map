""" 

Script for identifying clusters of activity in longitudinal recordings
of biological data. Used and tested on bioluminescence and
fluorescence data from rhythmic reporters in organotypic brain slices. 


Can handle multiple samples in batch 

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

from data_handling_functions import get_files, get_timeseries_for_clustering, save_cluster_TS
from plotting_functions import reset_cmap, save_fig, plot_clustered_ts_results
from map_functions import Cluster_Timeseries


#set matplotlib parameters
plt.rcParams.update({'font.size': 20,
                     'axes.labelsize': 22,
                     'axes.linewidth': 1,
                     'xtick.major.size': 10,
                     'ytick.major.size': 10})

#%% Data import and basic parameters and clustering 

# parameters to set for analysis 
# for multiple files in batch, add to list, also works with just 1 value
#filenames of time series data extracted from ImageJ , 1 per sample.
TPs = ['example_ts1','example_ts2'] #csv suffix not needed
#Names of files for saving 
names = ['Example_1','Example_2']
#whether to invert the time series before clustering
inverts = [False,False]


interval = 0.5 #acquisition interval in hours 
no_k = 6 #number of clusters
detrend = False #detrend time series first



#%% K-Means Clustering 

cmap = reset_cmap() #get or reset colormap 

#iterate through samples 
for h in range(len(TPs)):
    TP = TPs[h] #get filename for import of sample
    name = names[h] #name for saving
    invert = inverts[h] 

    #import sample data
    (df, n_rows, n_cols, n_rois, n_slices) = get_files(TP) 
    #get time points
    T = np.arange(0,(n_slices/(1/interval)),interval)
    
    #remove dips in time series
    (TS_list, ts_perslice, X) = get_timeseries_for_clustering(T, df, n_rois, n_slices, no_k, cluster_norm=True)
    
        
    #perform clustering
    kmeans = Cluster_Timeseries(no_k, TS_list)
    #plot and save results
    (clustered_TS, lab) = plot_clustered_ts_results(T, kmeans, no_k, n_rows, n_cols, n_slices, ts_perslice, cmap_name, colors)
    save_cluster_TS(name, X, T, lab, ts_perslice, clustered_TS, no_k, n_slices, n_rows, n_cols, invert, detrend)
