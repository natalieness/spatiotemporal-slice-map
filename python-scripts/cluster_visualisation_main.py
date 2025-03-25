
""" 
This code is used to visualise the clusters of 
activity across the organotypic slice. 

Data input is clusters and time series from the previous analysis, 
and the output of data analysed in BioDare. If analysing circadian 
parameters a different way, the function import_BioDare_results will 
need to be modified to handle different data formats. 


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
import sys

from data_handling_functions import get_files, import_BioDare_results_finalcluster_ksmall
from plotting_functions import reset_cmap, save_fig, plot_clustered_ts_results, plot_final_map_k6
from map_functions import check_and_format_BioDare_results_finalcluster

#set matplotlib parameters
plt.rcParams.update({'font.size': 20,
                     'axes.labelsize': 22,
                     'axes.linewidth': 1,
                     'xtick.major.size': 10,
                     'ytick.major.size': 10})

#%% 
cmap = reset_cmap() #reset colormap 


phasetype = 'Abs Phase To Zero'

#filenames of time series data extracted from ImageJ , 1 per sample.
TPs = ['example_ts2'] #csv suffix not needed
#Names of files for saving 
names = ['Example_2']
#whether to invert the time series before clustering
inverts = [False,False]


no_k = 6
detrend = False
interval = 0.5 #acquisition in hours 


TP = TPs[0] 

#check if file exists 
if not os.path.isfile('%s.csv'%TP):
    #navigate to subdirectory if needed
    os.chdir('example_data_plot')

    # Verify the current working directory
    print("Current Directory:", os.getcwd())
    
    #check if file exists in subdirectory
    if not os.path.isfile('%s.csv'%TP):
        print('File not found. Check the filename and path.')
        sys.exit()

(df, n_rows, n_cols, n_rois, n_slices) = get_files(TP)
T = np.arange(0,(n_slices/(1/interval)),interval)

for h in range(0,len(names)):
    cmap = reset_cmap() #reset colormap
    name = '%s_k6'%names[h]
    print(name)
    invert=inverts[h]

    (biod_res, GFP_lab, GFP_clustered_TS) = import_BioDare_results_finalcluster_ksmall(name, 6)
    (cmap, labels_ordered, reclustered_phases_str, GFP_lab_ordered, grand_order, per_str) = check_and_format_BioDare_results_finalcluster(biod_res, cmap, invert, phasetype, GFP_lab)
    plot_final_map_k6(T, GFP_clustered_TS, cmap, labels_ordered, reclustered_phases_str, GFP_lab_ordered, grand_order, interval, per_str, name)