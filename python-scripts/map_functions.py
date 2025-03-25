# -*- coding: utf-8 -*-
"""
Functions for data processing and mapping 

"""
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from matplotlib.colors import LinearSegmentedColormap

def plot_avg_groupby_pd(df):
    a = np.array(df.groupby('ROI')['Mean'])
    arr = np.zeros([a[0,1].shape[0], a.shape[0]])
    for h in range(a.shape[0]):
        arr[:,h] = a[h,1][:]
    arr_m = np.mean(arr, axis=1)
    #plot mean trace 
    plt.plot(np.arange(0,arr_m.shape[0]), arr_m)
    plt.show()
    
    return arr_m


def get_timeseries_for_clustering(T, df, n_rois, n_slices, no_k, cluster_norm=True):
    # set up X values 
    X = np.array(T)
    X = X.reshape(-1,1)
    
    #normalisation
    TS_list = []
    TS_list_norm = []
    
    #exclude background clusters 
    #exclude ROIs here based on the standard deviation of the timeseries
    
    for i in range(1,n_rois+1):
        means = df[df['ROI'] == i]['Mean']
        means_max = np.amax(means)
        means_min = np.amin(means)
        means_range = means_max - np.amin(means)
        if means_range > 0:
            means_norm = [(i - means_min)/means_range for i in means]
        elif means_range == 0:
            means_norm = [i-means_max for i in means]
        means_nonnorm = [i for i in means]
        if cluster_norm == True:
            TS_list.append(means_norm)
        elif cluster_norm == False:
            TS_list.append(means_nonnorm)
        else:
            print('Error')
        TS_list_norm.append(means_norm)
    


    ts_perslice = []
    for k in range(n_slices):
        ts_ = []
        for k_ in range(len(TS_list)):
            ts_.append(TS_list_norm[k_][k])
        ts_perslice.append(ts_)
        
    return (TS_list, ts_perslice, X)

def Cluster_Timeseries(no_k, TS_list):
    #KMeans clustering
    kmeans = KMeans(n_clusters=no_k, init = 'random', random_state=None, algorithm='lloyd').fit(TS_list)
    #full is lloyd in newer versions
    return kmeans

def check_and_format_BioDare_results_finalcluster(biod_res, cmap, invert, phasetype, GFP_lab):
    
    print(biod_res[phasetype])
    
    #find background cluster 
    bg_found = False
    for i in range(len(biod_res[phasetype])):
        if 'FAILED' in str(biod_res['Status'][i]):
            pos_BG = i
            bg_found = True
    if bg_found == False:
        c_ignored = 0 
        for i in range(len(biod_res[phasetype])):
            if 'IGNORED' in str(biod_res['Status'][i]):
                c_ignored += 1
        if c_ignored == 1: 
            for i in range(len(biod_res[phasetype])):
                if 'IGNORED' in str(biod_res['Status'][i]):
                    pos_BG=i
                    bg_found = True
    if bg_found == False:
        #if no clusters failed period estimation, assume lowest amplitude
        pos_BG = biod_res['Amplitude'].argmin()
    
    #positions excluding bg 
    pos_av = [x for x in list(np.arange(6)) if x != pos_BG ]


    #find non-rhythmic regions - ignored in biodare dataset
    NR = []
    for i in pos_av:
        if 'IGNORED' in str(biod_res['Status'][i]):
            NR.append(i)
    
    #update ignored regions 
    pos_av = [x for x in pos_av if x not in NR]
    #get number of non-rhythmic clusters in total incl bg 
    n_nan = 6-len(pos_av)
        
    #get average period if traces need inversion
    if invert == True:
        period = float(biod_res['Period'][pos_av].mean())
        print("Period: %f"%period)
    
    #get list of detected phases 
    reclustered_phases = []
    for i in pos_av:
        a = float(biod_res[phasetype][i])
        if (invert == True) & (np.isnan(biod_res[phasetype][i]) == False):
            a = a-(period/2)
            if a <0:
                a = a+24
        b = np.floor(a)
        c = (a-b)*0.6+b
        d = round(c,2)
        reclustered_phases.append(d)
        
    #repeat for period 
    pers = []
    for i in pos_av:
        a = float(biod_res['Period'][i])
        b = np.floor(a)
        c = (a-b)*0.6+b
        d = round(c,2)
        pers.append(d)
    
    #initiate list to keep track of orders
    grand_order = list(np.arange(6))
    grand_order = [pos_BG] + NR + pos_av
    
    #sort through rhythmic clusters 
    if len(reclustered_phases) > 1:
        reclustered_phases_O = list(np.argsort(reclustered_phases))
        print('Original positions:')
        for g in range(len(reclustered_phases)):
            print(g, reclustered_phases[g] )
        print('New order:', np.array(reclustered_phases)[reclustered_phases_O])
    
        correct_order = int(input('Is the order of phases correct? Answer 1 if yes, 0 if no '))
        if correct_order == False: 
            new_order = input("Please enter correct order of ROIs, with no spaces or commas. Example: 012345. ")
            new_order = list(map(int,new_order))

        else:
            new_order = reclustered_phases_O
    
        grand_order = [pos_BG] + NR + list(np.array(pos_av)[new_order])
        labels_ordered = ['NR'] + ['NR']*len(NR) + list(np.array(reclustered_phases)[new_order])
        pers_labels_ordered = ['NR'] + ['NR']*len(NR) + list(np.array(pers)[new_order])
    else:
        grand_order = [pos_BG] + NR + pos_av
        labels_ordered = ['NR'] + ['NR']*len(NR) + reclustered_phases
        pers_labels_ordered = ['NR'] + ['NR']*len(NR) + pers

        
    #initiate ordered lists and map
    GFP_lab_ordered = np.zeros(GFP_lab.shape)
    reclustered_phases_ordered = []
    pers_ordered=[]
    
    #order rhythmic clusters and assign label to map 
    for e,i in enumerate(grand_order):
        idx = np.where(GFP_lab == i)
        GFP_lab_ordered[idx] = e
        reclustered_phases_ordered.append(labels_ordered[e])
        pers_ordered.append(pers_labels_ordered[e])
    

    reclustered_phases_str = ['%.2f'%i  if e>(n_nan-1) else i for e,i in enumerate(reclustered_phases_ordered)]
    per_str = ['%.2f'%i  if e>(n_nan-1) else i for e,i in enumerate(pers_ordered)]
    #stop indent
        
    if n_nan > 1: 
        print('entered multi-non rhythmic loop')
        white = (0.7, 0.7, 0.7,1.0)
        n_c = -(6-n_nan) #number of normal color clusters, negative to use as index
        colors = [(0.0, 0.0, 0.0, 1.0),(0.0, 0.0, 0.8, 1.0),(0.46875, 0.0, 1.0, 1.0),(1.0, 0.3600000000000002, 0.6399999999999999, 1.0),(1.0, 0.7600000000000001, 0.24, 1.0)]  
        colorsn = [*list([colors[0]]), *list([white]*(n_nan-1)), *list(colors[n_c:])]
        if n_nan ==6:
            colorsn = [*list([colors[0]]), *list([white]*(n_nan-1))]
        n_bin = 6  # Discretizes the interpolation into bins
        cmap_name = 'gnuplot_nowhite'
        cmap = LinearSegmentedColormap.from_list(cmap_name, colorsn, N=n_bin)
    else:
        cmap=cmap
    
    return (cmap, labels_ordered, reclustered_phases_str, GFP_lab_ordered, grand_order, per_str)






