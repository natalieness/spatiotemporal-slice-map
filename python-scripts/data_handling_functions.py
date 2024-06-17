# -*- coding: utf-8 -*-
"""

Functions to handle data import and export

"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def get_files(TP, dimensions_file='dimensions.csv'):
    try:
        df = pd.read_csv(f"{TP}.csv")
        dimdf = pd.read_csv(dimensions_file)
    except FileNotFoundError as e:
        raise FileNotFoundError(f"File not found: {e.filename}")
    
    try:
        n_rows = int(dimdf.columns[0])
        n_cols = int(dimdf.columns[1])
    except ValueError:
        raise ValueError("Dimension columns should contain integer values.")

    n_rois = n_rows * n_cols
    n_slices = df['Slice'].max()
    
    rois = np.arange(1, n_rois + 1).tolist() * n_slices
    if len(rois) != len(df):
        raise ValueError("Mismatch between number of ROIs and dataframe rows.")
    
    df['ROI'] = rois

    return df, n_rows, n_cols, n_rois, n_slices

def detrend_time_series(X, time_series):
    for i in range(time_series.shape[0]):
        TS = time_series[i, :]
        model = np.polyfit(X, TS, 3)
        predicted = np.polyval(model, X)
        time_series[i, :] = TS - predicted
    return time_series

def invert_time_series(time_series):
    for i in range(time_series.shape[0]):
        time_series[i, :] = time_series[i, :] * (-1) + np.amax(time_series[i, :])
    return time_series

def save_time_series_to_excel(filepath, clustered_TS, time_points):
    biodare_arr = np.zeros([clustered_TS.shape[0], clustered_TS.shape[1] + 1])
    biodare_arr[:, 1:] = clustered_TS
    biodare_arr[:, 0] = np.arange(1, clustered_TS.shape[0] + 1)

    biodare_f = pd.DataFrame(biodare_arr)
    col_names = [''] + list(time_points)
    biodare_f.columns = col_names
    biodare_f.to_excel(filepath, index=False)

def save_cluster_TS(name, X, T, lab, ts_perslice, clustered_TS, invert=False, detrend=False):
    X = X.flatten()

    if detrend:
        clustered_TS = detrend_time_series(X, clustered_TS)

    if clustered_TS.shape[0] < 2:
        raise ValueError("Clustered time series array has less than 2 rows, cannot plot the second row.")

    plt.figure()
    plt.plot(np.arange(clustered_TS.shape[1]), clustered_TS[1, :])

    save_time_series_to_excel(f'{name}_TS.xlsx', clustered_TS, T)
    np.savez_compressed(f'{name}_clusterdata.npz', lab=lab, clustered_TS=clustered_TS, ts_perslice=ts_perslice)

    if invert:
        clustered_TS = invert_time_series(clustered_TS)
        if detrend:
            clustered_TS = detrend_time_series(X, clustered_TS)

        plt.figure()
        plt.plot(np.arange(clustered_TS.shape[1]), clustered_TS[1, :])
        
        save_time_series_to_excel(f'{name}_TS_inv.xlsx', clustered_TS, T)
        np.savez_compressed(f'{name}_clusterdata_inv.npz', lab=lab, clustered_TS=clustered_TS, ts_perslice=ts_perslice)

def import_BioDare_results_finalcluster_ksmall(name, new_k):
    try:
        biod_res = pd.read_csv(f'{name}_BioDare.csv', skiprows=range(0, 22), on_bad_lines='skip')
    except FileNotFoundError as e:
        raise FileNotFoundError(f"File not found: {e.filename}")

    biod_res = biod_res.iloc[:new_k, :]

    with np.load(f'{name}_clusterdata.npz') as data:
        print(data.files)
        GFP_lab = data['lab']
        GFP_clustered_TS = data['clustered_TS'][:, :]
        
        print(len(GFP_clustered_TS[1, :]))

    return biod_res, GFP_lab, GFP_clustered_TS
