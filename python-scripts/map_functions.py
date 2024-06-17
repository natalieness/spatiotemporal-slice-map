# -*- coding: utf-8 -*-
"""
Functions for data processing and mapping 

"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from matplotlib.colors import LinearSegmentedColormap


def plot_avg_groupby_pd(df: pd.DataFrame) -> np.ndarray:
    """
    Plot the average of 'Mean' values grouped by 'ROI' and return the averaged array.

    Parameters:
    df (pd.DataFrame): DataFrame containing 'ROI' and 'Mean' columns.

    Returns:
    np.ndarray: Averaged array of 'Mean' values.
    """
    grouped_means = df.groupby('ROI')['Mean'].apply(list).values
    mean_matrix = np.vstack(grouped_means)
    avg_means = mean_matrix.mean(axis=0)
    
    plt.plot(np.arange(len(avg_means)), avg_means)
    plt.show()
    
    return avg_means


def get_timeseries_for_clustering(T: np.ndarray, df: pd.DataFrame, n_rois: int, n_slices: int, cluster_norm: bool = True):
    """
    Prepare timeseries data for clustering.

    Parameters:
    T (np.ndarray): Array of time points.
    df (pd.DataFrame): DataFrame containing 'ROI' and 'Mean' columns.
    n_rois (int): Number of ROIs.
    n_slices (int): Number of slices.
    cluster_norm (bool): Flag to determine if normalization should be applied.

    Returns:
    tuple: Timeseries list, timeseries per slice, and reshaped time points.
    """
    X = T.reshape(-1, 1)
    
    TS_list = []
    TS_list_norm = []
    
    for roi in range(1, n_rois + 1):
        means = df[df['ROI'] == roi]['Mean'].values
        means_range = means.max() - means.min()
        if means_range > 0:
            means_norm = (means - means.min()) / means_range
        else:
            means_norm = means - means.max()
        
        if cluster_norm:
            TS_list.append(means_norm)
        else:
            TS_list.append(means)
        
        TS_list_norm.append(means_norm)
    
    ts_perslice = np.array(TS_list_norm).T.tolist()
    
    return TS_list, ts_perslice, X


def cluster_timeseries(no_k: int, TS_list: list) -> KMeans:
    """
    Perform KMeans clustering on the timeseries data.

    Parameters:
    no_k (int): Number of clusters.
    TS_list (list): List of timeseries data.

    Returns:
    KMeans: Fitted KMeans object.
    """
    kmeans = KMeans(n_clusters=no_k, init='random', algorithm='full').fit(TS_list)
    return kmeans


def find_non_rhythmic_regions(biod_res: pd.DataFrame, pos_av: list) -> list:
    """
    Find non-rhythmic regions in the BioDare dataset.

    Parameters:
    biod_res (pd.DataFrame): BioDare results DataFrame.
    pos_av (list): List of available positions.

    Returns:
    list: List of non-rhythmic regions.
    """
    return [i for i in pos_av if 'IGNORED' in str(biod_res['Status'][i])]


def calculate_period(biod_res: pd.DataFrame, pos_av: list) -> float:
    """
    Calculate the average period.

    Parameters:
    biod_res (pd.DataFrame): BioDare results DataFrame.
    pos_av (list): List of available positions.

    Returns:
    float: Average period.
    """
    
    period = biod_res['Period'][pos_av].mean()

    return period


def recluster_phases(biod_res: pd.DataFrame, pos_av: list, period: float, invert: bool, phasetype: str) -> list:
    """
    Recluster the phases based on the BioDare results.

    Parameters:
    biod_res (pd.DataFrame): BioDare results DataFrame.
    pos_av (list): List of available positions.
    period (float): Average period.
    invert (bool): Flag to determine if inversion is needed.
    phasetype (str): Type of phase to consider.

    Returns:
    list: Reclustered phases.
    """
    reclustered_phases = []
    for i in pos_av:
        phase = float(biod_res[phasetype][i])
        if invert and not np.isnan(phase):
            phase -= period / 2
            if phase < 0:
                phase += 24
        reclustered_phases.append(round((phase - np.floor(phase)) * 0.6 + np.floor(phase), 2))
    return reclustered_phases


def create_grand_order(pos_BG: int, NR: list, pos_av: list, reclustered_phases: list) -> list:
    """
    Create the grand order of clusters.

    Parameters:
    pos_BG (int): Background position.
    NR (list): List of non-rhythmic regions.
    pos_av (list): List of available positions.
    reclustered_phases (list): Reclustered phases.

    Returns:
    list: Grand order of clusters.
    """
    grand_order = [pos_BG] + NR + pos_av
    if len(reclustered_phases) > 1:
        reclustered_phases_O = np.argsort(reclustered_phases)
        correct_order = int(input('Is the order of phases correct? Answer 1 if yes, 0 if no '))
        if not correct_order:
            new_order = list(map(int, input("Please enter correct order of ROIs, with no spaces or commas. Example: 012345. ")))
        else:
            new_order = reclustered_phases_O
        grand_order = [pos_BG] + NR + list(np.array(pos_av)[new_order])
    return grand_order


def generate_ordered_labels(grand_order: list, GFP_lab: np.ndarray, biod_res: pd.DataFrame, pos_av: list, NR: list, reclustered_phases: list, pers: list, n_nan: int) -> tuple:
    """
    Generate ordered labels and period strings.

    Parameters:
    grand_order (list): Grand order of clusters.
    GFP_lab (np.ndarray): Array of GFP labels.
    biod_res (pd.DataFrame): BioDare results DataFrame.
    pos_av (list): List of available positions.
    NR (list): List of non-rhythmic regions.
    reclustered_phases (list): Reclustered phases.
    pers (list): List of periods.
    n_nan (int): Number of non-rhythmic clusters.

    Returns:
    tuple: Ordered labels, reclustered phases strings, GFP lab ordered, and period strings.
    """
    GFP_lab_ordered = np.zeros_like(GFP_lab)
    labels_ordered = ['NR'] + ['NR'] * len(NR) + reclustered_phases
    pers_labels_ordered = ['NR'] + ['NR'] * len(NR) + pers

    for e, i in enumerate(grand_order):
        idx = np.where(GFP_lab == i)
        GFP_lab_ordered[idx] = e

    reclustered_phases_str = [f'{i:.2f}' if e >= n_nan else i for e, i in enumerate(labels_ordered)]
    pers_str = [f'{i:.2f}' if e >= n_nan else i for e, i in enumerate(pers_labels_ordered)]
    
    return labels_ordered, reclustered_phases_str, GFP_lab_ordered, pers_str


def check_and_format_BioDare_results_finalcluster(
    biod_res: pd.DataFrame, pos_av: list, pos_BG: int, cmap: LinearSegmentedColormap, invert: bool, phasetype: str
) -> tuple:
    """
    Check and format BioDare results for final clustering.

    Parameters:
    biod_res (pd.DataFrame): BioDare results DataFrame.
    pos_av (list): List of available positions.
    pos_BG (int): Background position.
    cmap (LinearSegmentedColormap): Colormap.
    invert (bool): Flag to determine if inversion is needed.
    phasetype (str): Type of phase to consider.

    Returns:
    tuple: Updated colormap, ordered labels, reclustered phases, GFP lab ordered, grand order, and period strings.
    """
    print(biod_res[phasetype])
    
    NR = find_non_rhythmic_regions(biod_res, pos_av)
    pos_av = [x for x in pos_av if x not in NR]
    n_nan = 6 - len(pos_av)
    
    period = calculate_period(biod_res, pos_av)
    reclustered_phases = recluster_phases(biod_res, pos_av, period, invert, phasetype)
    pers = [round((float(biod_res['Period'][i]) - np.floor(float(biod_res['Period'][i]))) * 0.6 + np.floor(float(biod_res['Period'][i])), 2) for i in pos_av]
    
    grand_order = create_grand_order(pos_BG, NR, pos_av, reclustered_phases)
    labels_ordered, reclustered_phases_str, GFP_lab_ordered, pers_str = generate_ordered_labels(grand_order, GFP_lab, biod_res, pos_av, NR, reclustered_phases, pers, n_nan)
    
    if n_nan > 1:
        white = (0.7, 0.7, 0.7, 1.0)
        colors = [(0.0, 0.0, 0.0, 1.0), (0.0, 0.0, 0.8, 1.0), (0.46875, 0.0, 1.0, 1.0), (1.0, 0.36, 0.64, 1.0), (1.0, 0.76, 0.24, 1.0)]
        colorsn = [colors[0]] + [white] * (n_nan - 1) + colors[-(6 - n_nan):]
        if n_nan == 6:
            colorsn = [colors[0]] + [white] * (n_nan - 1)
        cmap = LinearSegmentedColormap.from_list('gnuplot_nowhite', colorsn, N=6)
    
    return cmap, labels_ordered, reclustered_phases_str, GFP_lab_ordered, grand_order, pers_str



