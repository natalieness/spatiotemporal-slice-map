# -*- coding: utf-8 -*-
"""

Functions for plotting 

"""
import os
import csv
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


def reset_cmap(n_bin: int = 6) -> LinearSegmentedColormap:
    """Reset colormap to predefined colors."""
    colors = [
        (0.0, 0.0, 0.0, 1.0),
        (0.0, 0.0, 0.8, 1.0),
        (0.46875, 0.0, 1.0, 1.0),
        (1.0, 0.36, 0.64, 1.0),
        (1.0, 0.76, 0.24, 1.0)
    ]
    cmap_name = 'gnuplot_nowhite'
    cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)
    return cmap


def save_fig(fig: plt.Figure, name: str):
    """Save a matplotlib figure to a file."""
    fig.savefig(f'{name}.svg', format='svg', bbox_inches='tight')


def plot_clustered_ts_results(T: np.ndarray, kmeans: KMeans, no_k: int, n_rows: int, n_cols: int,
                              n_slices: int, ts_perslice: list, cmap_name: str, colors: list):
    """Plot clustered timeseries results and return clustered time series and labels."""
    cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=no_k)
    lab = kmeans.labels_.reshape([n_rows, n_cols])

    fig, ax = plt.subplots()
    ax.imshow(lab, cmap=cmap)
    plt.axis('off')

    clustered_TS = np.zeros([no_k, n_slices])
    plt.figure()
    for i in np.unique(kmeans.labels_):
        indices = np.where(kmeans.labels_ == i)
        mean_ts_i = [np.mean(np.array(ts_perslice[j])[indices[0]]) for j in range(n_slices)]
        clustered_TS[i, :] = mean_ts_i
        plt.plot(T, mean_ts_i, color=cmap(i))
        ax.set_ylabel('Normalized fluorescence (A.U.)')
        ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
    plt.xlabel('Time (h)')
    plt.show()
    return clustered_TS, lab


def imshow_specific_cluster(lab: np.ndarray, cmap: LinearSegmentedColormap):
    """Display specific clusters based on user input."""
    while int(input("Show specific cluster? Answer 1 to continue, 0 to exit: ")):
        choice = np.unique(lab)
        print(choice)
        choice_input = int(input("Choose a cluster from the shown choices: "))
        wanted_cluster = (lab == choice_input) * np.ones(lab.shape)
        plt.figure()
        plt.imshow(wanted_cluster, cmap=cmap)
        plt.axis('off')
        plt.show()


def plot_final_map_k6(T: np.ndarray, GFP_clustered_TS: np.ndarray, cmap: LinearSegmentedColormap,
                      labels_ordered: list, reclustered_phases_str: list, GFP_lab_ordered: np.ndarray,
                      grand_order: list, interval: float, per_str: list, name: str):
    """Plot final map for k=6 clusters and save figures based on user input."""
    if len(T) != len(GFP_clustered_TS[1, :]):
        print('X values and length of timeseries do not match. Adjusting X values...')
        new_nslices = len(GFP_clustered_TS[0, :])
        T = np.arange(0, (new_nslices / (1 / interval)), interval)
        print(len(T))

    possible_cycles = (len(GFP_clustered_TS[0]) - 10) / (24 / interval)
    lindtr = True

    sf = int(input('Should figures be automatically saved in current folder? Press 1 for yes, 0 otherwise: '))

    fig = plt.figure(figsize=(8, 5))
    ax = plt.subplot(111)

    for e, i in enumerate(grand_order):
        ts = detrend_linear(GFP_clustered_TS[i, :]) if lindtr else GFP_clustered_TS[i, :]
        ax.plot(T, ts, color=cmap(e), linewidth=3)

    ax.set_ylabel('Normalized \n fluorescence (A.U.)', fontsize=40)
    ax.set_xlabel('Time (h)', fontsize=40)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
    ax.tick_params(labelsize=35)
    plt.show()

    if sf == 1:
        save_fig(fig, f'{name}_TS')

    while int(input("Press 1 to replot using different range, 0 otherwise: ")):
        print(f"Possible number of cycles to plot: {possible_cycles}")
        n_cyc = int(input("Enter number of cycles to plot: "))
        skip_at_start = float(input("Enter number of cycles to skip at the start: "))
        n_X = int(n_cyc * (24 / interval))
        SP = int(skip_at_start * (24 / interval))
        EP = SP + n_X

        fig = plt.figure(figsize=(8, 5))
        ax = plt.subplot(111)

        for e, i in enumerate(grand_order):
            ts = detrend_linear(GFP_clustered_TS[i, SP:EP]) if lindtr else GFP_clustered_TS[i, SP:EP]
            ax.plot(T[:n_X], ts, color=cmap(e), linewidth=3)

        ax.set_ylabel('Normalized \n fluorescence (A.U.)', fontsize=40)
        ax.set_xlabel('Time (h)', fontsize=40)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.xaxis.set_major_locator(ticker.MultipleLocator(24))
        ax.tick_params(labelsize=35)
        plt.show()
        if sf == 1:
            save_fig(fig, f'{name}_TS{n_cyc}')

    fig, ax = plt.subplots()
    ax1 = ax.imshow(GFP_lab_ordered, cmap=cmap)
    num_ticks = 6
    cbar = fig.colorbar(ax1, ticks=np.linspace(0.5, 4.5, num_ticks))
    cbar.set_ticklabels(reclustered_phases_str)
    cbar.ax.invert_yaxis()
    plt.axis('off')
    if sf == 1:
        save_fig(fig, f'{name}_imshow')

    fig, ax = plt.subplots()
    ax1 = ax.imshow(GFP_lab_ordered, cmap=cmap)
    cbar = fig.colorbar(ax1, ticks=np.linspace(0.5, 4.5, num_ticks))
    cbar.set_ticklabels(per_str)
    cbar.ax.invert_yaxis()
    cbar.ax.set_title('Period (h)', loc='left')
    plt.axis('off')
    if sf == 1:
        save_fig(fig, f'{name}_period')
