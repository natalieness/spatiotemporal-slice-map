# -*- coding: utf-8 -*-
"""
Helper functions to convert time formats, angle formats, etc. 

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


def time_to_dec(tm):
    '''
    converts all elements of a list from real time to a decimal

    '''
    tm = [float(i) for i in tm] #make sure all elements are floats
    dec = [(i-np.floor(i))/0.6+np.floor(i) for i in tm]
    return dec 

def dec_to_time(dec):
    '''
    converts all elements of a list from a decimal to real time

    '''
    dec = [float(i) for i in dec] #make sure all elements are floats
    tm = [(i-np.floor(i))*0.6+np.floor(i) for i in dec]
    return tm 


