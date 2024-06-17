# -*- coding: utf-8 -*-
"""
Helper functions to convert time formats, angle formats, etc. 

"""

import numpy as np

def time_to_dec(times: list[float]) -> list[float]:
    """
    Converts all elements of a list from real time (HH.MM) to decimal hours.
    
    Parameters:
    times (list of float): List of times in HH.MM format.
    
    Returns:
    list of float: List of times in decimal hours.
    """
    times = np.array(times, dtype=float)
    decimal_hours = np.floor(times) + (times % 1) / 0.6
    return decimal_hours.tolist()


def dec_to_time(decimals: list[float]) -> list[float]:
    """
    Converts all elements of a list from decimal hours to real time (HH.MM).
    
    Parameters:
    decimals (list of float): List of times in decimal hours.
    
    Returns:
    list of float: List of times in HH.MM format.
    """
    decimals = np.array(decimals, dtype=float)
    real_times = np.floor(decimals) + (decimals % 1) * 0.6
    return real_times.tolist()

