# spatiotemporal-slice-map
Tools to map spatiotemporal circadian activity across brain tissue slices from longitudinal recordings

Includes macro to get grid-like ROI time series from image stacks in ImageJ 
Scripts to cluster ROI time series in Python and then visualise with circadian phase information
Current implementation involves getting circadian parameters using biodare2.com before visualising the clusters with circadian phase labels

To cluster ROI time series in Python, run clustering_main.py and then cluster_visualisation_main.py to visualise phase maps.
