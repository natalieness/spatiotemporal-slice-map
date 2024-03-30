# -*- coding: utf-8 -*-
"""

Functions to handle data import and export

"""

def get_files(TP):
    df = pd.read_csv("%s.csv"%TP) 
    #df = df[~df.Label.str.contains("LeftSCN")]
    dimdf = pd.read_csv("dimensions.csv")
    n_rows = int(dimdf.columns[0])
    n_cols = int(dimdf.columns[1]) 
    n_rois = n_rows*n_cols
    n_slices = max(df['Slice'])
    rois = np.arange(1,n_rois+1)
    rois = list(rois)
    rois = rois*n_slices
    df['ROI'] = rois
    #to exclude certain frames from analysis
    #df = df.iloc[(n_rois*210):]
    #n_slices = n_slices-210
    
    return (df, n_rows, n_cols, n_rois, n_slices)

def save_cluster_TS(name, X, T, lab, ts_perslice, clustered_TS, no_k, n_slices, n_rows, n_cols, invert, detrend):
    X = X.flatten()
    if detrend == True:
        for i in range(clustered_TS.shape[0]):
            TS = clustered_TS[i,:]
            model = np.polyfit(X, TS, 3)
            predicted = np.polyval(model, X)
            dtr = TS-predicted
            clustered_TS[i,:] = dtr

    plt.figure()
    plt.plot(np.arange(0,n_slices), clustered_TS[1,:])

    biodare_arr = np.zeros([no_k,n_slices+1])
    biodare_arr[:,1:] = clustered_TS
    biodare_arr[:,0] = np.arange(1, no_k+1)

    biodare_f = pd.DataFrame(biodare_arr) 

    col_names = list(T) #frame time /2 and 0.5 for 30min
    col_names.insert(0,'')
    biodare_f.columns = col_names

    filepath = '%s_TS.xlsx'%name
    biodare_f.to_excel(filepath, index=False)
    np.savez_compressed('%s_clusterdata.npz'%name, lab = lab, clustered_TS = clustered_TS, ts_perslice = ts_perslice)

    if invert == True:
        for i in range(clustered_TS.shape[0]):
            inverted_series = clustered_TS[i,:]*(-1) + np.amax(clustered_TS[i,:])
            clustered_TS[i,:] = inverted_series
            
    X = X.flatten()
    if detrend == True:
        for i in range(clustered_TS.shape[0]):
            TS = clustered_TS[i,:]
            model = np.polyfit(X, TS, 3)
            predicted = np.polyval(model, X)
            dtr = TS-predicted
            clustered_TS[i,:] = dtr

    plt.figure()
    plt.plot(np.arange(0,n_slices), clustered_TS[1,:])

    biodare_arr = np.zeros([no_k,n_slices+1])
    biodare_arr[:,1:] = clustered_TS
    biodare_arr[:,0] = np.arange(1, no_k+1)

    biodare_f = pd.DataFrame(biodare_arr) 

    col_names = list(T) 
    col_names.insert(0,'')
    biodare_f.columns = col_names

    if invert == True:
        filepath = '%s_TS_inv.xlsx'%name
        biodare_f.to_excel(filepath, index=False)
        np.savez_compressed('%s_clusterdata_inv.npz'%name, lab = lab, clustered_TS = clustered_TS, ts_perslice = ts_perslice)
        
        

def import_BioDare_results_finalcluster_ksmall(name, new_k):
    biod_res = pd.read_csv('%s_BioDare.csv'%name, skiprows = range(0, 22),error_bad_lines=False)
    
    biod_res = biod_res.iloc[:new_k,:]

    # get npy files 
    with np.load('%s_clusterdata.npz'%name) as data:
        print(data.files)
        GFP_lab = data['lab']
        GFP_clustered_TS = data['clustered_TS'][:,:]
            
        print(len(GFP_clustered_TS[1,:]))
        
    return (biod_res, GFP_lab, GFP_clustered_TS) 


