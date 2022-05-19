
# coding: utf-8

# In[ ]:


import os 
import h5py 
import sklearn 
import datetime 
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
from os.path import join 
from functools import partial 
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

#We input the preprocessed data. 
rawdata=pd.read_csv('./data/',sep='\t',header=None,index_col=None)
data_norm = rawdata.values  
data_norm1 = data_norm.copy()

#perform PCA
num_components=30 #The number of principal components at least explains 40% of the variance in the data
pca = PCA(n_components=num_components,svd_solver = "randomized")
pca_data = pca.fit_transform(data_norm)
print(pca_data.shape)

#K-means
clusters=len(label_num)#set the number of cluster by calculating the silhouette coefficient.
kmeans = KMeans(n_clusters=clusters, random_state=0).fit(pca_data)
label_pr =  kmeans.labels_

def find_cluster_cell_idx(l, label):
    '''
        find cells whose label == l
        return bool_idx
    '''
    return label==l 
    
#The identification of dropout events
st = datetime.datetime.now()
def identify_dropout(cluster_cell_idxs, X):
    for idx in cluster_cell_idxs:
        dropout=(X[:,idx]==0).sum(axis=1)/(X[:,idx].shape[1])
        dropout_thr=0.5
        dropout_upper_thr,dropout_lower_thr = np.nanquantile(dropout,q=dropout_thr),np.nanquantile(dropout,q=0)
        gene_index1 = (dropout<=dropout_upper_thr)&(dropout>=dropout_lower_thr)
        #print(gene_index1)
        cv = X[:,idx].std(axis=1)/X[:, idx].mean(axis=1)
        cv_thr=0.5
        cv_upper_thr,cv_lower_thr = np.nanquantile(cv,q=cv_thr),np.nanquantile(cv,q=0)
        #print(cv_upper_thr,cv_lower_thr)
        gene_index2 = (cv<=cv_upper_thr)&(cv>=cv_lower_thr)
        #print(gene_index2)
        #include_faslezero_gene= list(np.intersect1d(gene_index1,gene_index2))
        include_faslezero_gene = np.logical_and(gene_index1, gene_index2)
        #print(list(include_faslezero_gene).count(True))
        tmp = X[:, idx]
        tmp[include_faslezero_gene] = tmp[include_faslezero_gene]+(tmp[include_faslezero_gene]==0)*-1
        X[:, idx] = tmp
    return X

label_set = np.unique(label_pr)
cluster_cell_idxs = list(map(partial(find_cluster_cell_idx,label=label_pr), label_set))
data_identi=identify_dropout(cluster_cell_idxs, X=data_norm.T)

ed = datetime.datetime.now()   
print('identify_dropout ：', (ed - st).total_seconds())






