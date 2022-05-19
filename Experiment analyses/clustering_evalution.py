
# coding: utf-8

import pandas as pd
import numpy as np
import time
start = time.time()

def validate_origin(original_data,true_label):
    print('start evaluating original_data....')
    lab_num=len(np.unique(true_label))
    print('Number of categories：',lab_num)
    pearson_sim = np.corrcoef(original_data)
    ##--------------louvain clustering-----------
    louvain = Louvain()
    lv_pre = louvain.fit_transform(pearson_sim)
    nmi_louvain = format(normalized_mutual_info_score(true_label, lv_pre), '.5f')
    ari_louvain = format(adjusted_rand_score(true_label, lv_pre), '.5f')
    print('nmi and ari of original data with louvain clustering: ',
          str(nmi_louvain), str(ari_louvain))
    ##------------kmeans clustering------------
    estimators = KMeans(n_clusters=lab_num)
    est = estimators.fit(original_data)
    kmeans_pred = est.labels_
    nmi_kmeans = format(normalized_mutual_info_score(true_label, kmeans_pred), '.5f')
    ari_kmeans = format(adjusted_rand_score(true_label, kmeans_pred), '.5f')
    print('nmi and ari of original data with kmeans clustering: ',
          str(nmi_kmeans), str(ari_kmeans))
    ##-------------pearson spectral clustering---------------
    sc_pred = sklearn.cluster.SpectralClustering(n_clusters=lab_num, affinity='precomputed').fit_predict(pearson_sim)
    nmi_spectral = format(normalized_mutual_info_score(true_label, sc_pred), '.5f')
    ari_spectral = format(adjusted_rand_score(true_label, sc_pred), '.5f')

    print('nmi and ari of original data with pearson spectral clustering: ',
          str(nmi_spectral), str(ari_spectral))

    return nmi_louvain, nmi_kmeans, nmi_spectral, ari_louvain, ari_kmeans, ari_spectral

end = time.time()
print ('time:',end-start)


start = time.time()
def validate_imputation(imputed_data, true_label):
    print('start evaluating imputation....')
    lab_num = len(np.unique(true_label))
    pearson_sim = np.corrcoef(imputed_data) 
    if np.min(np.min(pearson_sim))<0:
        pearson_sim[pearson_sim<0]=0
    ##--------------louvain clustering-----------
    louvain = Louvain()
    lv_pre = louvain.fit_transform(pearson_sim)
    nmi_louvain = format(normalized_mutual_info_score(true_label, lv_pre), '.5f')
    ari_louvain = format(adjusted_rand_score(true_label, lv_pre), '.5f')
    print('nmi and ari of imputed data with louvain clustering: ',
          str(nmi_louvain), str(ari_louvain))
    ##------------kmeans clustering------------
    estimators = KMeans(n_clusters=lab_num)
    est = estimators.fit(imputed_data)
    kmeans_pred = est.labels_
    nmi_kmeans = format(normalized_mutual_info_score(true_label, kmeans_pred), '.5f')
    ari_kmeans = format(adjusted_rand_score(true_label, kmeans_pred), '.5f')
    print('nmi and ari of imputed data with kmeans clustering: ',
          str(nmi_kmeans), str(ari_kmeans))
    ##-------------pearson spectral clustering---------------
    sc_pred = sklearn.cluster.SpectralClustering(n_clusters=lab_num, affinity='precomputed').fit_predict(pearson_sim)
    nmi_spectral = format(normalized_mutual_info_score(true_label, sc_pred), '.5f')
    ari_spectral = format(adjusted_rand_score(true_label, sc_pred), '.5f')

    print('nmi and ari of imputed data with pearson spectral clustering: ',
          str(nmi_spectral), str(ari_spectral))
    
    
    return nmi_louvain, nmi_kmeans, nmi_spectral, ari_louvain, ari_kmeans, ari_spectral

end = time.time()
print ('time:',end-start)  
