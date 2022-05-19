
# coding: utf-8

import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

dataset = ['sc_10x_5cl']

#输入真实标签
label_path = "./label/sc_10x_5cl_labelName.txt"
labels = pd.read_csv(label_path, delimiter='\t',header=None)
labels=labels.values
labels=labels.flatten()
labels=list(labels)
print(len(labels),type(labels))

impute_mtd=["Original","Bubble","SAVER","scImpute","VIPER","bayNorm","ALRA","DCA","autoimpute"]
# Create the scatter plot 
#plt.figure(figsize=(30, 10))
# for i in range(0,len(impute_mtd)):
for i in range(0,len(impute_mtd)):
    if impute_mtd[i]=="Original":
        imputedData=pd.read_csv('./processed_data/'+"scranNorm"+'/'+ dataset[0] +'_genebycell.txt',sep='\t',header=0,index_col=0)
        imputedData=imputedData.T
    elif impute_mtd[i] =="Bubble":
        imputedData=pd.read_csv('./imputed_data/'+"autoimpute"+'/'+ dataset[0] +'_20.txt',sep='\t',header=0,index_col=0)
        imputedData=imputedData.T
    elif impute_mtd[i]=="scImpute":
        imputedData=pd.read_csv('./imputed_data/'+"scimpute"+'/'+ dataset[0] +'_genebycell_count.txtscimpute_count.txt',sep=' ',header=0,index_col=0)
        imputedData=imputedData.T
    elif impute_mtd[i]=="autoimpute":
        imputedData=pd.read_csv('./imputed_data/'+impute_mtd[i]+'/'+ dataset[0] +'_genebycell.txt',sep=' ',header=0,index_col=0)
        imputedData=imputedData.T
    else:
        imputedData=pd.read_csv('./imputed_data/'+impute_mtd[i]+'/'+ dataset[0] +'_genebycell.txt',sep='\t',header=0,index_col=0)
        imputedData=imputedData.T
    print(imputedData.shape,type(imputedData))
          
    ####plot TSNE using scanpy ######
    ann = sc.AnnData(imputedData)
    ann.obs_names = list(imputedData.index)
    ann.var_names = list(imputedData.columns)
    ann.obs[impute_mtd[i]] = labels
    sc.tl.tsne(ann)
    sc.pl.tsne(ann, color=impute_mtd[i],legend_loc='on data',legend_fontsize=10)  
    plt.show() 
    print(impute_mtd[i]+"finish plot")