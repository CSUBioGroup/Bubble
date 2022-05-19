# coding: utf-8
import pandas as pd
import numpy as np
import scipy 
from scipy import io
from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler()
####plot TSNE using scanpy ######
dataset = ["sim_data3",'sim_data1.1','sim_data1.2',"sim_data2",'sim_data2.1','sim_data2.2',"sim_data1",'sim_data3.1','sim_data3.2',"sim_data4",'sim_data4.1','sim_data4.2']
impute_mtd=["Dropout","Bubble","SAVER","scImpute","VIPER","bayNorm","ALRA","AutoImpute","DCA","SCRABBLE"]
for j in range(0,len(dataset)):
    print("****************************"+dataset[j]+" "+"start computation"+"****************************")
    error_list=[]
    if dataset[j] in ["sim_data1","sim_data2","sim_data3","sim_data4","sim_data5"]:
        TrueData=pd.read_csv("./dataset/"+ dataset[j] +'_true.txt',sep='\t',header=0,index_col=0)
        TrueData=TrueData.T
        TrueData=TrueData.values
        print(type(TrueData),TrueData.shape)
    else: 
        TrueData=pd.read_csv("./dataset/"+ dataset[j] +'_true.txt',sep='\t',header=0,index_col=0)
        TrueData=TrueData.values
        print(type(TrueData),TrueData.shape)
    for i in range(0,len(impute_mtd)):
        if impute_mtd[i]=="Dropout":
            imputedData=pd.read_csv("./dataset/"+ dataset[j] +'_row_col.txt',sep='\t',header=0,index_col=0)
            imputedData=imputedData.T
            imputedData=imputedData.values
            imputedData = scaler.fit_transform(imputedData) 
        elif impute_mtd[i] =="Bubble":
            imputedData=pd.read_csv('./imputed_data/'+"autoimpute"+'/'+ dataset[j] +'.txt',sep='\t',header=0,index_col=0)
            imputedData=imputedData.values
            imputedData=imputedData.T
            imputedData = scaler.fit_transform(imputedData)   
        elif impute_mtd[i]=="scImpute":
            imputedData=pd.read_csv('./imputed_data/'+"scimpute"+'/'+ dataset[j] +'.txt',sep=' ',header=0,index_col=0)
            imputedData=imputedData.values
            imputedData=imputedData.T
            imputedData = scaler.fit_transform(imputedData)          
        elif impute_mtd[i]=="SCRABBLE":
            imputedData=pd.read_csv('./imputed_data/'+"scrabble"+'/'+ dataset[j] +'.txt',sep=' ',header=0,index_col=0)
            imputedData=imputedData.values
            imputedData=imputedData.T
            imputedData = scaler.fit_transform(imputedData)            
        elif impute_mtd[i]=="AutoImpute":
            imputedata_path = "./imputed_data/" + "auto" + "/" + dataset[j] + '.mat'
            imputedData= scipy.io.loadmat(imputedata_path)
            imputedData= list(imputedData.values())[-1]
            imputedData=np.squeeze(imputedData)
            imputedData=pd.DataFrame(imputedData)
            imputedData=imputedData.values
            imputedData = scaler.fit_transform(imputedData)       
        elif impute_mtd[i]=="DCA":  
            imputedata_path = "./imputed_data/" + 'DCA' + "/" + dataset[j] + '_results'+"/" +"mean.tsv"
            imputedData=pd.read_csv(imputedata_path,sep='\t',header=0,index_col=0)
            imputedData=imputedData.values
            imputedData=imputedData.T
            imputedData = scaler.fit_transform(imputedData)         
        else:
            imputedData=pd.read_csv('./imputed_data/'+impute_mtd[i]+'/'+ dataset[j] +'.txt',sep='\t',header=0,index_col=0)
            imputedData=imputedData.values
            imputedData=imputedData.T
            imputedData = scaler.fit_transform(imputedData)
            
        print(type(imputedData),imputedData.shape)
        Datadiff=imputedData-TrueData
        print(type(Datadiff),Datadiff.shape)
        L2_error=np.linalg.norm(Datadiff)
        print(L2_error)
        error_list.append(L2_error)
        print("*******"+impute_mtd[i]+":"+" "+str(L2_error)+"*******")
        print("*******"+impute_mtd[i]+" "+"finish computation"+"*******")
    print(error_list)
    error_list=pd.DataFrame(error_list)
    print("****************************"+dataset[j]+" "+"finish computation"+"****************************")

