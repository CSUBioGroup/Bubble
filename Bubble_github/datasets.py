
# coding: utf-8

# In[ ]:


import torch
from torch.utils.data import Dataset

import os
import pandas as pd
import numpy as np
from os.path import join

from utils import preprocess

class SingleCell(Dataset):
    def __init__(self, data_root, dataset_name):
        self.data_root = data_root
        self.dataset_name = dataset_name 

    def __len__(self):
        return len(self.X)

    def load_data(self):
        # customized
        self.X = np.loadtxt(join(self.data_root, self.dataset_name+'.txt'))
        #print(self.X.shape)
        return self.X
    def load_identify_data(self):
        identify_Data=pd.read_csv(join(self.data_root, self.dataset_name+'_identify.txt'),sep=',',header=0,index_col=0)
        identify_Data=identify_Data.T
        identify_Data=identify_Data.values
        identify_Data=torch.Tensor(identify_Data)
        identify_Data=identify_Data.to('cuda')
        #print('=====identify_Data=====',identify_Data.size(),identify_Data.get_device())
        return identify_Data
    def load_bulk_data(self):
        bulkdata = np.loadtxt(join(self.data_root, self.dataset_name+'_bulk'+'.txt'))
        bulkdata=torch.Tensor(bulkdata)
        bulkdata=bulkdata.to('cuda')
        #print('=====bulkdata device=====',bulkdata.size(),bulkdata.get_device())
        return bulkdata
    def __getitem__(self, i):
        return self.X[i] 

