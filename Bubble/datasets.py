import torch
from torch.utils.data import Dataset

import os
import numpy as np
from os.path import join

from utils import preprocess

class SingleCell(Dataset):
    def __init__(self, data_root, dataset_name):
        self.data_root = data_root
        self.dataset_name = dataset_name
        #self.load_data()

    def __len__(self):
        return len(self.X)

    def load_data(self):
        # customized
        self.X = np.loadtxt(join(self.data_root, self.dataset_name+'.txt'))
        print(self.X.shape,type(self.X),np.min(self.X),np.max(self.X))
        return self.X

    def __getitem__(self, i):
        return self.X[i] 
