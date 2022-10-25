# Bubble

A single-cell RNA-seq imputation method using an autoencoder constrained by bulk RNA-seq data. 

## Overview

We propose Bubble, an autoencoder-based model that is constrained by the matched bulk RNA-seq data to identify and impute the scRNA-seq data. 

​         First, Bubble identifies the likely non-biological zeros from all zeros based on expression rate and coefficient of variation of genes within cell subpopulation, and only performs imputation on these values.  

​       Next, Bubble is based on an autoencoder framework and simultaneously treats the matched bulk RNA-seq data as prior knowledge and constraint to impute dropout events. 

## Installation

Clone this repository. The Bubble has been implemented in Python3.6 and Pytorch  1.7.0. 

```shell
git clone https://github.com/CSUBioGroup/Bubble.git
cd Bubble_github
```

## Datasets

All real datasets used in our paper can be found in:

- The RNAmix_CEL-seq2, sc_10x_5cl, and three cell mixture datasets from GSE118767, the corresponding bulk RNA-seq samples are available in the NCBI repository with the GEO accession number of GSE86337. 

- The liver, fetal brain, and matched bulk RNA-seq datasets from the Mouse Cell Atlas project and downloaded at https://doi.org/10.5281/zenodo.2585885. 
- The 10x_293t_jurkat dataset comes from the 10x Genomics data platform (https://support.10xgenomics.com/single-cell-geneexpression/datasets/1.1.0/293t, https://support.10xgenomics.com/singlecell-gene-expression/datasets/1.1.0/jurkat). The matched bulk RNA-seq dataset can be accessed with the GEO number GSE129240. 

- The HCA_10x_tissue dataset from sample MantonBM6 measured using 10x Genomics (https://data.humancellatlas.org/explore/projects/cc95ff89-2e68-4a08-a234-480eca21ce79). The matched bulk RNA-seq samples are available at GEO under accession code GSE74246. 

- The Tabula Sapiens atlas is downloaded from https://cellxgene.cziscience.com/e/5a11f879-d1ef-458a-910c-9b0bdfca5ebf.cxg/.

## Input data and format

- sc_dataset: the normalized (logarithmized) scRNA-seq data  (cell by gene);
- bulkdata: the normalized (logarithmized)  matched bulk RNA-seq data  (sample by gene);
- identify_Data: the identification dataset, you can get the identity dataset by running identification.py (The position corresponding to the dropouts will be marked as -1) (gene by cell).

The detailed format is as follows:

```python
import torch
from torch.utils.data import Dataset
import os
import pandas as pd
import numpy as np
from os.path import join

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
```



## Output data

Reconstructed scRNA-seq data (you can get the decoded result in ./results/datasets/exp_id/decoding.txt').

***Note**:* you need to correct the reconstructed dataset to get the final imputation result (ie: keep non-zero expression values and the biological zeros unchanged).

## For calling Bubble by Bash commands

**Step1:** identify non-biological zeros from all zeros to obtain the identification dataset;

- run the identification.py

**Step2:** training and inference, input scRNA-seq dataset, identification dataset, and matched bulk RNA-seq dataset;

- Command line: 


```shell
python3 train.py --exp_id experiment id --datasets dataset_name --eps epochs of training --bs batch_size of training
python3 infer.py --exp_id experiment id --datasets dataset_name
```

**Step3:** please make sure to keep true expression values unchanged.

## Look for more usage of Bubble via

```
python3 train.py --help
```

```shell
usage: train.py [-h] --exp_id EXP_ID --datasets DATASETS [--lr LR] [--eps EPS]
                [--bs BS] [--gpu_id GPU_ID]

optional arguments:
  -h, --help           show this help message and exit
  --exp_id EXP_ID      experiment id
  --datasets DATASETS  dataset_name
  --lr LR              learning rate
  --eps EPS            reset training epochs of training
  --bs BS              batch_size for training only
  --gpu_id GPU_ID      which gpu to use
```

