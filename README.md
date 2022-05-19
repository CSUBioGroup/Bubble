# Bubble

A single-cell RNA-seq imputation method using an autoencoder constrained by bulk RNA-seq data. 

## Overview

We propose Bubble, an autoencoder-based model that is constrained by the matched bulk RNA-seq data to identify and impute the scRNA-seq data. 

​         First, Bubble identifies the likely non-biological zeros from all zeros based on expression rate and coefficient of variation of genes within cell subpopulation, and only performs imputation on these values.  

​       Next, Bubble is based on an autoencoder framework and simultaneously treats the matched bulk RNA-seq data as prior knowledge and constraint to impute dropout events. 

## Installation

Clone this repository.

```shell
git clone https://github.com/CSUBioGroup/Bubble.git
cd Bubble
```

### Usage

The Bubble has been implemented in Python3.6 and Pytorch  1.7.0. 

Step1: identify non-biological zeros from all zeros to obtain the identification dataset;

- run the identification.py

Step2: training and inference, input scRNA-seq dataset, identification dataset, and matched bulk RNA-seq dataset;

Command line: 

```shell
python3 train.py [--exp_id experiment id] [--datasets dataset_name] [--eps epochs of training] [--bs batch_size of training] 
python3 infer.py [--exp_id experiment id] [--datasets dataset_name]
```

Step3: please make sure to keep true expression values unchanged.
