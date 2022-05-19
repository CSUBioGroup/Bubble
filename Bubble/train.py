# coding: utf-8
import psutil
import datetime

import pandas as pd
import argparse
from os.path import join
import os
import time
from tqdm import tqdm
import torch
import torch.nn as nn
import torch.optim as optim
import torch.utils.data.dataloader as dataloader
from torch.utils.data import RandomSampler,BatchSampler
from tensorboardX import SummaryWriter
from config import Config
from utils import *
from datasets import *
from model import *
from loss import *

config = Config()
# data_root = '/home/csq/AutoImpute/dataset'

# exp_id = 'v1.0'

# middle_layer_size = [256, 128, 256]  # [n_input_features] + middle_layer_size + [n_output_features]
start_time = time.time()
start_memory=psutil.Process(os.getpid()).memory_info().rss / 1024 
if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        '--exp_id',
        required=True,
        help='experiment id')
    argparser.add_argument(
        '--datasets',
        required=True,
        help='dataset_name')
    argparser.add_argument(
        '--lr',
        default=1e-4,
        type=float,
        help='learning rate')
    argparser.add_argument(
        '--eps',
        default=80,
        type=int,
        help='reset training epochs of training')
    argparser.add_argument(
        '--bs',
        default=64,
        type=int,
        help='batch_size for training only')
    argparser.add_argument(
        '--gpu_id',
        default='0',
        help='which gpu to use')

    args = argparser.parse_args()

    os.environ['CUDA_VISIBLE_DEVICES'] = args.gpu_id
    logs_dir = 'logs/%s/%s' % (args.datasets, args.exp_id)
    create_dirs([logs_dir])
    
    dataset = SingleCell(
        config.data_root,
        args.datasets,
    )
    dataset=dataset.load_data()
    #load identify_Data
    identify_Data=pd.read_csv("./identify_data/"+args.datasets+"_identify.txt",sep='\t',header=None,index_col=None)
    identify_Data=identify_Data.T
    identify_Data=identify_Data.values
    identify_Data=torch.Tensor(identify_Data)
    identify_Data=identify_Data.to('cuda')
    
    #load bulk data
    bulkdata = np.loadtxt(join(config.data_root, args.datasets+'_bulk'+'.txt'))
    bulkdata=torch.Tensor(bulkdata)
    bulkdata=bulkdata.to('cuda')

    #load normalized raw data
    scdata = torch.Tensor(dataset)
    scdata = scdata.to('cuda')
    
    data_loader = dataloader.DataLoader(
          dataset = scdata,
          sampler = BatchSampler(RandomSampler(dataset),batch_size=args.bs, drop_last=False)
    )
    
    data_iter = iter(data_loader)  
    sample_batch = next(data_iter) 

    in_dim = sample_batch.size(-1) 
    
    ## model configs
    model_layer_sizes = [in_dim] + config.middle_layer_size + [in_dim]

    #1. Construct the network structure, 2. Initialize parameters, define loss function, define update rules
    autoencoder = AutoEncoder(model_layer_sizes)
    
    if torch.cuda.is_available():
        autoencoder.to('cuda')  

    NonzeroLoss = nn.MSELoss()
    NonzeroLoss.to('cuda')
    RegularizeLoss = SquareRegularizeLoss(p=config.p)
    RegularizeLoss.to('cuda')
    bulkloss=nn.MSELoss()
    bulkloss.to('cuda')
    
    optimizer = optim.Adam(autoencoder.parameters(),lr=args.lr)

    best_loss = 1e10
    
    Loss_total=[]
    Loss_nonzero=[]
    Loss_regularize=[]
    Loss_bulk=[]
    for epoch in range(args.eps):
        epoch_loss,epoch_nonzero,epoch_regularize,epoch_bulkloss = 0,0,0,0
        print('=====epoch start=====',epoch)
        iter_num=0
        for sample_idx in data_loader.sampler:
            iter_num = iter_num+1
            index=sample_idx
            index=torch.Tensor(index)
            index=index.type(torch.LongTensor)
            index=index.to('cuda')
            #print('====index=====',index.size(),index.size(0),index.get_device())
            batch=scdata[index]
            #print('====batch=====',batch.size(),batch.get_device())
            if torch.cuda.is_available():
                batch = batch.to('cuda')
                #print('=====batch=====',batch.size(),batch.get_device())
                
            enc_out, dec_out = autoencoder(batch)
            optimizer.zero_grad()
            #define true signal loss 
            identify_batch=identify_Data[index]
            mask=torch.where(identify_batch==-1,torch.zeros_like(batch),torch.ones_like(batch)) 
            mask=mask.to('cuda')
            loss_nonzero = NonzeroLoss(dec_out.mul(mask), batch)
            #define regularization loss 
            loss_regularize = RegularizeLoss(enc_out)
            #define bulk loss 
            scdata_copy=scdata.clone()  
            scdata_copy[index]=dec_out    
            scdata_copy=torch.mean(scdata_copy,dim=0) 
            loss_bulk = bulkloss(scdata_copy,bulkdata)
            #define total loss
            loss = loss_nonzero + loss_bulk + loss_regularize 

            #loss backward
            loss.backward()

            #Initialize parameters
            optimizer.step()

            epoch_nonzero += loss_nonzero.data.cpu().numpy()*batch.size(0)
            epoch_bulkloss += loss_bulk.data.cpu().numpy()
            epoch_regularize += loss_regularize.data.cpu().numpy()*batch.size(0)
            epoch_loss += epoch_nonzero+epoch_bulkloss+epoch_regularize
        epoch_nonzero /= len(scdata)
        epoch_regularize /= len(scdata)
        epoch_bulkloss /= iter_num
        epoch_loss = epoch_nonzero+epoch_bulkloss+epoch_regularize
        
        Loss_total.append(epoch_loss)
        Loss_nonzero.append(epoch_nonzero)
        Loss_bulk.append(epoch_bulkloss)
        Loss_regularize.append(epoch_regularize)
 
        if epoch_loss<best_loss:
            best_loss = epoch_loss
            torch.save(autoencoder.state_dict(), join(logs_dir, 'model_{:03d}_loss_{:.3f}.pth'.format(epoch, best_loss)))

    
    print('=====Loss total=====',Loss_total.index(min(Loss_total)),min(Loss_total))
    print('=====Loss nonzero=====',Loss_nonzero.index(min(Loss_nonzero)),min(Loss_nonzero))
    print('=====Loss bulk=====',Loss_bulk.index(min(Loss_bulk)),min(Loss_bulk))
    print('=====Loss regularize=====',Loss_regularize.index(min(Loss_regularize)),min(Loss_regularize))
    
    end_memory=psutil.Process(os.getpid()).memory_info().rss/1024
    end_time = time.time()
    print('=====use memory=====',end_memory-start_memory,'kb')
    print('=====use time=====', end_time-start_time,'second')
    
    logs_dir = '/home/csq/AutoImpute/csq-autoimpute-main/loss_mousedata/%s/%s' % (args.datasets, args.exp_id)
    create_dirs([logs_dir])

    Loss_total=pd.DataFrame(Loss_total)
    Loss_total=Loss_total.T
    Loss_nonzero=pd.DataFrame(Loss_nonzero)
    Loss_nonzero=Loss_nonzero.T
    Loss_regularize=pd.DataFrame(Loss_regularize)
    Loss_regularize=Loss_regularize.T
    Loss_bulk=pd.DataFrame(Loss_bulk)
    Loss_bulk=Loss_bulk.T

    Loss_total.to_csv(logs_dir+"/Loss_total.txt",sep='\t',header=0,index=0)
    Loss_nonzero.to_csv(logs_dir+"/Loss_nonzero.txt",sep='\t',header=0,index=0)
    Loss_regularize.to_csv(logs_dir+"/Loss_regularize.txt",sep='\t',header=0,index=0)
    Loss_bulk.to_csv(logs_dir+"/Loss_bulk.txt",sep='\t',header=0,index=0)
    
    print('training finished')   