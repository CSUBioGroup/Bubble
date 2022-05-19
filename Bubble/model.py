import torch
import torch.nn as nn
import os
import time

class AutoEncoder(nn.Module):
    def __init__(self, layer_size=[]):  
        super(AutoEncoder, self).__init__()
        self.layer_size = layer_size
        self.n_layers   = len(layer_size)//2
        print(self.layer_size,self.n_layers)
        self.build()

    def build(self):
        
        self.encoders = nn.ModuleList()
        self.decoders = nn.ModuleList()
        
        n_layers = self.n_layers#2
        layer_size = self.layer_size#[12548, 256, 128, 256, 12548]
        
        #range(2):0,1
        for i in range(n_layers):
            
            self.encoders.append(nn.Sequential(
                    nn.Linear(in_features=layer_size[i], out_features=layer_size[i+1]),
                    nn.ReLU()
                )
            )
            
            self.decoders.append(nn.Sequential(
                    nn.Linear(in_features=layer_size[i+1], out_features=layer_size[i]),
                    nn.ReLU()
                ) 
            )
    def forward(self, input):
        enc_out = self.encoders[0](input)
        # encoder
        for i in range(1, self.n_layers):
            enc_out = self.encoders[i](enc_out)
        # decoder
        dec_out = self.decoders[-1](enc_out)
        for i in range(1, self.n_layers):
            dec_out = self.decoders[-(i+1)](dec_out)
        return enc_out, dec_out
    
    
    
    
    
    
    
