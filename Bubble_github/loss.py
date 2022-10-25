
# coding: utf-8

# In[ ]:

import torch
import torch.nn as nn
import torch.nn.functional as F

class SquareRegularizeLoss(nn.Module):
    def __init__(self, p=1):
        super().__init__()
        self.p = p
   ##input:samples*features
   ##regularized loss:1/n*(|1-(x1^2+...+xn^2)|^p)
    def forward(self, input):
        feature_num=input.size(1)
        input = torch.pow(input, 2).sum(dim=1)
        if self.p == 1:
            loss = torch.abs(1-input)
        else:
            loss = torch.pow(1-input, self.p)
        loss = loss.mean()/feature_num
       
        return loss

