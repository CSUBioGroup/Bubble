
# coding: utf-8

# In[ ]:


import numpy as np

class Config(object):
    # data directory root
    data_root = '/home/csq/AutoImpute/csq-autoimpute-main/Bubble_github/testdata/'

    # model configs
    middle_layer_size = [256, 128, 256]

    # regularized loss, (1-(x1^2+...+xn^2))^p
    p = 2

    # format for saving encoded results
    formt = 'txt'   # 'npy'

