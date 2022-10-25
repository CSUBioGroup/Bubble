
# coding: utf-8

# In[ ]:

import numpy as np
import os
from os.path import join
from sklearn.preprocessing import MinMaxScaler

def create_dirs(dirs):
    for _dir in dirs:
        os.makedirs(_dir, exist_ok=True)

def preprocess(x):
    scaler = MinMaxScaler()
    y = scaler.fit_transform(x.T).T
    return y

def normalize(x):
    pass
#ckpts = find_last(log_dir)
def find_last(log_dir):
    ckpts = next(os.walk(log_dir))[2]
    ckpts = sorted(filter(lambda x:x.endswith('.pth'), ckpts))

    if not ckpts:
        import errno
        raise FileNotFoundError(
            errno.ENOENT,
            "Could not find model ckpts"
            )

    ckpt = join(log_dir, ckpts[-1])
    return ckpt

