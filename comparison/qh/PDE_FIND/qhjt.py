# -*- coding: utf-8 -*-
"""
@author: Ye Yuan
utilizing PDE-FIND method to analyze the QH equation:            
          u{t} = 0.5iu{xx} - i(x^2/2)u
"""

import numpy as np
import scipy.io as sio
from PDE_FIND import *
import warnings
warnings.filterwarnings('ignore')

data = sio.loadmat('data1xt500.mat')
theta = data['theta1']
y = data['y1']

wall = np.ones((theta.shape[1],1))
name = 'data1-xt500-grid2'
d_tol_all = []
lam_all = []
for lam in np.arange(start=1e-4,step=1e-4,stop=1e-2):
    for d_tol_i in np.arange(9):
        for d_tol_j in np.arange(start=1,step=0.5,stop=9):

            d_tol = 10**(-5)*10**(d_tol_i)*d_tol_j
            d_tol_all.append(d_tol)
            lam_all.append(lam)
            w = TrainSTRidge(theta,y,lam,d_tol)
            print lam,d_tol
            wall = np.hstack([wall,w])
            print w

        sio.savemat(name,{'y':y,'theta':theta,'wall':wall,'d_tol_all':d_tol_all,'lam_all':lam_all})
sio.savemat(name,{'y':y,'theta':theta,'wall':wall,'d_tol_all':d_tol_all,'lam_all':lam_all})
