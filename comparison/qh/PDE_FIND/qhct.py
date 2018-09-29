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
name = 'data1-xt500-grid1'
d_tol_all = []
lam_all = []
for lam in [10**(-8),10**(-7), 10**(-6),10**(-5),10**(-4),10**(-3),10**(-2),10**(-1)]:
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


# d_tol = np.zeros((144,1))
# for i in range(2,20):
#     for j in range(8):
#         d_tol[i-2+j*18] = 0.5*i*10**(j-4)
#
# w_all = np.zeros((40,144),complex)
# err = np.zeros((2,144))
# for i in range(144):
#     w = TrainSTRidge(theta,y,10**-5,d_tol[i])
#     w_all[:,i] = w.ravel()
#     err[0,i] = np.mean(abs(np.array([(1 -  w[8].imag)*100, (.5 - w[11].imag)*100/0.5])))
#     err[1,i] = np.std(abs(np.array([(1 -  w[8].imag)*100, (.5 - w[11].imag)*100/0.5])))
#
#
# d_tol2 = np.zeros((3333,1))
# for i in range(3333):
#     d_tol2[i] = 1+(i+1)*3
#
# w_all2 = np.zeros((40,3333),complex)
# err2 = np.zeros((2,3333))
# for i in range(3333):
#     w = TrainSTRidge(theta,y,10**-5,d_tol2[i])
#     w_all2[:,i] = w.ravel()
#     err2[0,i] = np.mean(abs(np.array([(1 -  w[8].imag)*100, (.5 - w[11].imag)*100/0.5])))
#     err2[1,i] = np.std(abs(np.array([(1 -  w[8].imag)*100, (.5 - w[11].imag)*100/0.5])))
#
# minerr1 = min(err[0,:])
# minerr2 = min(err2[0,:])
