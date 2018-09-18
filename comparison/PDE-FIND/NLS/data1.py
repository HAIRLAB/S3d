# -*- coding: utf-8 -*-
"""

@author: Ye Yuan

utilizing PDE-FIND method to analyze the NLS equation:          
             u{t} = 0.5iu{xx} + iu|u|^2

"""

import numpy as np
import scipy.io as sio
from PDE_FIND import *

data = sio.loadmat('data1.mat')
theta = data['theta1']
y = data['y1']

d_tol = np.zeros((144,1))
for i in range(2,20):
    for j in range(8):
        d_tol[i-2+j*18] = 0.5*i*10**(j-4)
        
w_all = np.zeros((40,144),complex)     
err = np.zeros((2,144))    
for i in range(144):
    w = TrainSTRidge(theta,y,10**-5,d_tol[i])
    w_all[:,i] = w.ravel()
    err[0,i] = np.mean(abs(np.array([(1 -  w[8].imag)*100, (.5 - w[11].imag)*100/0.5])))
    err[1,i] = np.std(abs(np.array([(1 -  w[8].imag)*100, (.5 - w[11].imag)*100/0.5])))
    
    
d_tol2 = np.zeros((3333,1))
for i in range(3333):
    d_tol2[i] = 1+(i+1)*3

w_all2 = np.zeros((40,3333),complex)
err2 = np.zeros((2,3333))
for i in range(3333):
    w = TrainSTRidge(theta,y,10**-5,d_tol2[i])
    w_all2[:,i] = w.ravel()
    err2[0,i] = np.mean(abs(np.array([(1 -  w[8].imag)*100, (.5 - w[11].imag)*100/0.5])))
    err2[1,i] = np.std(abs(np.array([(1 -  w[8].imag)*100, (.5 - w[11].imag)*100/0.5])))
    
minerr1 = min(err[0,:])
minerr2 = min(err2[0,:])
