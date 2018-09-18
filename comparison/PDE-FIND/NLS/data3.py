# -*- coding: utf-8 -*-
"""

@author: Ye Yuan

utilizing PDE-FIND method to analyze the NLS equation:          
             u{t} = 0.5iu{xx} + iu|u|^2

"""

import numpy as np
import scipy.io as sio
from PDE_FIND import *

data = sio.loadmat('data3.mat')
theta = data['theta3']
y = data['y3']

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
    
minerr = min(err[0,:])   