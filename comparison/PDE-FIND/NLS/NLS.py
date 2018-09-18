# -*- coding: utf-8 -*-
"""

@author: Ye Yuan

utilizing PDE-FIND method to analyze the NLS equation:          
             u{t} = 0.5iu{xx} + iu|u|^2

"""

import numpy as np
import scipy.io as sio
from PDE_FIND import *

data1 = sio.loadmat('data1.mat')
theta1 = data1['theta1']
y1 = data1['y1']

data2 = sio.loadmat('data2.mat')
theta2 = data2['theta2']
y2 = data2['y2']

data3 = sio.loadmat('data3.mat')
theta3 = data3['theta3']
y3 = data3['y3']

w1 = TrainSTRidge(theta1,y1,10**-5,500)
w2 = TrainSTRidge(theta2,y2,10**-5,500)
w3 = TrainSTRidge(theta3,y3,10**-5,150)

print("\ndata1:\nu{t} = (%.4f +%.4f)u|u||u| +(%.4f+%.4f)u{xx}" %  (w1[11].real,w1[11].imag,w1[8].real,w1[8].imag))
err1 = abs(np.array([(1 -  w1[8].imag)*100, (.5 - w1[11].imag)*100/0.5]))
print("error1: mean: %.4f%% std: %.4f%% " % (np.mean(err1),np.std(err1)))

#failed
print("\ndata2:\nu{t} = (%.4f +%.4f)u|u||u| +(%.4f+%.4f)u{xx}" %  (w2[11].real,w2[11].imag,w2[8].real,w2[8].imag))
err2 = abs(np.array([(1 -  w2[8].imag)*100, (.5 - w2[11].imag)*100/0.5]))
print("error2: mean: %.4f%% std: %.4f%% " % (np.mean(err2),np.std(err2)))

#failed
print("\ndata3:\nu{t} = (%.4f +%.4f)u|u||u| +(%.4f+%.4f)u{xx}" %  (w3[11].real,w3[11].imag,w3[8].real,w3[8].imag))
err3 = abs(np.array([(1 -  w3[8].imag)*100, (.5 - w3[11].imag)*100/0.5]))
print("error3: mean: %.4f%% std: %.4f%% " % (np.mean(err3),np.std(err3)))

#err = (err1 + err3)/2
#print("\nerror: mean: %.4f%% std: %.4f%% " % (np.mean(err),np.std(err)))

