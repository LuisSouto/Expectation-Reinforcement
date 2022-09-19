#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 6/9/22 10:40

@Author: Luis Antonio Souto Arias

@Software: PyCharm
"""

## Load modules

import numpy as np
import matplotlib.pyplot as plt

## Observations and theoretical curve (Poisson)
output_dir   = '../figures/'
maxY = 85
mY   = 65
errT = np.array([4.18741e-5,1.72211e-5,6.06898e-6,2.24272e-6,7.6358e-7,1.8809e-7,3.46329e-08,6.50852e-09,9.48194e-10])
n    = errT.size
idx  = 2
v    = 85+np.arange(1+idx,n+1)
errt = errT[idx]*np.exp(-v*(v-v[0])/mY)

## Plot results
plt.figure()
plt.plot(v,errT[idx:],'r*')
plt.plot(v,errt,'b--')
plt.yscale("log")
plt.legend(("Error",r"$\exp(-y_M^2/\lambda_Y)$"),fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel(r"$y_M$",fontsize=16)
plt.ylabel(r"Error",fontsize=16)
plt.xticks(v)
plt.savefig(output_dir+'poisson_error.eps',format='eps')
plt.show()

##

