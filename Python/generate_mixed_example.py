#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 4/8/22 11:49

@Author: Luis Antonio Souto Arias

@Software: PyCharm
"""

import numpy as np
import pandas as pd
from scipy.stats import nhypergeom, betabinom, randint

## Write configuration in info file

ltx   = 40
la    = 40
lb    = 40
lc    = 40
sc    = 6
n     = 10000
eps0  = -5
leps  = 7
lcens = 2

distA     = 'Uniform '+str(la)
distB     = 'Beta-binomial '+str(lb)
distC     = 'Negative hypergeometric '+str(lc)
disttrunc = 'Poisson '+str(ltx)
disteps   = 'Poisson '+str(leps)+' '+str(eps0)
distcens  = 'Poisson '+str(lcens)

# Print info about the true solution
fname = '../data/mixed_ltrc'
finfo = open(fname+'_info.dat','w')
print('Distribution of A: '+distA,file=finfo)
print('Distribution of B: '+distB,file=finfo)
print('Distribution of C: '+distC,file=finfo)
print('Distribution of truncation of X: '+disttrunc,file=finfo)
print('Distribution of epsilon: '+disteps,file=finfo)
print('Distribution of censoring: '+distcens,file=finfo)

pA = randint(0,la+1)
pB = betabinom(lb,3.8,4.25)
pC = nhypergeom(lc,np.floor(0.8*lc),2)

## Write data in output file
fout = open(fname+'.dat','w')
print(n,file=fout)

X   = np.zeros((n,),dtype='int32')
Y   = np.zeros_like(X,dtype='int32')
TX  = np.zeros_like(X,dtype='int32')
TY  = np.zeros_like(X,dtype='int32')
eps = np.zeros_like(X,dtype='int32')
d   = np.ones_like(X,dtype='int32')
e   = np.ones_like(X,dtype='int32')

np.random.seed(10)

N = 100000  # Maximum memory allocation. More may be needed depending on the truncation.
A = pA.rvs(N)
B = pB.rvs(N)
C = pC.rvs(N)

iter  = 0
trunc = True
for i in range(n):
    if trunc:
        ai = 0.
        bi = 0.
        ci = 0.
        txi = ai+bi+1
        tyi = ai+ci+1
        while((txi>(ai+bi))|(tyi>(ai+ci))):
            ai   = A[iter]
            bi   = B[iter]
            ci   = C[iter]
            txi  = np.random.poisson(ltx)
            epsi = np.random.poisson(leps)+eps0
            tyi  = txi+epsi
            iter += 1
    else:
        ai   = A[i]
        bi   = B[i]
        ci   = C[i]
        txi  = np.random.poisson(ltx)
        epsi = np.random.poisson(leps)+eps0
        tyi  = txi+epsi
    cs     = np.random.poisson(lcens)
    X[i]   = min(ai+bi,txi+cs)
    Y[i]   = min(ai+ci,tyi+cs)
    TX[i]  = txi
    TY[i]  = tyi
    d[i]   = (X[i]==(ai+bi))
    e[i]   = (Y[i]==(ai+ci))
    eps[i] = epsi
    iter  += 1
    print(X[i],Y[i],TX[i],TY[i],d[i],e[i],eps[i],file=fout)

print('Number of iterations',iter)
finfo.close()
fout.close()

##
##

