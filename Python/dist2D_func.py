#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 6/8/22 18:08

@Author: Luis Antonio Souto Arias

@Software: PyCharm
"""

import numpy as np

def convolution2D(pA,pB,pC):
    m   = pA.size
    PXY = np.zeros((m,m))
    for i in range(m):
        for j in range(i):
            idx = j+1
            PXY[i,j] = np.sum(pA[:idx]*pB[i:i-idx:-1]*pC[j::-1])
        PXY[i,i] = np.sum(pA[:i+1]*pB[i::-1]*pC[i::-1])
        for j in range(i+1,m):
            idx = i+1
            PXY[i,j] = np.sum(pA[:idx]*pB[i::-1]*pC[j:j-idx:-1])

    return PXY

def pdf_to_surv2D(PXY):
    m   = PXY.shape[0]
    SXY = np.zeros_like(PXY)
    for i in range(m):
      for j in range(m):
        SXY[i,j]  = PXY[i+1:,j+1:].sum()

    return SXY

def pdf_to_mean(P):
    m = P.size
    v = np.arange(m)
    return np.sum(v*P)

def pdf_to_var(P):
    m = P.size
    v = np.arange(m)
    return np.sum(v**2*P)-pdf_to_mean(P)**2

def pdf_to_cov(P):
    m = P.shape[0]
    y = np.arange(m)
    x = np.arange(m)[:,np.newaxis]
    return np.sum(x*y*P)-pdf_to_mean(P.sum(1))*pdf_to_mean(P.sum(0))

def sample_from_cdf(F,n):
    U = np.random.random((n,1))
    return np.sum(U>=F,1).squeeze()
