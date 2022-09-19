#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 6/8/22 16:16

@Author: Luis Antonio Souto Arias

@Software: PyCharm
"""

import numpy as np

##

# Univariate Kaplan-Meier estimator for right-censored data
def kaplan_meier_rc(x,d):
    mx = x.max()
    Sx = np.ones((mx+1,))
    hx = np.zeros_like(Sx)
    sx = np.zeros_like(Sx)
    
    for i in range(mx+1):
        hx[i] = ((x==i)*d).sum()
        sx[i] = (x>=i).sum()

    Sx[0] = 1.-hx[0]/sx[0]
    for i in range(1,mx+1):
        if sx[i]>0:
            Sx[i]= Sx[i-1]*(1-hx[i]/sx[i])
        else:
            Sx[i] = 0

    return Sx

# Univariate Kaplan-Meier estimator for left-truncated right-censored data
def kaplan_meier_ltrc(x,Tx,d):
    mx = x.max()
    Sx = np.ones((mx+1,))
    hx = np.zeros_like(Sx)
    sx = np.zeros_like(Sx)
    
    for i in range(mx+1):
        hx[i] = ((x==i)*d).sum()
        sx[i] = ((x>=i)*(Tx<=i)).sum()

    if sx[0]>0:
        Sx[0] = 1.-hx[0]/sx[0]
    xm      = min(x)
    Sx[:xm] = 1.

    for i in range(max(xm,1),mx+1):
        if sx[i]>0:
            Sx[i] = Sx[i-1]*(1-hx[i]/sx[i])
        else:
            Sx[i] = 0

    return Sx

