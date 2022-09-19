#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 4/8/22 11:29

@Author: Luis Antonio Souto Arias

@Software: PyCharm
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import poisson, nhypergeom, betabinom, randint
from scipy.stats import permutation_test, ttest_ind

import dist2D_func as distf
from kaplan_meier_func import kaplan_meier_ltrc
from frees_copula import gompertz_pdf,frank_copula_dudv

## Probability of truncation
def prob_trunc(pXY,pT,pE,epsm):
    m      = pT.size
    FE     = np.concatenate((np.cumsum(pE),np.ones((epsm,))))
    pTE    = np.zeros_like(pXY)
    for x in range(m):
        for y in range(x-epsm+1):
            pTE[x,y] = np.dot(pT[:y+epsm+1],FE[y+epsm::-1])
        for y in range(x-epsm+1,m):
            pTE[x,y] = np.dot(pT[:x+1],FE[y+epsm:y+epsm-x-1:-1])

    ptrunc = 1.-(pXY*pTE).sum()
    return ptrunc

## Read input data

filename     = 'poisson'
nr           = ''
filename_res = filename+nr
fdata        = '../data/'+filename+'_ltrc.dat'
n            = pd.read_csv(fdata,nrows=1,header=None).to_numpy().squeeze()
data         = pd.read_csv(fdata,skiprows=1,header=None,sep=None,engine='python').to_numpy()
output_dir   = '../figures/'

if filename=='UOFM':
    x,y,Tx,Ty,d,e = data.T
    eps = Ty-Tx
else:
    x,y,Tx,Ty,d,e,eps = data.T

minT,maxT = np.min(Tx),np.max(Tx)
minE,maxE = np.min(eps),np.max(eps)
epsm      = -minE

## Read posterior data
fpostl    = '../results/'+filename_res+'_ltrc_res_l.dat'
mxy       = pd.read_csv(fpostl, nrows=1, header=None).to_numpy().squeeze()
post_dist = pd.read_csv(fpostl, skiprows=1, header=None, sep=None, engine='python').to_numpy()

pAl,pBl,pCl,pTl,pEl = post_dist.T

fposth    = '../results/'+filename_res+'_ltrc_res_h.dat'
mxy       = pd.read_csv(fposth, nrows=1, header=None).to_numpy().squeeze()
post_dist = pd.read_csv(fposth, skiprows=1, header=None, sep=None, engine='python').to_numpy()

pAh,pBh,pCh,pTh,pEh = post_dist.T

## Compute prior and reference distribution
xd   = x[d==1]
ye   = y[e==1]
xu   = x[d*e==1]
yu   = y[d*e==1]
minX,maxX = min(xd),max(xd)
minY,maxY = min(ye),max(ye)

# Empirical Kaplan-Meier
nx  = x.max()
ny  = y.max()
vx  = np.arange(nx+1)
vy  = np.arange(ny+1)
Sx  = kaplan_meier_ltrc(x,Tx,d)
Sy  = kaplan_meier_ltrc(y,Ty,e)
Px  = np.zeros_like(Sx)
Py  = np.zeros_like(Sy)
Px[0]  = 1.-Sx[0]
Px[1:] = Sx[0:-1]-Sx[1:]
Px    /= np.sum(Px)
Py[0]  = 1.-Sy[0]
Py[1:] = Sy[0:-1]-Sy[1:]
Py    /= np.sum(Py)

# Prior distribution
if nr=='2':
    lap,lbp,lcp,ltp,lep = 40,40,45,70,40
else:
    lap,lbp,lcp,ltp,lep = 20,20,20,50,10
v      = np.arange(mxy)
priorA = poisson.pmf(v,lap)
priorB = poisson.pmf(v,lbp)
priorC = poisson.pmf(v,lcp)
priorT = poisson.pmf(v,ltp)
priorE = poisson.pmf(v,lep)

# Analytical solution
if filename=='mixed':
    la,lb,lc,lt,le = 40,40,40,40,7
    pA_an = randint.pmf(v,0,la+1)
    pB_an = betabinom.pmf(v,lb,3.8,4.25)
    pC_an = nhypergeom.pmf(v,lc,np.floor(0.8*lc),2)
else:
    la,lb,lc,lt,le = 40,20,25,70,7
    pA_an = poisson.pmf(v,la)
    pB_an = poisson.pmf(v,lb)
    pC_an = poisson.pmf(v,lc)
pT_an = poisson.pmf(v,lt)
pE_an = poisson.pmf(v,le)
SA_an = 1.-np.cumsum(pA_an)
SB_an = 1.-np.cumsum(pB_an)
SC_an = 1.-np.cumsum(pC_an)

# Frank copula model
mx,sx,my,sy,alpha = 84.809,9.926,87.575,7.792,-4.081
PXc = gompertz_pdf(v,mx,sx)
PYc = gompertz_pdf(v,my,sy)

PXY     = distf.convolution2D(pA_an,pB_an,pC_an)
postlXY = distf.convolution2D(pAl,pBl,pCl)
posthXY = distf.convolution2D(pAh,pBh,pCh)
priorXY = distf.convolution2D(priorA,priorB,priorC)
PXYc    = PXc[:,np.newaxis]*PYc*frank_copula_dudv(np.cumsum(PXc),np.cumsum(PYc),alpha)

# Condition on the minimum and maximum observed values
PXY[:minX]     = 0.
PXY[:,:minY]   = 0.
PXY[maxX+1:]   = 0.
PXY[:,maxY+1:] = 0.
PXY           /= PXY.sum()
pT_an[:minT]   = 0.
pT_an[maxT:]   = 0.
pE_an[:minE+5] = 0.
pE_an[maxE+5:] = 0.
pE_an         /= pE_an.sum()
pT_an         /= pT_an.sum()
ST_an          = 1.-np.cumsum(pT_an)
SE_an          = 1.-np.cumsum(pE_an)

priorXY[:minX]     = 0.
priorXY[:,:minY]   = 0.
priorXY[maxX+1:]   = 0.
priorXY[:,maxY+1:] = 0.
priorXY           /= priorXY.sum()
priorT[:minT]      = 0.
priorT[maxT:]      = 0.
priorE[maxE-minE:] = 0.
priorE            /= priorE.sum()
priorT            /= priorT.sum()
STp                = 1.-np.cumsum(priorT)
SEp                = 1.-np.cumsum(priorE)

postlXY[:minX]     = 0.
postlXY[:,:minY]   = 0.
postlXY[maxX+1:]   = 0.
postlXY[:,maxY+1:] = 0.
postlXY           /= postlXY.sum()

posthXY[:minX]     = 0.
posthXY[:,:minY]   = 0.
posthXY[maxX+1:]   = 0.
posthXY[:,maxY+1:] = 0.
posthXY           /= posthXY.sum()

PXYc[:minX]     = 0.
PXYc[:,:minY]   = 0.
PXYc[maxX+1:]   = 0.
PXYc[:,maxY+1:] = 0.
PXYc           /= PXYc.sum()

# Survival functions and marginals
PX     = PXY.sum(1)
PY     = PXY.sum(0)
PXp    = priorXY.sum(1)
PYp    = priorXY.sum(0)
PXpol  = postlXY.sum(1)
PYpol  = postlXY.sum(0)
PXpoh  = posthXY.sum(1)
PYpoh  = posthXY.sum(0)
PXc    = PXYc.sum(1)
PYc    = PXYc.sum(0)
SX     = 1.-np.cumsum(PX)
SY     = 1.-np.cumsum(PY)
SXp    = 1.-np.cumsum(PXp)
SYp    = 1.-np.cumsum(PYp)
SXpol  = 1.-np.cumsum(PXpol)
SYpol  = 1.-np.cumsum(PYpol)
SXpoh  = 1.-np.cumsum(PXpoh)
SYpoh  = 1.-np.cumsum(PYpoh)
SXc    = 1.-np.cumsum(PXc)
SYc    = 1.-np.cumsum(PYc)

# Compute moments and correlation
mX = distf.pdf_to_mean(PX)
sX = distf.pdf_to_var(PX)
mY = distf.pdf_to_mean(PY)
sY = distf.pdf_to_var(PY)
corrXY = distf.pdf_to_cov(PXY)/np.sqrt(sX*sY)

mXp = distf.pdf_to_mean(PXp)
sXp = distf.pdf_to_var(PXp)
mYp = distf.pdf_to_mean(PYp)
sYp = distf.pdf_to_var(PYp)
corrXYp = distf.pdf_to_cov(priorXY)/np.sqrt(sXp*sYp)

mXpol = distf.pdf_to_mean(PXpol)
sXpol = distf.pdf_to_var(PXpol)
mYpol = distf.pdf_to_mean(PYpol)
sYpol = distf.pdf_to_var(PYpol)
corrXYpol = distf.pdf_to_cov(postlXY)/np.sqrt(sXpol*sYpol)

mXpoh = distf.pdf_to_mean(PXpoh)
sXpoh = distf.pdf_to_var(PXpoh)
mYpoh = distf.pdf_to_mean(PYpoh)
sYpoh = distf.pdf_to_var(PYpoh)
corrXYpoh = distf.pdf_to_cov(posthXY)/np.sqrt(sXpoh*sYpoh)

mXkm = distf.pdf_to_mean(Px)
sXkm = distf.pdf_to_var(Px)
mYkm = distf.pdf_to_mean(Py)
sYkm = distf.pdf_to_var(Py)

mXc = distf.pdf_to_mean(PXc)
sXc = distf.pdf_to_var(PXc)
mYc = distf.pdf_to_mean(PYc)
sYc = distf.pdf_to_var(PYc)
corrXYc = distf.pdf_to_cov(PXYc)/np.sqrt(sXc*sYc)

ptrunc     = prob_trunc(PXY,pT_an,pE_an,5)
ptrunc_p   = prob_trunc(priorXY,priorT,priorE,epsm)
ptrunc_pol = prob_trunc(postlXY,pTl,pEl,epsm)
ptrunc_poh = prob_trunc(posthXY,pTh,pEh,epsm)

## Print results
print('Data statistics:')
print('Mean(X) \t Mean(Y) \t Var(X) \t Var(Y) \t Corr(X,Y) \t P(A)')
print('%.4f %.4f %.4f %.4f %.4f' % (xd.mean(),ye.mean(),xd.var(),ye.var(),np.corrcoef(xu,yu)[0,1]))
print('%.4f %.4f %.4f %.4f %.4f %.4f' % (mX,mY,sX,sY,corrXY,ptrunc))
print('%.4f %.4f %.4f %.4f %.4f %.4f' % (mXp,mYp,sXp,sYp,corrXYp,ptrunc_p))
print('%.4f %.4f %.4f %.4f %.4f %.4f' % (mXpol,mYpol,sXpol,sYpol,corrXYpol,ptrunc_pol))
print('%.4f %.4f %.4f %.4f %.4f %.4f' % (mXpoh,mYpoh,sXpoh,sYpoh,corrXYpoh,ptrunc_poh))
print('%.4f %.4f %.4f %.4f %.4f' % (mXc,mYc,sXc,sYc,corrXYc))
print('%.4f %.4f %.4f %.4f' % (mXkm,mYkm,sXkm,sYkm))

print('True statistics:')
print('Mean(X) \t Mean(Y) \t Var(X) \t Var(Y) \t Corr(X,Y)')
print('%.4f \t %.4f \t %.4f \t %.4f \t %.4f\n' % (mX,mY,sX,sY,corrXY))

print('Prior statistics:')
print('Mean(X) \t Mean(Y) \t Var(X) \t Var(Y) \t Corr(X,Y)')
print('%.4f \t %.4f \t %.4f \t %.4f \t %.4f\n' % (mXp,mYp,sXp,sYp,corrXYp))

print('Posterior statistics (low):')
print('Mean(X) \t Mean(Y) \t Var(X) \t Var(Y) \t Corr(X,Y)')
print('%.4f \t %.4f \t %.4f \t %.4f \t %.4f\n' % (mXpol, mYpol, sXpol, sYpol, corrXYpol))

print('Posterior statistics (high):')
print('Mean(X) \t Mean(Y) \t Var(X) \t Var(Y) \t Corr(X,Y)')
print('%.4f \t %.4f \t %.4f \t %.4f \t %.4f\n' % (mXpoh, mYpoh, sXpoh, sYpoh, corrXYpoh))

print('KM statistics:')
print('Mean(X) \t Mean(Y) \t Var(X) \t Var(Y)')
print('%.4f \t %.4f \t %.4f \t %.4f\n' % (mXkm,mYkm,sXkm,sYkm))

print('Copula statistics:')
print('Mean(X) \t Mean(Y) \t Var(X) \t Var(Y) \t Corr(X,Y)')
print('%.4f \t %.4f \t %.4f \t %.4f \t %.4f\n' % (mXc,mYc,sXc,sYc,corrXYc))

## Plot the results
xmin = minX-5
ymin = minY-5

if filename=='UOFM':
    xmax = 120
    ymax = 120
else:
    xmax = 90
    ymax = 90

hx = plt.figure()
if nr!='':
    ERl = r'(ER$_'+nr+r'^l$)'
    ERh = r'(ER$_'+nr+r'^h$)'
else:
    ERl = r'(ER$^l$)'
    ERh = r'(ER$^h$)'
if filename=='UOFM':
    legendlist   = ('Kaplan-Meier', 'Copula', 'Prior', 'Posterior '+ERl, 'Posterior '+ERh)
    legendlist_t = legendlist[2:]
else:
    legendlist   = ('Kaplan-Meier', 'Reference', 'Prior', 'Posterior '+ERl, 'Posterior '+ERh)
    legendlist_t = legendlist[1:]


# Marginal X
if filename=='UOFM':
    plt.plot(vx,1.-Sx,'r^',v,1.-SXc,'bo',v,1.-SXp,'g*',v,1.-SXpol,'kx',v,1.-SXpoh,'c+')
else:
    plt.plot(vx,1.-Sx,'r^',v,1.-SX,'bo',v,1.-SXp,'g*',v,1.-SXpol,'kx',v,1.-SXpoh,'c+')
plt.xlim(xmin,xmax)
plt.xlabel('X',fontsize=16)
plt.ylabel('Cumulative probability',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(legendlist,loc='lower right',fontsize=12)
plt.savefig(output_dir+'marginalX_'+filename_res+'.eps',format='eps')
plt.show()

# Marginal Y
hy = plt.figure()
if filename=='UOFM':
    plt.plot(vy,1.-Sy,'r^',v,1.-SYc,'bo',v,1.-SYp,'g*',v,1.-SYpol,'kx',v,1.-SYpoh,'c+')
else:
    plt.plot(vy,1.-Sy,'r^',v,1.-SY,'bo',v,1.-SYp,'g*',v,1.-SYpol,'kx',v,1.-SYpoh,'c+')
plt.xlim(xmin,xmax)
plt.xlabel('Y',fontsize=16)
plt.ylabel('Cumulative probability',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(legendlist,loc='lower right',fontsize=12)
plt.savefig(output_dir+'marginalY_'+filename_res+'.eps',format='eps')
plt.show()

# Truncation X
htx = plt.figure()
if filename=='UOFM':
    plt.plot(v,1.-STp,'g*',v,np.cumsum(pTl),'kx',v,np.cumsum(pTh),'c+')
else:
    plt.plot(v,1.-ST_an,'bo',v,1.-STp,'g*',v,np.cumsum(pTl),'kx',v,np.cumsum(pTh),'c+')
plt.xlim(xmin-10,xmax)
plt.xlabel(r'$T^X$',fontsize=16)
plt.ylabel('Cumulative probability',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(legendlist_t,loc='lower right',fontsize=12)
plt.savefig(output_dir+'marginalTx_'+filename_res+'.eps',format='eps')
plt.show()

# Epsilon
he = plt.figure()
if filename=='UOFM':
    plt.plot(v-epsm,1.-SEp,'g*',v-epsm,np.cumsum(pEl),'kx',v-epsm,np.cumsum(pEh),'c+')
else:
    plt.plot(v-5,1.-SE_an,'bo',v-epsm,1.-SEp,'g*',v-epsm,np.cumsum(pEl),'kx',v-epsm,np.cumsum(pEh),'c+')
plt.xlim(-epsm,maxE)
plt.xlabel(r'$\epsilon$',fontsize=16)
plt.ylabel('Cumulative probability',fontsize=16)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(legendlist_t,loc='lower right',fontsize=12)
plt.savefig(output_dir+'marginalEps_'+filename_res+'.eps',format='eps')
plt.show()

# Joint distribution
levels = np.linspace(1e-4,5e-3,8)
X,Y    = np.meshgrid(v,v)

nrows = 2
ncols = 2
fig, axes = plt.subplots(nrows,ncols,sharex='all',sharey='all')
if filename=='UOFM':
    CS = axes[0][0].contour(X,Y,PXYc.T,cmap='hot',linewidths=1)
    axes[0][0].set_title('Frank copula',fontsize=10)
else:
    CS = axes[0][0].contour(X,Y,PXY.T,cmap='hot',linewidths=1)
    axes[0][0].set_title('Reference distribution',fontsize=10)
axes[0][1].contour(X,Y,priorXY.T,CS.levels,cmap='hot',linewidths=1)
axes[0][1].set_title('Prior distribution',fontsize=10)
axes[1][0].contour(X,Y,postlXY.T,CS.levels,cmap='hot',linewidths=1)
axes[1][0].set_title(r'Posterior distribution '+ERl,fontsize=10)
axes[1][1].contour(X,Y,posthXY.T,CS.levels,cmap='hot',linewidths=1)
axes[1][1].set_title('Posterior distribution '+ERh,fontsize=10)
for i in range(nrows):
    axes[i][0].set_ylabel('Y',fontsize=12)
for i in range(ncols):
    axes[-1][i].set_xlabel('X',fontsize=12)
for i in range(nrows):
    for j in range(ncols):
        axes[i][j].set_xlim(xmin,xmax)
        axes[i][j].set_ylim(ymin,ymax)
        axes[i][j].tick_params(axis='both', which='major', labelsize=10)
        axes[i][j].tick_params(axis='both', which='minor', labelsize=10)
plt.savefig(output_dir+'dist2D_'+filename_res+'.eps',format='eps')
plt.show()

## Permutation tests

if filename=='poisson':
    N      = 500
    Xs_or  = distf.sample_from_cdf(1.-SX,N)
    Ys_or  = distf.sample_from_cdf(1.-SY,N)
    Xs_po  = distf.sample_from_cdf(1.-SXpol,N)
    Ys_po  = distf.sample_from_cdf(1.-SYpol,N)
    f_mean = lambda x,y,axis: np.mean(x,axis=axis)-np.mean(y,axis=axis)
    f_var  = lambda x,y,axis: np.var(x,axis=axis)-np.var(y,axis=axis)


    stat_valX = 0
    stat_valY = 0
    pvalX     = 0
    pvalY     = 0
    stat_varX = 0
    stat_varY = 0
    pvalVX    = 0
    pvalVY    = 0
    nperm     = 1
    nsamples  = 100000
    for i in range(nperm):
        res        = permutation_test((Xs_or[:,np.newaxis],Xs_po[:,np.newaxis]),f_mean,
                                     n_resamples=nsamples,vectorized=True)
        stat_valX += res.statistic
        pvalX     += res.pvalue

        res        = permutation_test((Ys_or[:,np.newaxis],Ys_po[:,np.newaxis]),f_mean,
                                     n_resamples=nsamples,vectorized=True)
        stat_valY += res.statistic
        pvalY     += res.pvalue

        res        = permutation_test((Xs_or[:,np.newaxis],Xs_po[:,np.newaxis]),f_var,
                                     n_resamples=nsamples,vectorized=True)
        stat_varX += res.statistic
        pvalVX    += res.pvalue

        res        = permutation_test((Ys_or[:,np.newaxis],Ys_po[:,np.newaxis]),f_var,
                                     n_resamples=nsamples,vectorized=True)
        stat_varY += res.statistic
        pvalVY    += res.pvalue

    stat_valX /= nperm
    pvalX     /= nperm
    stat_valY /= nperm
    pvalY     /= nperm
    stat_varX /= nperm
    pvalVX    /= nperm
    stat_varY /= nperm
    pvalVY    /= nperm
    print(stat_valX, pvalX)
    print(stat_valY, pvalY)
    print(stat_varX, pvalVX)
    print(stat_varY, pvalVY)

## QQ-plots

    Xss_or = np.sort(Xs_or)
    Yss_or = np.sort(Ys_or)
    Xss_po = np.sort(Xs_po)
    Yss_po = np.sort(Ys_po)

    hqqx = plt.figure()
    plt.plot(Xss_or,Xss_or,'r--')
    plt.plot(Xss_po,Xss_or,'bx')
    plt.plot(mX*np.ones((100,)),np.linspace(Xss_or.min(),Xss_or.max(),100),'k--')
    plt.text(mX*1.01,Xss_or.max()*0.99,r'$\mu_X$',fontsize=12)
    plt.plot((mX+np.sqrt(sX))*np.ones((100,)),np.linspace(Xss_or.min(),Xss_or.max(),100),'k--')
    plt.text((mX+np.sqrt(sX))*1.01,Xss_or.max()*0.99,r'$\mu_X+\sigma_X$',fontsize=12)
    plt.plot((mX-np.sqrt(sX))*np.ones((100,)),np.linspace(Xss_or.min(),Xss_or.max(),100),'k--')
    plt.text((mX-np.sqrt(sX))*0.88,Xss_or.max()*0.99,r'$\mu_X-\sigma_X$',fontsize=12)
    plt.xlabel('Quantiles of the posterior distribution',fontsize=16)
    plt.ylabel('Quantiles of the reference distribution',fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.savefig(output_dir+'QQplotX_poisson.eps',format='eps')
    plt.show()

    hqqy = plt.figure()
    plt.plot(Yss_or,Yss_or,'r--')
    plt.plot(Yss_po,Yss_or,'bx')
    plt.plot(mY*np.ones((100,)),np.linspace(Yss_or.min(),Yss_or.max(),100),'k--')
    plt.text(mY*1.01,Yss_or.max()*0.99,r'$\mu_Y$',fontsize=12)
    plt.plot((mY+np.sqrt(sY))*np.ones((100,)),np.linspace(Yss_or.min(),Yss_or.max(),100),'k--')
    plt.text((mY+np.sqrt(sY))*1.01,Yss_or.max()*0.99,r'$\mu_Y+\sigma_Y$',fontsize=12)
    plt.plot((mY-np.sqrt(sY))*np.ones((100,)),np.linspace(Yss_or.min(),Yss_or.max(),100),'k--')
    plt.text((mY-np.sqrt(sY))*0.88,Yss_or.max()*0.99,r'$\mu_Y-\sigma_Y$',fontsize=12)
    plt.xlabel('Quantiles of the posterior distribution',fontsize=16)
    plt.ylabel('Quantiles of the reference distribution',fontsize=16)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.savefig(output_dir+'QQplotY_poisson.eps',format='eps')
    plt.show()

##


##

