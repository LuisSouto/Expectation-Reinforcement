/* 
	 Implementation of the Expectation Maximization to calibrate the B-RUP under the pressence of
   bivariate censoring and truncation.
*/
 
#include "brup_em.h"

template<typename RealType>
void BRUP_ApplyPrior(const int &m, RealType* pA, RealType* pB, RealType *pC, RealType* pTx, RealType* pE,
										 const RealType* pAe, const RealType* pBe, const RealType* pCe, const RealType* pTxe,
										 const RealType* pEe, const RealType* pAp, const RealType* pBp, const RealType* pCp,
										 const RealType* pTxp, const RealType* pEp, const RealType &cA, const RealType &cB,
										 const RealType &cC, const RealType &cT, const RealType &cE, const RealType &r){
										 
	unique_ptr<RealType[]> BA(new RealType[m]),  BB(new RealType[m]),  BC(new RealType[m]),  BTx(new RealType[m]),  BE(new RealType[m]);
	unique_ptr<RealType[]> NA(new RealType[m]),  NB(new RealType[m]),  NC(new RealType[m]),  NTx(new RealType[m]),  NE(new RealType[m]);
	unique_ptr<RealType[]> BAo(new RealType[m]), BBo(new RealType[m]), BCo(new RealType[m]), BTxo(new RealType[m]), BEo(new RealType[m]);
	unique_ptr<RealType[]> NAo(new RealType[m]), NBo(new RealType[m]), NCo(new RealType[m]), NTxo(new RealType[m]), NEo(new RealType[m]);			
			
	/* Compute the urn compositions for the priors */
	urn_init(m,BAo.get(),NAo.get(),pAp,cA);
	urn_init(m,BBo.get(),NBo.get(),pBp,cB);
	urn_init(m,BCo.get(),NCo.get(),pCp,cC);
	urn_init(m,BTxo.get(),NTxo.get(),pTxp,cT);
	urn_init(m,BEo.get(),NEo.get(),pEp,cE);	
	
	/* Compute the urn compositions for the final estimates */
	urn_init(m,BA.get(),NA.get(),pAe,r);
	urn_init(m,BB.get(),NB.get(),pBe,r);
	urn_init(m,BC.get(),NC.get(),pCe,r);
	urn_init(m,BTx.get(),NTx.get(),pTxe,r);
	urn_init(m,BE.get(),NE.get(),pEe,r);

	add_u_to_v(m,BA.get(),BAo.get());
	add_u_to_v(m,NA.get(),NAo.get());
	add_u_to_v(m,BB.get(),BBo.get());
	add_u_to_v(m,NB.get(),NBo.get());
	add_u_to_v(m,BC.get(),BCo.get());
	add_u_to_v(m,NC.get(),NCo.get());
	add_u_to_v(m,BTx.get(),BTxo.get());
	add_u_to_v(m,NTx.get(),NTxo.get());
	add_u_to_v(m,BE.get(),BEo.get());
	add_u_to_v(m,NE.get(),NEo.get());
	
	urn_pdf(m,pA,BA.get(),NA.get());
	urn_pdf(m,pB,BB.get(),NB.get());
	urn_pdf(m,pC,BC.get(),NC.get()); 
	urn_pdf(m,pTx,BTx.get(),NTx.get());
	urn_pdf(m,pE,BE.get(),NE.get());
}

template <typename RealType>
void BRUP_EMopt(const int &m, RealType *pA, RealType *pB, RealType *pC, RealType *pTx, RealType *pE,
								const RealType *pA0, const RealType *pB0, const RealType *pC0, const RealType *pTx0, 
								const RealType *pE0, const int &n, const int *X, const int *Y, const int *TX,
								const int *TY, const int *dx, const int *dy, const int *eps, const RealType &cA,
								const RealType &cB, const RealType &cC, const RealType &cT, const RealType &cE, const RealType &r,
								const int &N, const bool &IF_TRUNC, const RealType &tol){

	/*-------------------------------------------------------------------*/
	/*----------------------DATA INITIALIZATION--------------------------*/
	/*-------------------------------------------------------------------*/

	/* Select upper bounds for each variable */
	int maxX   = *max_element(X,X+n),     maxY  = *max_element(Y,Y+n),     maxXY = max(maxX,maxY);
	int minE   = *min_element(eps,eps+n), maxE  = *max_element(eps,eps+n);	
	int ma     = maxXY+2,                 mb    = maxX+2;
	int mc     = maxY+2,                  mtx   = *max_element(TX,TX+n)+1, meps  = maxE-minE+1;
	int mty    = mtx+meps-1;
	
	/* Find repeated values */
	int *Xu = new int[n], *Yu = new int[n], *dxu = new int[n], *dyu = new int[n];
	int  nu = 0;
	find_unique_LTRC(n,nu,Xu,Yu,dxu,dyu,X,Y,dx,dy);
	int *uni_counts = new int[nu];
	count_unique_LTRC(n,nu,uni_counts,Xu,Yu,dxu,dyu,X,Y,dx,dy);
	
	int	*minxy = new int[nu];
	for(int i=0;i<nu;i++){
		minxy[i] = min(Xu[i],Yu[i]);
	}
	
	/* Types of censoring */
	int *censtype = new int[nu];
	for(int i=0;i<nu;i++){
		if((dxu[i]==1)&&(dyu[i]==1)){
			censtype[i] = 0;
		}
		else if((dxu[i]==0)&&(dyu[i]==1)){
			censtype[i] = 1;
		}
		else if((dxu[i]==1)&&(dyu[i]==0)){
			censtype[i] = 2;
		}
		else if((dxu[i]==0)&&(dyu[i]==0)){
			censtype[i] = 3;
		}		
	}
	
	/* Shift epsilon to make it positive */
	int *epsp = new int[n];
	int  eps0 = -minE;
	for(int i=0;i<n;i++){
		epsp[i] = eps[i]+eps0;
	}	

	RealType *BA   = new RealType[ma],    *BB   = new RealType[mb],   *BC   = new RealType[mc],   *BTx   = new RealType[mtx], *BE   = new RealType[meps];
	RealType *NA   = new RealType[ma],    *NB   = new RealType[mb],   *NC   = new RealType[mc],   *NTx   = new RealType[mtx], *NE   = new RealType[meps];
	RealType *SA   = new RealType[m],     *SB   = new RealType[m],    *SC   = new RealType[m],    *STx   = new RealType[m],   *SE   = new RealType[m]; 
	RealType *pAt  = new RealType[ma],    *pBt  = new RealType[mb],   *pCt  = new RealType[mc],   *pTxt  = new RealType[mtx], *pEt  = new RealType[meps];
	RealType *pA0_ = new RealType[ma],    *pB0_ = new RealType[mb],   *pC0_ = new RealType[mc],   *pTx0_ = new RealType[mtx], *pE0_ = new RealType[meps];	
	RealType *SAt  = new RealType[ma],    *SBt  = new RealType[mb],   *SCt  = new RealType[mc],   *STxt  = new RealType[mtx], *SEt  = new RealType[meps];
	RealType *pABC = new RealType[m*m*m], *P    = new RealType[m*m],  *SXY  = new RealType[m*m],  *PTE   = new RealType[mty*mtx];	
	RealType *PX   = new RealType[m],     *SX   = new RealType[m],    *PY   = new RealType[m],    *SY    = new RealType[m];
	RealType *pASB = new RealType[ma*m],  *pASC = new RealType[ma*m], *pApB = new RealType[ma*m], *pApC  = new RealType[ma*m];		
	RealType *pAc  = new RealType[ma];
	
	/* Normalize initial distributions */
	copy(pA0,pA0+ma,pA0_);
	copy(pB0,pB0+mb,pB0_);
	copy(pC0,pC0+mc,pC0_);
	copy(pTx0,pTx0+mtx,pTx0_);
	copy(pE0,pE0+meps,pE0_);
	div(ma,pA0_,sum(ma,pA0_));
	div(mb,pB0_,sum(mb,pB0_));
	div(mc,pC0_,sum(mc,pC0_));
	div(mtx,pTx0_,sum(mtx,pTx0_));
	div(meps,pE0_,sum(meps,pE0_));					

	/* Compute the urn compositions for the initial estimates */
	urn_init(ma,BA,NA,pA0_,(RealType)1.);
	urn_init(mb,BB,NB,pB0_,(RealType)1.);
	urn_init(mc,BC,NC,pC0_,(RealType)1.);
	urn_init(mtx,BTx,NTx,pTx0_,(RealType)1.);
	urn_init(meps,BE,NE,pE0_,(RealType)1.);
	
	/* Hazard functions of truncation variables */
	int *mtxd = new int[m], *stxd = new int[m], *mtyd = new int[m], *styd = new int[m]; 
	mt(mtxd,m-1,n,TX);
	st(stxd,m-1,n,TX);
	mt(mtyd,m-1,n,epsp);
	st(styd,m-1,n,epsp);
	
	/* Initialize vectors with zeros */
	zeros(m*m*m,pABC);
	zeros(m*m,P);
	zeros(mtx*mty,PTE);
	zeros(m,SA);
	zeros(m,SB);
	zeros(m,SC);
	zeros(m,STx);
	zeros(m,SE);
	zeros(m,pA);
	zeros(m,pB);
	zeros(m,pC);
	zeros(m,pTx);
	zeros(m,pE);
	zeros(ma*m,pASB);
	zeros(ma*m,pASC);
	
	RealType loglikelihood, loglikelihood_old = 1.e5;
	
	/*-------------------------------------------------------------------------------------*/
	/*----------------------------EXPECTATION-MAXIMIZATION ROUTINE-------------------------*/
	/*-------------------------------------------------------------------------------------*/
	
	for(int k=0;k<N;k++){
		/* Precompute all distributions from urn compositions */
		urn_surv(ma,SA,BA,NA);
		urn_surv(mb,SB,BB,NB);
		urn_surv(mc,SC,BC,NC); 
		urn_surv(mtx,STx,BTx,NTx);
		urn_surv(meps,SE,BE,NE);
		pdf_surv(ma,pA,SA);
		pdf_surv(mb,pB,SB);
		pdf_surv(mc,pC,SC);
		pdf_surv(mtx,pTx,STx);
		pdf_surv(meps,pE,SE);
		
		Compute_pXY(pABC,P,m,ma,mb,mc,pA,pB,pC);
		Compute_pTrunc(PTE,mtx,mty,meps,eps0,pTx,STx,SE);
		Convolution(PX,m,ma,mb,pA,pB);
		Convolution(PY,m,ma,mc,pA,pC);
		surv_pdf(m,SX,PX);
		surv_pdf(m,SY,PY);
		surv_pdf2(m,SXY,P,PX,PY);
		
		/* Precompute censoring terms */
		DoubleCens_Conv(pASB,pA,SB,m,ma,mb);
		DoubleCens_Conv(pASC,pA,SC,m,ma,mc);
		SingleCens_Conv(pApB,pA,pB,m,ma,mb);
		SingleCens_Conv(pApC,pA,pC,m,ma,mc);
		
		/* Restart estimates */
		zeros(ma,BA);
		zeros(mb,BB);
		zeros(mc,BC);
		zeros(ma,NA);
		zeros(mb,NB);
		zeros(mc,NC);
		copy(mtxd,mtxd+mtx,BTx);
		copy(stxd,stxd+mtx,NTx);
		copy(mtyd,mtyd+meps,BE);
		copy(styd,styd+meps,NE);
			
		/* Truncation component */
		RealType pTX = 0., M = 0., errT = 0., errE = 0.;
		if(IF_TRUNC){
			pTxt[0] = (eps0>=0)?1.-SE[eps0]:0.;
			for(int e=eps0+1;e<min(m+eps0,meps);e++){
				pTxt[0] += pE[e]*SY[e-eps0-1];
			}
			pTxt[0] = pTx[0]*(1-pTxt[0]);
			pTX += pTxt[0];
			for(int t=1;t<mtx;t++){
				if(t<=eps0){
					pTxt[t] = SX[t-1]*(1.-SE[eps0-t]);
				}
				else{
					pTxt[t] = SXY[t-eps0-1+(t-1)*m]*pE[0];
				}
				for(int e=max(eps0-t,0)+1;e<min(m+eps0-t,meps);e++){
					pTxt[t] += pE[e]*SXY[t+e-eps0-1+(t-1)*m];						
				}
				pTxt[t] = pTx[t]*(1-pTxt[t]);
				pTX    += pTxt[t];
				if(t>m-meps+eps0){
					errT += pTx[t]*SXY[m-1+(t-1)*m]*SE[m-1-t+eps0];
				}
			}
			M = n/(1.-pTX);
			
			for(int e=0;e<meps;e++){
				pEt[e] = pTx[0]*((e>eps0)?SY[e-eps0-1]:1.);
				for(int t=1;t<eps0+1-e;t++){
					pEt[e] += pTx[t]*SX[t-1];
				}
				for(int t=max(eps0-e,0)+1;t<min(m-e+eps0,mtx);t++){
					pEt[e] += pTx[t]*SXY[t+e-eps0-1+(t-1)*m];
				}				
				pEt[e] = pE[e]*(1-pEt[e]);
				for(int t=m+eps0-e;t<mtx;t++){
					errE += pE[e]*pTx[t]*SXY[m-1+(t-1)*m];
				}
			}
			
			zeros(ma,pAt);
			zeros(mc,pCt);
			for(int b=0;b<mb;b++){		
				pBt[b] = 0.;
				for(int c=0;c<mc;c++){
					RealType sumA = 0.;
					int idbc = max(mtx-b,mty-c);
					for(int a=0;a<min(ma,idbc);a++){
						sumA    = pA[a]*pB[b]*pC[c]*PTE[min(a+c,mty-1)+(min(a+b,mtx-1))*mty];
						pAt[a] += sumA;
						pBt[b] += sumA;
						pCt[c] += sumA;
					}
					for(int a=idbc;a<ma;a++){
						pAt[a] += pA[a]*pB[b]*pC[c];
					}
					if(idbc<ma){				
						sumA    = pC[c]*pB[b]*((idbc>0)?SA[idbc-1]:1.);
						pCt[c] += sumA;
						pBt[b] += sumA;
					}
				}
			}
			for(int a=0;a<ma;a++){
				pAt[a] = pA[a]-pAt[a];
			}			
			for(int b=0;b<mb;b++){
				pBt[b] = pB[b]-pBt[b];
			}			
			for(int c=0;c<mc;c++){
				pCt[c] = pC[c]-pCt[c];
			}

			surv_pdf(ma,SAt,pAt,pTX);
			surv_pdf(mb,SBt,pBt,pTX);
			surv_pdf(mc,SCt,pCt,pTX);
			surv_pdf(mtx,STxt,pTxt,pTX);
			surv_pdf(meps,SEt,pEt,pTX);
		}

	  /* Censoring component */
		for(int i=0;i<nu;i++){
			int idXY = Yu[i]+Xu[i]*m;
			switch(censtype[i]){
				case 0:{
					if(P[idXY]>0){
						for(int j=0;j<min(minxy[i]+1,ma);j++){
							RealType pabc = uni_counts[i]*pABC[j+Yu[i]*m+Xu[i]*m*m]/P[idXY];
							BA[j]        += pabc;
							BB[Xu[i]-j]  += pabc;
							BC[Yu[i]-j]  += pabc;
						}
					}
					break;
				}
				case 1:{
					RealType px = SXY[idXY-1]-SXY[idXY];
					SingleCensoring(pAc,BA,BB,BC,pA,pB,pC,SB,pApC,Xu[i],Yu[i],px,ma,mb,mc,uni_counts[i]);					
					break;
				}
				case 2:{
					RealType py = SXY[idXY-m]-SXY[idXY];
					SingleCensoring(pAc,BA,BC,BB,pA,pC,pB,SC,pApB,Yu[i],Xu[i],py,ma,mb,mc,uni_counts[i]);
					break;
				}
				case 3:{
					RealType p = SXY[idXY];
					DoubleCensoring(pAc,BA,BB,BC,pA,pB,pC,SB,SC,pASB,pASC,Xu[i],Yu[i],p,ma,mb,mc,uni_counts[i]);
					break;
				}
			}
		}
		
		/* Update the urns with new estimates */
		NA[ma-1] = BA[ma-1];
		for(int i=ma-2;i>=0;i--){
			NA[i] = NA[i+1]+BA[i];
		}

		NB[mb-1] = BB[mb-1];		
		for(int i=mb-2;i>=0;i--){
			NB[i] = NB[i+1]+BB[i];
		}		

		NC[mc-1] = BC[mc-1];
		for(int i=mc-2;i>=0;i--){
			NC[i] = NC[i+1]+BC[i];
		}		
				
		if(IF_TRUNC){
			add_u_to_v(ma,BA,pAt,M);
			add_u_to_v(ma,NA,SAt,M);
			add_u_to_v(ma,NA,pAt,M);		
			add_u_to_v(mb,BB,pBt,M);
			add_u_to_v(mb,NB,SBt,M);
			add_u_to_v(mb,NB,pBt,M);
			add_u_to_v(mc,BC,pCt,M);
			add_u_to_v(mc,NC,SCt,M);
			add_u_to_v(mc,NC,pCt,M);
			add_u_to_v(mtx,BTx,pTxt,M);
			add_u_to_v(mtx,NTx,STxt,M);
			add_u_to_v(mtx,NTx,pTxt,M);		
			add_u_to_v(meps,BE,pEt,M);
			add_u_to_v(meps,NE,SEt,M);
			add_u_to_v(meps,NE,pEt,M);
			loglikelihood = loglike_ltrc_2D(n,X,Y,TX,TY,dx,dy,m,SXY);
		}
		else{
			loglikelihood = loglike_rc_2D(n,X,Y,dx,dy,m,SXY);
		}
		cout<<k<<" "<<sum(ma,pA)<<" "<<sum(mb,pB)<<" "<<sum(mc,pC)<<" "<<sum(mtx,pTx)<<" "<<sum(meps,pE)<<" "<<sum(mtx*mty,PTE)<<endl;
		cout<<k<<" "<<sum(ma,SA)<<" "<<sum(mb,SB)<<" "<<sum(mc,SC)<<" "<<sum(mtx,STx)<<" "<<sum(meps,SE)<<" "<<PTE[mtx*mty-1]<<endl;
		cout<<k<<" "<<loglikelihood<<" "<<abs((loglikelihood-loglikelihood_old)/loglikelihood_old)<<" "<<M<<endl;
				
		if(abs((loglikelihood-loglikelihood_old)/loglikelihood_old)<=tol){
			break;
		}
		loglikelihood_old = loglikelihood;
	}
}

/* Template instantiations */

template void BRUP_EMopt(const int &m, float *pA, float *pB, float *pC, float *pTx, float *pE,
												 const float *pA0, const float *pB0, const float *pC0, const float *pTx0, 
												 const float *pE0, const int &n, const int *X, const int *Y, const int *TX,
												 const int *TY, const int *dx, const int *dy, const int *eps, const float &cA,
												 const float &cB, const float &cC, const float &cT, const float &cE,
												 const float &r, const int &N, const bool &IF_TRUNC, const float &tol);
												 
template void BRUP_ApplyPrior(const int &m, float* pA, float* pB, float *pC, float* pTx, float* pE,
										          const float* pAe, const float* pBe, const float* pCe, const float* pTxe,
				 										  const float* pEe, const float* pAp, const float* pBp, const float* pCp,
										 					const float* pTxp, const float* pEp, const float &cA, const float &cB,
				  										const float &cC, const float &cT, const float &cE, const float &r);												 
												 
template void BRUP_EMopt(const int &m, double *pA, double *pB, double *pC, double *pTx, double *pE,
												 const double *pA0, const double *pB0, const double *pC0, const double *pTx0, 
												 const double *pE0, const int &n, const int *X, const int *Y, const int *TX,
												 const int *TY, const int *dx, const int *dy, const int *eps, const double &cA,
												 const double &cB, const double &cC, const double &cT, const double &cE,												 
												 const double &r, const int &N, const bool &IF_TRUNC, const double &tol);	
												 
template void BRUP_ApplyPrior(const int &m, double* pA, double* pB, double *pC, double* pTx, double* pE,
										          const double* pAe, const double* pBe, const double* pCe, const double* pTxe,
				 										  const double* pEe, const double* pAp, const double* pBp, const double* pCp,
										 					const double* pTxp, const double* pEp, const double &cA, const double &cB,
										 					const double &cC, const double &cT, const double &cE, const double &r);														 											 
