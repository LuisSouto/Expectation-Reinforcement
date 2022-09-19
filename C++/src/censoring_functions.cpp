/* Source code for censoring routines related to the B-RUP calibration via EM */

#include "censoring_functions.h"


void find_unique_LTRC(const int &n, int &nu, int* Xu, int *Yu, int *dxu, int *dyu,
								      const int* X, const int *Y, const int *dx, const int *dy){
	nu     = 1;
	Xu[0]  = X[0];
	Yu[0]  = Y[0];
	dxu[0] = dx[0];
	dyu[0] = dy[0];
	for(int i=1;i<n;i++){
		int uni_id = 0;
		for(int j=0;j<nu;j++){
			if((X[i]!=Xu[j])||(Y[i]!=Yu[j])||(dx[i]!=dxu[j])||(dy[i]!=dyu[j])){
				uni_id++;
			}
		}
		if(uni_id==nu){
			Xu[nu]  = X[i];
			Yu[nu]  = Y[i];
			dxu[nu] = dx[i];
			dyu[nu] = dy[i];
			nu++;
		}
	}
}

void count_unique_LTRC(const int &n, const int &nu, int* uni_counts, const int* Xu,
										   const int* Yu, const int* dxu, const int* dyu, const int* X,
										   const int* Y, const int* dx, const int* dy){
	for(int i=0;i<n;i++){
		for(int j=0;j<nu;j++){
			if((X[i]==Xu[j])&&(Y[i]==Yu[j])&&(dx[i]==dxu[j])&&(dy[i]==dyu[j])){
				uni_counts[j]++;
				break;
			}
		}
	}										   
}

template <typename RealType>
void kernel_Compute_pXY(RealType *pABC, RealType *P, const RealType *pA, const RealType *pB, const RealType *pC,
												const int i, const int j, const int m, const int mxy2, const int lmin, const int lmax){
	for(int l=lmin;l<=lmax;l++){
		pABC[l+j*m+i*mxy2] = pA[l]*pB[i-l]*pC[j-l];
		P[j+i*m] += pABC[l+j*m+i*mxy2];            // P(X = i,Y = j)
	}
}

template <typename RealType>
void Compute_pXY(RealType *pABC,RealType *P,const int m,const int ma,const int mb,const int mc,const RealType *pA,
								 const RealType *pB,const RealType *pC, const int minX, const int minY){
	int mxy2 = m*m;								 
	int mab  = min(ma,mb), mabxy = min(ma+mb-1,m), macxy = min(ma+mc-1,m);
	
	for(int i=minX;i<mab;i++){
		for(int j=minY;j<=i;j++){
			P[j+i*m] = 0.;
			kernel_Compute_pXY(pABC,P,pA,pB,pC,i,j,m,mxy2,max(0,j-mc+1),j);
		}
		for(int j=max(minY,i+1);j<min(i+mc,m);j++){
			P[j+i*m] = 0.;
			kernel_Compute_pXY(pABC,P,pA,pB,pC,i,j,m,mxy2,max(0,j-mc+1),i);			
		}
	}
	
	for(int i=mab;i<ma;i++){
		for(int j=max(minY,i-mb+1);j<=i;j++){
			P[j+i*m] = 0.;
			kernel_Compute_pXY(pABC,P,pA,pB,pC,i,j,m,mxy2,max(i-mb+1,j-mc+1),j);
		}
		for(int j=i+1;j<min(i+mc,m);j++){
			P[j+i*m] = 0.;
			kernel_Compute_pXY(pABC,P,pA,pB,pC,i,j,m,mxy2,max(i-mb+1,j-mc+1),i);			
		}		
	}
	
	for(int i=ma;i<mabxy;i++){
		for(int j=max(minY,i-mb+1);j<ma;j++){
			P[j+i*m] = 0.;
			kernel_Compute_pXY(pABC,P,pA,pB,pC,i,j,m,mxy2,max(max(0,i-mb+1),j-mc+1),j);			
		}
		for(int j=ma;j<macxy;j++){
			P[j+i*m] = 0.;
			kernel_Compute_pXY(pABC,P,pA,pB,pC,i,j,m,mxy2,max(max(0,i-mb+1),j-mc+1),ma-1);
		}		
	}	
}

template <typename RealType>
void Compute_pTrunc(RealType *P,const int mtx,const int mty,const int meps,const int eps0,const RealType *pT,
									  const RealType *ST, const RealType *Se, const RealType cumTx){
	for(int y=0;y<meps-eps0;y++){
		P[y] = pT[0]*(1.-Se[y+eps0]);
		for(int x=1;x<=min(y+eps0,mtx-1);x++){
			P[y+x*mty] = P[y+(x-1)*mty]+pT[x]*(1.-Se[y+eps0-x]); // P[Tx<=x, Ty<=y]
		}
		for(int x=y+eps0+1;x<mtx;x++){
			P[y+x*mty] = P[y+(x-1)*mty];
		}
	}
	for(int y=meps-eps0;y<mty;y++){
		P[y] = pT[0];
		for(int x=1;x<=min(y+eps0-meps,mtx-1);x++){
			P[y+x*mty] = P[y+(x-1)*mty]+pT[x];
		}
		for(int x=y+eps0-meps+1;x<=min(y+eps0,mtx-1);x++){
			P[y+x*mty] = P[y+(x-1)*mty]+pT[x]*(1.-Se[y+eps0-x]);
		}
		for(int x=y+eps0+1;x<mtx;x++){
			P[y+x*mty] = P[y+(x-1)*mty];
		}
	}
}

template <typename RealType>
void Convolution(RealType *P,const int m,const int ma,const int mb,const RealType *pA,const RealType *pB){
	for(int i=0;i<m;i++){
		P[i] = 0.;
		for(int l=max(0,i-mb+1);l<min(i+1,ma);l++){
			P[i] += pA[l]*pB[i-l]; // P(X = i)
		}
	}
}

template <typename RealType>
void SingleCens_Conv(RealType *P, const RealType *pA, const RealType *pB,const int m,const int ma,const int mb){
	int ma1 = ma-1;
	for(int i=min(m,ma1+mb)-1;i>=ma1;i--){
		P[ma1+i*ma] = pA[ma1]*pB[i-ma1];
		for(int j=ma1-1;j>=0;j--){
			P[j+i*ma] = P[j+1+i*ma]+pA[j]*pB[i-j];
		}
	}
	for(int i=ma1-1;i>=0;i--){
		for(int j=i;j>=0;j--){
			P[j+i*ma] = P[j+1+i*ma]+pA[j]*pB[i-j];
		}
	}	
}

template <typename RealType>
void SingleCensoring(RealType *pAc, RealType *mA, RealType *mB, RealType *mC, const RealType *pA, const RealType *pB,
										 const RealType *pC, const RealType *SB, const RealType *pApC, const int x, const int y,
										 const RealType px, const int ma, const int mb, const int mc, const int &counts){
	if(px>0){
		int Myc = max(0,y-mc+1), mya = min(y+1,ma);
		
		for(int j=Myc;j<mya;j++){
			pAc[j] = pA[j]*pC[y-j]/px;
		}
		for(int j=max(Myc,x-mb+1);j<min(mya,x+1);j++){
			pAc[j] *= SB[x-j];
		}
    for(int j=Myc;j<mya;j++){
    	mA[j] += counts*pAc[j];
    }
    for(int j=0;j<min(mb,x+1);j++){
    	mB[j] += counts*pB[j]/px*pApC[x-j+1+y*ma]; 
    }
    for(int j=x+1;j<mb;j++){
    	mB[j] += counts*pB[j]/px*pApC[y*ma]; 
    }    
    for(int j=Myc;j<mya;j++){
    	mC[y-j] += counts*pAc[j];
    }
	}		
}

template <typename RealType>
void DoubleCens_Conv(RealType *P, const RealType *pA, const RealType *SB,const int m,const int ma,const int mb){
	int ma1 = ma-1;
	for(int i=min(m,ma1+mb)-1;i>=ma1;i--){
		P[ma1+i*ma] = pA[ma1]*SB[i-ma1];
		for(int j=ma1-1;j>=0;j--){
			P[j+i*ma] = P[j+1+i*ma]+pA[j]*SB[i-j];
		}
	}
	for(int i=ma1-1;i>=0;i--){
		P[ma1+i*ma] = pA[ma1];	
		for(int j=ma1-1;j>i;j--){
			P[j+i*ma] = P[j+1+i*ma]+pA[j];
		}		
		for(int j=i;j>=0;j--){
			P[j+i*ma] = P[j+1+i*ma]+pA[j]*SB[i-j];
		}
	}	
}

template <typename RealType>
void DoubleCensoringBC(RealType *mB, const RealType *pB, const RealType *pASC, const int x, const int y, const int ma,
											 const int mb, const RealType p, const int &counts){
  for(int j=0;j<min(mb,x+1);j++){
  	mB[j] += counts*pB[j]/p*pASC[x-j+1+y*ma];
  }
  for(int j=x+1;j<mb;j++){
  	mB[j] += counts*pB[j]/p*pASC[y*ma];
  }  
}

template <typename RealType>
void DoubleCensoring(RealType *pAc, RealType *mA, RealType *mB, RealType *mC, const RealType *pA, const RealType *pB, const RealType *pC,
									   const RealType *SB, const RealType *SC, const RealType *pASB, const RealType *pASC, const int x,
									   const int y, const RealType p, const int ma, const int mb, const int mc, const int &counts){
	if(p>0){
		int mya = min(y+1,ma), mxa = min(x+1,ma);
	
		for(int j=0;j<ma;j++){
			pAc[j] = pA[j]/p;
		}
		for(int j=max(0,x-mb+1);j<mxa;j++){
			pAc[j] *= SB[x-j];
		}
		for(int j=max(0,y-mc+1);j<mya;j++){
			pAc[j] *= SC[y-j];
		}
		for(int j=0;j<ma;j++){
			mA[j] += counts*pAc[j];
		}		
		DoubleCensoringBC(mB,pB,pASC,x,y,ma,mb,p,counts);
		DoubleCensoringBC(mC,pC,pASB,y,x,ma,mc,p,counts);
	}
}

template <typename RealType>
void pABC_cond(const int &mx, const int &my, RealType* pAco, RealType* pBco, RealType* pCco, const RealType* pA,
							 const RealType* pB, const RealType* pC, const RealType* SB, const RealType* SC, const RealType &pXY){
	RealType minXY = min(mx,my);
	for(int a=0;a<=minXY;a++){
		pAco[a] = pA[a]*(1.-SB[mx-a])*(1.-SC[my-a])/pXY;
	}
	for(int b=0;b<=mx;b++){
		pBco[b] = 0.;
		for(int a=0;a<=min(mx-b,my);a++){
			pBco[b] += pA[a]*(1.-SC[my-a]);
		}
		pBco[b] *= pB[b]/pXY;
	}
	for(int c=0;c<=my;c++){
		pCco[c] = 0.;
		for(int a=0;a<=min(my-c,mx);a++){
			pCco[c] += pA[a]*(1.-SB[mx-a]);
		}
		pCco[c] *= pC[c]/pXY;
	}	
}

/* Float instantiation */
template void kernel_Compute_pXY(float *pABC, float *P, const float *pA, const float *pB, const float *pC,
																 const int i, const int j, const int m, const int mxy2, const int lmin, const int lmax);

template void Compute_pXY(float *pABC,float *P,const int m,const int ma,const int mb,const int mc,const float *pA,
			    								const float *pB,const float *pC, const int minX, const int minY);

template void Compute_pTrunc(float *P,const int mtx,const int mty,const int meps,const int eps0,const float *pT,
									  				 const float *ST, const float *Se, const float cumTx=1.);

template void Convolution(float *P,const int m,const int ma,const int mb,const float *pA,const float *pB);

template void SingleCens_Conv(float *P, const float *pA, const float *pB,const int m,const int ma,const int mb);

template void SingleCensoring(float *pAc, float *mA, float *mB, float *mC, const float *pA, const float *pB,
										 					const float *pC, const float *SB, const float *pApC, const int x, const int y,
										 					const float px, const int ma, const int mb, const int mc, const int &counts);

template void DoubleCens_Conv(float *P, const float *pA, const float *SB,const int m,const int ma,const int mb);

template void DoubleCensoringBC(float *mB, const float *pB, const float *pASC, const int x, const int y, const int ma,
											 					const int mb, const float p, const int &counts);

template void DoubleCensoring(float *pAc, float *mA, float *mB, float *mC, const float* pA, const float* pB, const float *pC,
									   					const float *SB, const float *SC, const float *pASB, const float *pASC, const int x,
									   					const int y, const float p, const int ma, const int mb, const int mc, const int &counts);
									   					
template void pABC_cond(const int &mx, const int &my, float* pAco, float* pBco, float* pCco, const float* pA,
							 const float* pB, const float* pC, const float* SB, const float* SC, const float &pXY);									   					


/* Double instantiation */
template void kernel_Compute_pXY(double *pABC, double *P, const double *pA, const double *pB, const double *pC,
																 const int i, const int j, const int m, const int mxy2, const int lmin, const int lmax);

template void Compute_pXY(double *pABC,double *P,const int m,const int ma,const int mb,const int mc,const double *pA,
			    								const double *pB,const double *pC, const int minX, const int minY);

template void Compute_pTrunc(double *P,const int mtx,const int mty,const int meps,const int eps0,const double *pT,
									  				 const double *ST, const double *Se, const double cumTx=1.);

template void Convolution(double *P,const int m,const int ma,const int mb,const double *pA,const double *pB);

template void SingleCens_Conv(double *P, const double *pA, const double *pB,const int m,const int ma,const int mb);

template void SingleCensoring(double *pAc, double *mA, double *mB, double *mC, const double *pA, const double *pB,
										 					const double *pC, const double *SB, const double *pApC, const int x, const int y,
										 					const double px, const int ma, const int mb, const int mc, const int &counts);

template void DoubleCens_Conv(double *P, const double *pA, const double *SB,const int m,const int ma,const int mb);

template void DoubleCensoringBC(double *mB, const double *pB, const double *pASC, const int x, const int y, const int ma,
											 					const int mb, const double p, const int &counts);

template void DoubleCensoring(double *pAc, double *mA, double *mB, double *mC, const double* pA, const double* pB, const double *pC,
									   					const double *SB, const double *SC, const double *pASB, const double *pASC, const int x,
									   					const int y, const double p, const int ma, const int mb, const int mc, const int &counts);

template void pABC_cond(const int &mx, const int &my, double* pAco, double* pBco, double* pCco, const double* pA,
							 const double* pB, const double* pC, const double* SB, const double* SC, const double &pXY);									   					

