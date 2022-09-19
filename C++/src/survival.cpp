// survival related functions

#include "survival.h"

// Compute the survival distribution from the pdf
template <typename RealType>
void surv_pdf(const int &n, RealType *S, const RealType *p, const RealType &S0){
	S[0] = S0 - p[0];
	for(int i=1;i<n;i++){
		S[i] = S[i-1] - p[i];
	}
}

// Compute the pdf from the survival distribution
template <typename RealType>
void pdf_surv(const int &n, RealType *p, const RealType *S){
	p[0] = 1.-S[0];
	for(int i=1;i<n;i++){
		p[i] = S[i-1]-S[i];
	}
}

// Compute bivariate survival from bivariate pdf
template <typename RealType>
void surv_pdf2(const int &n, RealType *S, const RealType *P, const RealType *PX, const RealType* PY){
	RealType sumy = P[0]-PX[0], Sy = 1.-PY[0];
	S[0] = Sy+sumy;
	for(int y=1;y<n;y++){
		sumy += P[y];
		Sy   -= PY[y];
		S[y]  = Sy+sumy;
	}
	for(int x=1;x<n;x++){
		RealType sumy = -PX[x];
		for(int y=0;y<n;y++){
			sumy    += P[y+x*n];
			S[y+x*n] = S[y+(x-1)*n]+sumy;
		}
	}
}

template <typename RealType>
void pdf2_marginals(const int &n, RealType *pX, RealType* pY, const RealType* P){
	for(int x=0;x<n;x++){
		pX[x] = 0.;
		for(int y=0;y<n;y++){
			pX[x] += P[y+x*n];
		}
	}
	for(int y=0;y<n;y++){
		pY[y] = 0.;
		for(int x=0;x<n;x++){
			pY[y] += P[y+x*n];
		}
	}	
}

template <typename RealType>
void compute_FXSY(const int &m, RealType* FXSY, const RealType* pX, const RealType* P){
	FXSY[0] = pX[0]-P[0];
	for(int y=1;y<m;y++){
		FXSY[y] = FXSY[y-1]-P[y];
	}
	for(int x=1;x<m;x++){
		RealType pXSY = pX[x];
		for(int y=0;y<m;y++){
			pXSY       -= P[y+x*m];		
			FXSY[y+x*m] = FXSY[y+(x-1)*m]+pXSY;
		}
	}
}

template <typename RealType>
void compute_SXFY(const int &m, RealType* SXFY, const RealType* pY, const RealType* P){
	SXFY[0] = pY[0]-P[0];
	for(int x=1;x<m;x++){
		SXFY[x*m] = SXFY[(x-1)*m]-P[x*m];
	}
	for(int y=1;y<m;y++){
		RealType SXpY = pY[y];
		for(int x=0;x<m;x++){
			SXpY       -= P[y+x*m];		
			SXFY[y+x*m] = SXFY[y-1+x*m]+SXpY;
		}
	}
}

/* Float instantiation */
template void surv_pdf(const int &n, float *S, const float *p, const float &S0);
template void pdf_surv(const int &n, float *p, const float *S);
template void surv_pdf2(const int &n, float *S, const float *P, const float *PX, const float* PY);
template void pdf2_marginals(const int &n, float *pX, float* pY, const float* P);
template void compute_FXSY(const int &m, float* FXSY, const float* pX, const float* P);
template void compute_SXFY(const int &m, float* SXFY, const float* pY, const float* P);

/* Double instantiation */
template void surv_pdf(const int &n, double *S, const double *p, const double &S0);
template void pdf_surv(const int &n, double *p, const double *S);
template void surv_pdf2(const int &n, double *S, const double *P, const double *PX, const double* PY);
template void pdf2_marginals(const int &n, double *pX, double* pY, const double* P);
template void compute_FXSY(const int &m, double* FXSY, const double* pX, const double* P);
template void compute_SXFY(const int &m, double* SXFY, const double* pY, const double* P);
