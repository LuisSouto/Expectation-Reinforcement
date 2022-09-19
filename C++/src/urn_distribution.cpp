// Urn distribution 

#include <urn_distribution.h>

template <typename RealType>
void urn_init(const int &n, RealType *B, RealType *N, const RealType *P, const RealType &c){
	B[0] = c*P[0];
	N[0] = c;
	for(int i=1;i<n;i++){
		B[i] = c*P[i];
		N[i] = N[i-1]-B[i-1];
	}
}

template <typename RealType>
void urn_pdf(const int &n, RealType *P, const RealType *B, const RealType *N){

	P[0] = 1.;
	for(int i=1;i<n;i++)
		P[i] = (1.-B[i-1]/N[i-1])*P[i-1];
		
	for(int i=0;i<n;i++){
		if(N[i]==0){
			P[i] = 0.;
		}
		else{
			P[i] *= B[i]/N[i];
		}
	}
}

template <typename RealType>
void urn_surv(const int &n, RealType *P, const RealType *B, const RealType *N, const RealType &P0, const RealType &stp){

	P[0] = (1.-B[0]/N[0])*P0;
	for(int i=1;i<n;i++){
		if(N[i]<=1e-16){
			P[i] = 0.;
		}
		else{
			P[i] = (1.-B[i]/N[i])*P[i-1];
//			if(P[i-1]-P[i]>stp) P[i] = P[i-1];
		}
	}
}

template <typename RealType>
void urn_surv_trunc(const int &n, RealType *P, const int &maxT, const RealType *B, const RealType *N){

	P[maxT] = 1.;
	for(int i=n-1;i>=0;i--){
		if(i>=maxT){
			P[i] = 1.;
		}
		else{
			P[i] = (N[i+1]-B[i+1])/N[i+1]*P[i+1];
		}
	}
	for(int i=0;i<n;i++) P[i] = 1. - P[i];
}

template <typename RealType>
void update_urn(const int &n, RealType *B, RealType *N, const RealType *Bo, const RealType *No, const RealType *m,
								const RealType *s, const RealType &r){
	for(int i=0;i<n;i++){
		B[i] = Bo[i]+r*m[i];
		N[i] = No[i]+r*(m[i]+s[i]);
	}
}

/* Float instantiation */
template void urn_init(const int &n, float *B, float *N, const float *P, const float &c);
template void urn_pdf(const int &n, float *P, const float *B, const float *N);
template void urn_surv(const int &n, float *P, const float *B, const float *N, const float &P0, const float &stp);
template void urn_surv_trunc(const int &n, float *P, const int &maxT, const float *B, const float *N);
template void update_urn(const int &n, float *B, float *N, const float *Bo, const float *No, const float *m,
								         const float *s, const float &r);
								         
/* Double instantiation */
template void urn_init(const int &n, double *B, double *N, const double *P, const double &c);
template void urn_pdf(const int &n, double *P, const double *B, const double *N);
template void urn_surv(const int &n, double *P, const double *B, const double *N, const double &P0, const double &stp);
template void urn_surv_trunc(const int &n, double *P, const int &maxT, const double *B, const double *N);
template void update_urn(const int &n, double *B, double *N, const double *Bo, const double *No, const double *m,
								         const double *s, const double &r);								         
