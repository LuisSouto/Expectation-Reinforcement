/* Likelihood functions */

#include "likelihood_functions.h"

template <typename RealType>
RealType loglike_ltrc_2D(const int &n, const int* X, const int* Y, const int* TX, const int* TY,
									 		 const int* dx, const int* dy, const int &m, const RealType* S){
	RealType loglikelihood = 0.;									 
	for(int i=0;i<n;i++){
		loglikelihood -= log((S[Y[i]-1+(X[i]-1)*m]-dx[i]*S[Y[i]-1+X[i]*m]-dy[i]*S[Y[i]+(X[i]-1)*m]+dx[i]*dy[i]*S[Y[i]+X[i]*m])
											   /S[TY[i]-1+(TX[i]-1)*m]);
	}
	
	return loglikelihood/n;					 
}

template <typename RealType>
RealType loglike_rc_2D(const int &n, const int* X, const int* Y, const int* dx, const int* dy,
										 const int &m, const RealType* S){
	RealType loglikelihood = 0.;									 
	for(int i=0;i<n;i++){
		loglikelihood -= log(S[Y[i]-1+(X[i]-1)*m]-dx[i]*S[Y[i]-1+X[i]*m]-dy[i]*S[Y[i]+(X[i]-1)*m]+dx[i]*dy[i]*S[Y[i]+X[i]*m]);
	}
	
	return loglikelihood/n;
}


/* Float instantiation */
template float loglike_ltrc_2D(const int &n, const int* X, const int* Y, const int* TX, const int* TY,
									 		 				 const int* dx, const int* dy, const int &m, const float* S);
									 		 						
template float loglike_rc_2D(const int &n, const int* X, const int* Y, const int* dx, const int* dy,
										 				 const int &m, const float* S);									 		 						
										 				 
/* Double instantiation */
template double loglike_ltrc_2D(const int &n, const int* X, const int* Y, const int* TX, const int* TY,
									 		 				 const int* dx, const int* dy, const int &m, const double* S);
									 		 						
template double loglike_rc_2D(const int &n, const int* X, const int* Y, const int* dx, const int* dy,
										 				 const int &m, const double* S);									 		 																 				 
