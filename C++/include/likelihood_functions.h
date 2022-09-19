
#ifndef likelihood_functions_h__
#define likelihood_functions_h__

#include <math.h>

template <typename RealType>
RealType loglike_ltrc_2D(const int &n, const int* X, const int* Y, const int* TX, const int* TY,
									 		 const int* dx, const int* dy, const int &m, const RealType* S);
									 		 
template <typename RealType>									 		 
RealType loglike_rc_2D(const int &n, const int* X, const int* Y, const int* dx, const int* dy,
										 const int &m, const RealType* S);
									 
#endif // likelihood_functions_h__									 
