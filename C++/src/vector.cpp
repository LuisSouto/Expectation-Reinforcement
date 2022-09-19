// Basic pointer operations

#include "vector.h"

template<typename RealType>
void zeros(const int &n, RealType* v){
	for(int i=0;i<n;i++) v[i] = 0.;
}

template<typename RealType>
RealType sum(const int &n,const RealType* v){
	RealType s = 0.;
	for(int i=0;i<n;i++) s+= v[i];
	return s;
}

template<typename RealType>
void mult(const int &n, RealType* v, const RealType &c){
	for(int i=0;i<n;i++) v[i] *= c;
}

template<typename RealType>
void div(const int &n, RealType* v, const RealType &c){
	for(int i=0;i<n;i++) v[i] /= c;
}

template<typename RealType>
void add_u_to_v(const int &n, RealType* v, const RealType* u, const RealType &c){
	for(int i=0;i<n;i++) v[i] += c*u[i];
}

/* Integer instantiation */
template int sum(const int &n, const int* v);
template void zeros(const int &n, int* v);

/* Float instantiation */
template void zeros(const int &n, float* v);
template float sum(const int &n, const float* v);
template void mult(const int &n, float* v, const float &c);
template void div(const int &n, float* v, const float &c);
template void add_u_to_v(const int &n, float* v, const float* u, const float &c);


/* Double instantiation */
template void zeros(const int &n, double* v);
template double sum(const int &n, const double* v);
template void mult(const int &n, double* v, const double &c);
template void div(const int &n, double* v, const double &c);
template void add_u_to_v(const int &n, double* v, const double* u, const double &c);
