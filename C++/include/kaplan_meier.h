#ifndef kaplan_meier_h__
#define kaplan_meier_h__

#include <algorithm>

using namespace std;

/* Kaplan-Meier functions (header) */

template <typename RealType>
RealType mt(const int n,const int *a,const int a_n);

template <typename RealType>	
RealType rt(const int n,const int *a,const int a_n);

template <typename RealType>
RealType st(const int n,const int *a,const int a_n);

template <typename RealType>
RealType st(const int n,const int *a,const int *t,const int a_n);

template <typename RealType>
RealType mtc(const int n,const int *a,const int *d,const int a_n);

void mt(int *m,const int nd,const int n,const int *a);

void rt(int *m,const int nd,const int n,const int *a);

void lt(int *d,const int nd,const int n,const int *a);

void st(int *m,const int nd,const int n,const int *a);

void rt(int *m,const int nd,const int n,const int *a,const int *t);

void st(int *m,const int nd,const int n,const int *a,const int *t);

void mtc(int *m,const int nd,const int n,const int *a,const int *d);

template <typename RealType>
void surv_km(RealType *Sx,const int mx,const int n,const int *X);

template <typename RealType>
void surv_km(RealType *Sx,const int mx,const int n,const int *X,const int *d);

template <typename RealType>
void km_ltrc(const int &m, RealType *Sx, const int &n, const int *X, const int *T, const int *d);

#endif // kaplan_meier_h__
