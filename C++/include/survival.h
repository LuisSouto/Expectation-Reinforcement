#ifndef survival_h__
#define survival_h__

#include <iostream>

using namespace std;

template <typename RealType>
void surv_pdf(const int &n, RealType *S, const RealType *p, const RealType &S0=1.);

template <typename RealType>
void pdf_surv(const int &n, RealType *p, const RealType *S);

template <typename RealType>
void surv_pdf2(const int &n, RealType *S, const RealType *P, const RealType *PX, const RealType* PY);

template <typename RealType>
void pdf2_marginals(const int &n, RealType *pX, RealType* pY, const RealType* P);

template <typename RealType>
void compute_FXSY(const int &m, RealType* FXSY, const RealType* pX, const RealType* P);

template <typename RealType>
void compute_SXFY(const int &m, RealType* SXFY, const RealType* pY, const RealType* P);

#endif // survival_h__
