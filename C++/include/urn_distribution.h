#ifndef urn_distribution_h__
#define urn_distribution_h__

// Urn distribution (header)

template <typename RealType>
void urn_init(const int &n, RealType *B, RealType *N, const RealType *P, const RealType &c);

template <typename RealType>
void urn_pdf(const int &n, RealType *P, const RealType *B, const RealType *N);

template <typename RealType>
void urn_surv(const int &n, RealType *P, const RealType *B, const RealType *N, const RealType &P0 = 1., const RealType &stp = 1.);

template <typename RealType>
void urn_surv_trunc(const int &n, RealType *P, const int &maxT, const RealType *B, const RealType *N);

template <typename RealType>
void update_urn(const int &n, RealType *B, RealType *N, const RealType *Bo, const RealType *No, const RealType *m,
								const RealType *s, const RealType &r);
								
#endif // urn_distribution_h__								
