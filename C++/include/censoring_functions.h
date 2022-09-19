#ifndef censoring_functions_h__
#define censoring_functions_h__

/* Headers for censoring routines related to the B-RUP calibration via EM */

#include <iostream>

using namespace std;

void find_unique_LTRC(const int &n, int &nu, int* Xu, int *Yu, int *dxu, int *dyu,
								      const int* X, const int *Y, const int *dx, const int *dy);
								      
void count_unique_LTRC(const int &n, const int &nu, int* uni_counts, const int* Xu,
										   const int* Yu, const int* dxu, const int* dyu, const int* X,
										   const int* Y, const int* dx, const int* dy);								      

template <typename RealType>
void kernel_Compute_pXY(RealType *pABC, RealType *P, const RealType *pA, const RealType *pB, const RealType *pC,
												const int i, const int j, const int mxy, const int mxy2, const int lmin, const int lmax);

template <typename RealType>
void Compute_pXY(RealType *pABC,RealType *P,const int mxy,const int ma,const int mb,const int mc,const RealType *pA,
								 const RealType *pB,const RealType *pC, const int minX = 0, const int minY = 0);

template <typename RealType>
void Compute_pTrunc(RealType *P,const int mtx,const int mty,const int meps,const int eps0,const RealType *pT,
									  const RealType *ST, const RealType *Se, const RealType cumTx=1.);

template <typename RealType>
void Convolution(RealType *P,const int mxy,const int ma,const int mb,const RealType *pA,const RealType *pB);

template <typename RealType>
void SingleCens_Conv(RealType *P, const RealType *pA, const RealType *pB,const int mxy,const int ma,const int mb);

template <typename RealType>
void SingleCensoring(RealType *pAc, RealType *mA, RealType *mB, RealType *mC, const RealType *pA, const RealType *pB,
										 const RealType *pC, const RealType *SB, const RealType *pApC, const int x, const int y,
										 const RealType px, const int ma, const int mb, const int mc, const int &counts);

template <typename RealType>
void DoubleCens_Conv(RealType *P, const RealType *pA, const RealType *SB,const int mxy,const int ma,const int mb);

template <typename RealType>
void DoubleCensoringBC(RealType *mB, const RealType *pB, const RealType *pASC, const int x, const int y, const int ma,
											 const int mb, const RealType p, const int &counts);

template <typename RealType>
void DoubleCensoring(RealType *pAc, RealType *mA, RealType *mB, RealType *mC, const RealType* pA, const RealType* pB, const RealType *pC,
									   const RealType *SB, const RealType *SC, const RealType *pASB, const RealType *pASC, const int x,
									   const int y, const RealType p, const int ma, const int mb, const int mc, const int &counts);
									   
template <typename RealType>
void pABC_cond(const int &mx, const int &my, RealType* pAco, RealType* pBco, RealType* pCco, const RealType* pA,
							 const RealType* pB, const RealType* pC, const RealType* SB, const RealType* SC, const RealType &pXY);									   
									   
#endif // censoring_functions_h__									   
