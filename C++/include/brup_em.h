#ifndef brup_em_h__
#define brup_em_h__

/* Header for the BRUP-EM optimization routine */

#include <iostream>
#include <random>
#include <algorithm>
#include <fstream>
#include <string>
#include <memory>

#include "urn_distribution.h"
#include "survival.h"
#include "vector.h"
#include "kaplan_meier.h"
#include "censoring_functions.h"
#include "likelihood_functions.h"

using namespace std;


template <typename RealType>
void BRUP_EMopt(const int &m, RealType *pA, RealType *pB, RealType *pC, RealType *pTx, RealType *pE,
								const RealType *pA0, const RealType *pB0, const RealType *pC0, const RealType *pTx0, 
								const RealType *pE0, const int &n, const int *X, const int *Y, const int *TX,
								const int *TY, const int *dx, const int *dy, const int *eps, const RealType &cA,
								const RealType &cB, const RealType &cC, const RealType &cT, const RealType &cE,
								const RealType &r, const int &N=5000, const bool &is_trunc=true, const RealType &tol=1.e-4);
								
template<typename RealType>
void BRUP_ApplyPrior(const int &m, RealType* pA, RealType* pB, RealType *pC, RealType* pTx, RealType* pE,
										 const RealType* pAe, const RealType* pBe, const RealType* pCe, const RealType* pTxe,
										 const RealType* pEe, const RealType* pAp, const RealType* pBp, const RealType* pCp,
										 const RealType* pTxp, const RealType* pEp, const RealType &cA, const RealType &cB,
										 const RealType &cC, const RealType &cT, const RealType &cE, const RealType &r);								
								
#endif // brup_em_h__
