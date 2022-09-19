/* Main script. Calls the Expectation-Reinforcement routine to compute the nonparametric */
/* probability distribution under bivariate left-truncation and right-censoring.         */

#include <iostream>
#include <random>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <sys/time.h>
#include <string>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/normal.hpp>
#include <typeinfo>

#include "brup_em.h"

#define RealType double

using namespace std;
using namespace boost;

int find_maxlim(RealType *P,const int m,const RealType thres=1e-12){
	for(int i=m-1;i>=0;i--){
		if(P[i]>thres){
			return i;
		}
	}
	return 0;
}

int main(){

/*-------------------------------------------------------------------*/
/*----------------------DATA INITIALIZATION--------------------------*/
/*-------------------------------------------------------------------*/

	/* Input and output files */
	
	string filename = "poisson", nr = "";
	ifstream fin("../data/"+filename+"_ltrc.dat");
	ofstream fout("../results/"+filename+nr+"_ltrc_res_l.dat");
	ofstream fout_h("../results/"+filename+nr+"_ltrc_res_h.dat");	

	struct timeval start, end;   // Time profiling 
	
	// Sample size
	int n = 0;
	
	/* Ignores lines with text until it finds a line with a single number */
	string str; int istring = 0;
 	while(istring==0){
		getline(fin,str);
		istringstream iss(str);
		iss >> istring;
	}
	n = istring;   // Number of observations
	
  int *X   = new int[n], *Y  = new int[n], *TX = new int[n], *TY = new int[n];  // Target and truncation variables
	int *dx  = new int[n], *dy = new int[n];                                      // Censoring indicators
	int *eps = new int[n];                                                        // Difference in age
	
  /* Read input data */
	for(int i=0;i<n;i++){
		if(filename=="UOFM"){ // UOFM is the name we have given to the empirical data file. Change accordingly!
			fin>>X[i]>>Y[i]>>TX[i]>>TY[i]>>dx[i]>>dy[i];
			eps[i] = TY[i]-TX[i];
		}	
		else{
			fin>>X[i]>>Y[i]>>TX[i]>>TY[i]>>dx[i]>>dy[i]>>eps[i];
		}
	}
	fin.close();

	/* Select upper bounds for each variable */	
	int mx = *max_element(X,X+n), my = *max_element(Y,Y+n), m = max(mx,my)+10;

/*-------------------------------------------------------------------*/
/*----------------------CALL ER ROUTINE------------------------------*/
/*-------------------------------------------------------------------*/
	
	/* EM priors and settings */
	RealType *pAe = new RealType[m], *pBe = new RealType[m], *pCe = new RealType[m], *pTxe = new RealType[m], *pEe = new RealType[m];	
	RealType *pA0 = new RealType[m], *pB0 = new RealType[m], *pC0 = new RealType[m], *pTx0 = new RealType[m], *pE0 = new RealType[m];
	RealType *pAl = new RealType[m], *pBl = new RealType[m], *pCl = new RealType[m], *pTxl = new RealType[m], *pEl = new RealType[m];
	RealType *pAh = new RealType[m], *pBh = new RealType[m], *pCh = new RealType[m], *pTxh = new RealType[m], *pEh = new RealType[m];
	
	const RealType la = 20., lb = 20., lc = 20., ltx = 50, leps = 10;  // Poisson parameters 
	int N = 10000;                                                     // Max number of iterations
	RealType c = 1.e0, cA = c, cB = c, cC = c, cT = c, r = 1.e4;       // Belief and reinforcement parameters

	math::poisson_distribution<> dist_A(la), dist_B(lb), dist_C(lc);
	math::poisson_distribution<> dist_Tx(ltx), dist_E(leps);

	for(int i=0;i<m;i++){
		pA0[i]  = math::pdf(dist_A,i);
		pB0[i]  = math::pdf(dist_B,i);
		pC0[i]  = math::pdf(dist_C,i);
		pTx0[i] = math::pdf(dist_Tx,i);
		pE0[i]  = math::pdf(dist_E,i);
	}
		
	gettimeofday(&start, NULL); // Start timing
	
	/* Low belief scenario */
	BRUP_EMopt(m,pAe,pBe,pCe,pTxe,pEe,pA0,pB0,pC0,pTx0,pE0,n,X,Y,TX,TY,dx,dy,eps,cA,cB,cC,cT,cT,r,N,true,(RealType)1e-9);
	BRUP_ApplyPrior(m,pAl,pBl,pCl,pTxl,pEl,pAe,pBe,pCe,pTxe,pEe,pA0,pB0,pC0,pTx0,pE0,cA,cB,cC,cT,cT,r);
	
	/* High belief scenario */
	c  = 1.e1;
	cA = c;
	cB = c;
	cC = c;
	cT = c;
	r  = 1.e0;
	BRUP_ApplyPrior(m,pAh,pBh,pCh,pTxh,pEh,pAe,pBe,pCe,pTxe,pEe,pA0,pB0,pC0,pTx0,pE0,cA,cB,cC,cT,cT,r);	
	
	gettimeofday(&end, NULL); // End timing
	
  float delta = ((end.tv_sec-start.tv_sec)*1e6+end.tv_usec-start.tv_usec)/1.e6;
         
  cout<<"Elapsed time: "<<delta<<" s."<<endl;
 
  
	/*------------------------------------------------------------------------------------------------------------*/
	/*-----------------------------------PRINT RESULTS TO OUTPUT FILE---------------------------------------------*/
	/*------------------------------------------------------------------------------------------------------------*/
	
	fout<<m<<endl;
  for(int i=0;i<m;i++){
  	fout<<pAl[i]<<" "<<pBl[i]<<" "<<pCl[i]<<" "<<pTxl[i]<<" "<<pEl[i]<<endl;
  }
  
	fout_h<<m<<endl;
  for(int i=0;i<m;i++){
  	fout_h<<pAh[i]<<" "<<pBh[i]<<" "<<pCh[i]<<" "<<pTxh[i]<<" "<<pEh[i]<<endl;
  }
  
  fout.close();
  fout_h.close();

	return 0;
}
