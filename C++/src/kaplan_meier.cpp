/* Death and survivor functions for censored and non-censored data */

#include "kaplan_meier.h"

template <typename RealType>
RealType mt(const int n,const int *a,const int a_n){
	RealType d = 0.;
	for(int i=0;i<n;i++)
		d += (a[i]==a_n);
	
	return d;
}
	
void mt(int *d,const int nd,const int n,const int *a){
	for(int i=0;i<=nd;i++)
		d[i] = 0;
	for(int i=0;i<n;i++)
		d[a[i]]++;
}

template <typename RealType>	
RealType rt(const int n,const int *a,const int a_n){
	RealType d = 0.;
	for(int i=0;i<n;i++)
		d += (a[i]>a_n);
	
	return d;
}

void rt(int *d,const int nd,const int n,const int *a){
	for(int i=0;i<=nd;i++)
		d[i] = 0;
	for(int i=0;i<n;i++){
		for(int j=0;j<a[i];j++)
			d[j]++;
	}
}

void lt(int *d,const int nd,const int n,const int *a){
	for(int i=0;i<=nd;i++)
		d[i] = 0;
	for(int i=0;i<n;i++){
		for(int j=a[i]+1;j<=nd;j++)
			d[j]++;
	}
}

void rt(int *d,const int nd,const int n,const int *a,const int *t){
	for(int i=0;i<=nd;i++)
		d[i] = 0;
	for(int i=0;i<n;i++){
		for(int j=t[i];j<a[i];j++)
			d[j]++;
	}
}

template <typename RealType>
RealType st(const int n,const int *a,const int a_n){
	RealType d = 0.;
	for(int i=0;i<n;i++)
		d += (a[i]>=a_n);
	
	return d;
}

template <typename RealType>
RealType st(const int n,const int *a,const int *t,const int a_n){
	RealType d = 0.;
	for(int i=0;i<n;i++)
		d += ((a[i]>=a_n)&&(t[i]<=a_n));
		
	return d;
}

void st(int *d,const int nd,const int n,const int *a){
	for(int i=0;i<=nd;i++)
		d[i] = 0;
	for(int i=0;i<n;i++){
		for(int j=0;j<=a[i];j++){d[j]++;}
	}
}

void st(int *d,const int nd,const int n,const int *a,const int *t){
	for(int i=0;i<=nd;i++)
		d[i] = 0;
	for(int i=0;i<n;i++){
		for(int j=t[i];j<=a[i];j++)
			d[j]++;
	}
}

template <typename RealType>
RealType mtc(const int n,const int *a,const int *d,const int a_n){
	RealType m = 0.;
	for(int i=0;i<n;i++)
		m += (a[i]==a_n)&&(d[i]==1);
	
	return m;
}

void mtc(int *m,const int nd,const int n,const int *a,const int *d){
	for(int i=0;i<=nd;i++)
		m[i] = 0;	
	for(int i=0;i<n;i++)
		m[a[i]] += d[i];	
}

template <typename RealType>
void surv_km(RealType *Sx,const int mx,const int n,const int *X){
  Sx[0] = (1.-mt<RealType>(n,X,0)/st<RealType>(n,X,0));
  for(int i=1;i<=mx;i++){
  	if(st<RealType>(n,X,i)>0){
	  	Sx[i] = Sx[i-1]*(1.-mt<RealType>(n,X,i)/st<RealType>(n,X,i));
	  }
	  else{
	  	Sx[i] = 0.;
	  }
  }		
}

template <typename RealType>
void km_rc(RealType *Sx,const int mx,const int n,const int *X,const int *d){
  Sx[0] = (1.-mtc<RealType>(n,X,d,0)/st<RealType>(n,X,0));
  for(int i=1;i<=mx;i++){
  	if(st<RealType>(n,X,i)>0){
	  	Sx[i] = Sx[i-1]*(1.-mtc<RealType>(n,X,d,i)/st<RealType>(n,X,i));
	  }
	  else{
	  	Sx[i] = 0.;
	  }
  }		
}

template <typename RealType>
void km_ltrc(const int &m, RealType *Sx, const int &n, const int *X, const int *T, const int *d){
	RealType Ni   = st<RealType>(n,X,T,0);
	int      minT = *min_element(T,T+n);
	int      mx   = *max_element(X,X+n);
	
	if(Ni>0){
	  Sx[0] = (1.-mtc<RealType>(n,X,d,0)/st<RealType>(n,X,T,0));
	}
	else{
		Sx[0] = 1.;
	}
  for(int i=0;i<minT;i++) Sx[i] = 1.;
  
  for(int i=max(1,minT);i<=mx;i++){
  	Ni = st<RealType>(n,X,T,i);
  	if(Ni>0){
	  	Sx[i] = Sx[i-1]*(1.-mtc<RealType>(n,X,d,i)/Ni);
	  }
	  else{
	  	Sx[i] = 0.;
	  }
  }		
}

/* Float instantiation */
template float mt(const int n,const int *a,const int a_n);
template float rt(const int n,const int *a,const int a_n);
template float st(const int n,const int *a,const int a_n);
template float st(const int n,const int *a,const int *t,const int a_n);
template float mtc(const int n,const int *a,const int *d,const int a_n);
template void surv_km(float *Sx,const int mx,const int n,const int *X);
template void km_rc(float *Sx,const int mx,const int n,const int *X,const int *d);
template void km_ltrc(const int &m, float *Sx, const int &n, const int *X, const int *T, const int *d);

/* Double instantiation */
template double mt(const int n,const int *a,const int a_n);
template double rt(const int n,const int *a,const int a_n);
template double st(const int n,const int *a,const int a_n);
template double st(const int n,const int *a,const int *t,const int a_n);
template double mtc(const int n,const int *a,const int *d,const int a_n);
template void surv_km(double *Sx,const int mx,const int n,const int *X);
template void km_rc(double *Sx,const int mx,const int n,const int *X,const int *d);
template void km_ltrc(const int &m, double *Sx, const int &n, const int *X, const int *T, const int *d);
