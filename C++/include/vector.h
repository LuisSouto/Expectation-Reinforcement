#ifndef vector_h__
#define vector_h__


template <typename RealType>
void zeros(const int &n, RealType *v);

template <typename RealType>
RealType sum(const int &n,const RealType *v);

template<typename RealType>
void mult(const int &n, RealType* v, const RealType &c);

template<typename RealType>
void div(const int &n, RealType* v, const RealType &c);

template <typename RealType>
void add_u_to_v(const int &n, RealType *v, const RealType *u, const RealType &c=1.);

#endif // vector_h__
