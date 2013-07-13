template<typename T>
Calc<T>::Calc()
{
	Init();
}

template<typename T>
Calc<T>::~Calc()
{

}

template<typename T>
void Calc<T>::Init()
{

}

template<typename T> 
void Calc<T>::FieldOfTriangle(const TriangleType& tr)
{
	ComplexType ji(0,1);
	ComplexType e1;
	ComplexType m1;
	ComplexType cn[3];
	ComplexType cdr3[3];
	ComplexType pr1[3];

	T tol(0);
	T z[3][3];
	T r0[3];
	T p0[3];
	T p1[3];
	T n0[3][3];

	ComplexType et[3];
	ComplexType ht[3];
	ComplexType a;
	ComplexType bi[3];
	ComplexType c;
	ComplexType c0;
	ComplexType c1;
	ComplexType p1t[3];
	ComplexType c2;
	ComplexType ep;
		
	ComplexType tes5;
	T a1(0);
	T a2(0);
	T na(0);
	T nb(0);
	T nc(0);
	T nd(0);
	T h(0);

	T p0t[3];
	T rr1[3];
	T r0t[3];

	T r0p[3];
	T dr3[3];
	T d(0);
	T n[3];
	
	ComplexType f[3];
	T ad(0);
	T k0(0);
	T dd(0);
	T tes(0);

	VecType r1;
	VecType rr;

	VecType dr0(r1+rr);
	VecType dr1(tr.p2 - tr.p3);
	VecType dr2(tr.p1 - tr.p3);
	
	a1 = -k0*(dr0*dr1).ManhattanLength();
	a2 = -k0*(dr0*dr2).ManhattanLength();

	return;
}