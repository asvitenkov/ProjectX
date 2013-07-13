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
	T tes(0);

	VecType r1; // TODO: Требуется правильная инициализация
	VecType rr; // TODO: Требуется правильная инициализация
	
	/*
	da=(z(2,2)-z(3,2))*(z(1,3)-z(3,3))-(z(1,2)-z(3,2))*(z(2,3)-z(3,3)) ! da=z(,)*z(,)-z(,)*z(,)
	db=(z(2,3)-z(3,3))*(z(1,1)-z(3,1))-(z(1,3)-z(3,3))*(z(2,1)-z(3,1)) ! db=z(,)*z(,)-z(,)*z(,) 
	dg=(z(2,1)-z(3,1))*(z(1,2)-z(3,2))-(z(2,2)-z(3,2))*(z(1,1)-z(3,1)) ! dg=z(,)*z(,)-z(,)*z(,)
	dd=dsqrt(da*da+db*db+dg*dg)  
	*/
	T dd = tr.Square();

	/*
	a1=0
	a2=0
	do i6=1,3 ! КООРДИНАТА ТРЕУГОЛЬНИКА(X=1),(Y=2),(Z=3) 
		dr0(i6)=rr(i6)+r1(i6)   ! (Line 1155:) rr(i6)=-r1(i6)
		dr1(i6)=z(2,i6)-z(3,i6) ! массив из 3-х расстояний (по x,y,z) между 2 и 3 точками треугольника
		dr2(i6)=z(1,i6)-z(3,i6) ! массив из 3-х расстояний (по x,y,z) между 1 и 3 точками треугольника
		a1=a1+dr0(i6)*dr1(i6)   ! сумма (произведении по компанентам x,y,z)
		a2=a2+dr0(i6)*dr2(i6)   ! сумма (произведении по компанентам x,y,z)
	end do

	a1=-a1*k0
	a2=-a2*k0
	*/
	VecType dr0(r1+rr);
	VecType dr1(tr.p2 - tr.p3); // вектор треугольника
	VecType dr2(tr.p1 - tr.p3);
	
	a1 = -k0*(dr0*dr1).ManhattanLength();
	a2 = -k0*(dr0*dr2).ManhattanLength();

	VecType n0[3]; // Вектора нормалей к веришнам
	/*
	do i2=1,3 ! НОМЕР ТОЧКИ ТРЕУГОЛЬНИКА
		tes=0     ! АЗИМУТ
		do I6=1,3 ! КООРДИНАТА ТРЕУГОЛЬНИКА(X=1),(Y=2),(Z=3) 
			tes=tes+rr(i6)*n0(i2,i6)
		end do
		if (tes.gt.0) then
			do i6=1,3 ! КООРДИНАТА ТРЕУГОЛЬНИКА(X=1),(Y=2),(Z=3) 
				n0(i2,i6)=-n0(i2,i6)
			end do
		end if
	end do
	*/
	T tes;
	for (int i = 0; i < 3 ; i++) // перебираем вершины треугольника
	{
		tes = 0;
		tes = (rr*n0[i]).ManhattanLength();

		if(tes > std::numeric_limits<T>::epsilon())
		{
			n0[i] = n0[i].Revert(); 
		}
	}
	
	return;
}