extern "C"
{
#include <ccomplex>
};


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
	ComplexType e1, m1;  //e1,m1 - относительные проницаемости
	ComplexType cn[3];
	ComplexType cdr3[3];
	ComplexType pr1[3];

	ElCondType::TYPE tol; // tol - признак (0-металл, -1 - диэлектрик)
	T z[3][3];
	T r0[3];

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

	T rr1[3];
	VecType r0t;

	VecType r0p;
	T dr3[3];
	T d(0);

	ComplexType f[3];
	T ad(0);
	T k0(0);

	VecType r1; // TODO: Требуется правильная инициализация
	VecType rr; // TODO: Требуется правильная инициализация
	VecType p0; // TODO: Требуется правильная инициализация
	VecType p1; // TODO: Требуется правильная инициализация
	
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
	
	/*
	! Вычисление полного рассеянного поля в вершинах треугольника 
	do i2=1,3 ! НОМЕР ТОЧКИ ТРЕУГОЛЬНИКА	
	*/
	PointType n;
	for (int i = 0; i < 3 ; i++)
	{
		/*
		do i6=1,3! КООРДИНАТА ТРЕУГОЛЬНИКА(X=1),(Y=2), (Z=3)  
			n(i6)=n0(i2,i6)
		end do
		tes=0
		*/		
		n = tr[i];

		/*
		do i6=1,3 ! НОМЕР ТОЧКИ ТРЕУГОЛЬНИКА
			tes=tes+p0(i6)*n(i6)
		end do
		p0t=p0-n*tes
		tes=0
		*/
		tes = (p0*n).ManhattanLength();
		VecType p0t = (p0 - n*tes);

		/*
		do i6=1,3 ! НОМЕР ТОЧКИ ТРЕУГОЛЬНИКА
			tes=tes+rr(i6)*n(i6)
		end do
		a=0
		*/
		tes = (rr*n).ManhattanLength();
		a = 0;

		/*
		if (abs(tes).gt.0.001) then
		*/
		if(abs(tes) > std::numeric_limits<float>::epsilon())
		{
			/*
			tes2=tes
			c0=csqrt(1.-(1-tes2*tes2)/(e1*m1))
			c2=c0
			*/
			T tes2 = tes;

			ComplexType sqrtRoot(T(1)-(T(1)-tes2*tes2));
			sqrtRoot = sqrtRoot/(e1*m1);

			c0 = std::sqrt(sqrtRoot);
			c2 = c0;

			/*
			! Расчет поля для металлического треугольника
			*/
			//if (tol.ge.0.0) then 
			if(tol == ElCondType::Metall)
			{
				//c1=k0*csqrt(e1*m1)*tol*c0;
				c1 = k0*std::sqrt(e1*m1)*ElCondType::Metall*c0; // Ноль?!												
				
				//c=csqrt(m1/e1)*c0*cdsin(c1)/cdcos(c1)
				c = std::sqrt(m1/e1)*c0*std::sin(c1)/std::cos(c1);	
				
				//c0=c*tes2
				c0 = c*tes2;

				//c1=(ji*c0-1)/(ji*c0+1)
				c1 = (ji*c0)/(ji*c0+1);

				//call wekt(n,rr,r0p,amod) 
				r0p = VecType(n).CrossProduct(rr);
				T amod = r0p.Length();

				//r0t=rr-n*tes2
				r0t = rr - n*tes2;

//				p1t=c1*p0t+2*ji*c/(ji*c0+1)*(r0t*(r0t(1)*p0(1)+r0t(2)*p0(2)
//#     +r0t(3)*p0(3))/(ji*c+tes2)+r0p*(r0p(1)*p0(1)+r0p(2)*p0(2)
//#     +r0p(3)*p0(3))/e1/m1/(ji*c+c2*c2/tes2))

				// (r0t*(r0t(1)*p0(1)+r0t(2)*p0(2)+r0t(3)*p0(3))
				VecType tmpMulti(r0t*(r0t*p0).ManhattanLength());
				p1t  = c1*p0t;
				p1t += 2*ji*c/(ji*c0+1)*((tmpMulti)/(ji*c+tes2)+tmpMulti/e1/m1/(ji*c+c2*c2/tes2));
			}
		}
	}

	return;
}