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

/*
 * tr - Треугольник (матрица z)
 * tol - признак (0-металл, -1 - диэлектрик)
 **/
template<typename T> 
void Calc<T>::FieldOfTriangle(
	const TriangleType& tr, //!z  - координаты точек 
	const VecType& p0, //!P0,P1 - поляризации передатчика и приемника
	const VecType& p1, //!P0,P1 - поляризации передатчика и приемника
	const T& k0, //!k0 - волновое число
	const ComplexType& e1, //e1,m1 - относительные проницаемости 
	const ComplexType& m1, //e1,m1 - относительные проницаемости
	const ElCondType::TYPE tol,  //!tol - признак (0-металл, -1 - диэлектрик)
	ComplexType& ep// !ep - расчитанное поле
	) 
{
	T a1(0);
	T a2(0);
	T na(0);
	T nb(0);
	T nc(0);
	T nd(0);
	T h(0);
	T d(0);
	T ad(0);

	CVecType cn;
	CVecType cdr3;
	CVecType pr1;
	CVecType bi;
	CVecType p1t;
	CVecType f;
	CVecType et;
	CVecType ht;

	ComplexType ji(0,1);
	ComplexType a;
	ComplexType c;
	ComplexType c0;
	ComplexType c1;
	ComplexType c2;
	ComplexType tes5;

	//VecType r0; // не изпользуется
	VecType rr1;
	VecType r0t;
	VecType r0p;
	VecType dr3;
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

	VecType n0[3]; // Вектора нормалей к вершинам. Должны расчитываться в треугольнике
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
		if(tes > std::numeric_limits<float>::epsilon())
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
				c1 = std::sqrt(e1*m1)*c0*(k0*ElCondType::Metall); // Ноль?!												
				
				//c=csqrt(m1/e1)*c0*cdsin(c1)/cdcos(c1)
				c = std::sqrt(m1/e1)*c0*std::sin(c1)/std::cos(c1);	
				
				//c0=c*tes2
				c0 = c*tes2;

				//c1=(ji*c0-1)/(ji*c0+1)
				c1 = (ji*c0-T(1))/((ji*c0)+T(1));

				//call wekt(n,rr,r0p,amod) 
				r0p = VecType(n).CrossProduct(rr);
				T amod = r0p.Length();

				//r0t=rr-n*tes2
				r0t = rr - n*tes2;

//				p1t=c1*p0t+2*ji*c/(ji*c0+1)*(r0t*(r0t(1)*p0(1)+r0t(2)*p0(2)
//					+r0t(3)*p0(3))/(ji*c+tes2)+r0p*(r0p(1)*p0(1)+r0p(2)*p0(2)
//					+r0p(3)*p0(3))/e1/m1/(ji*c+c2*c2/tes2))

				// (r0t*(r0t(1)*p0(1)+r0t(2)*p0(2)+r0t(3)*p0(3))
				VecType tmp0(r0t*(r0t*p0).ManhattanLength()); // Вектор на число = вектор
				
				// комплексное число на вектор ?! 
				// На фортране это умножение комплексного числа на массив из 3 чисел - в геометрии такого нет
				//const ComplexType tmp1 = c1*p0t; 

				const ComplexType tmp2 = T(2)*ji*c;
				const ComplexType tmp3 = ji*c0+T(1);

				// Вектор делить на комплексное число. Или массив из 3 чисел делить на вектор.
				//const ComplexType tmp4 = tmp0/(ji*c+tes2);
				const ComplexType tmp5 = (ji*c+c2*c2/tes2);
				// Вообще не понтяно что тут делают.
				//const ComplexType tmp6 = tmp0/e1/m1/tmp5;

				// Массив из комплексных чисел - не известно что тут имеется ввиду.
				//p1t = tmp1 + tmp2 / (tmp3 * (tmp4 + tmp6));
			}
			else
			{
				//!Расчет поля для диэлектрического треугольника
				//c=csqrt(m1/e1)*c0
				c=std::sqrt(m1/e1)*c0;

				//c0=c*tes2
				c0 = c*tes2;
				// c1=(-c0-1)/(-c0+1)
				c1=(-c0-T(1))/(-c0+T(1));
				//call wekt(n,rr,r0p,amod)    ! векторное произведение двух векторов вх-n,rr; вых-r0p и amod
				r0p = VecType(n).CrossProduct(rr);
				T amod = r0p.Length();

				//r0t=rr-n*tes2
				r0t=rr-n*tes2;

				// Теже не понятные вычисления
				//p1t=c1*p0t+2*c/(-c0+1)*(r0t*(r0t(1)*p0(1)+r0t(2)*p0(2)
				//#     +r0t(3)*p0(3))/(c-tes2)+r0p*(r0p(1)*p0(1)+r0p(2)*p0(2)
				//#     +r0p(3)*p0(3))/e1/m1/(c-c2*c2/tes2))
			}

			//tes5=0
			//do i6=1,3
			//	tes5=tes5+rr(i6)*p1t(i6)
			//end do
			CVecType tmpRR = vector_cast<T, ComplexType>(rr);
			ComplexType tes5 = (tmpRR*p1t).ManhattanLength();

			//pr1=p1t+n*tes5/tes // координаты точки на комлексное число?!
			pr1 = p1t + vector_cast<T, ComplexType>(VecType(n)) * tes5 / tes;
			
			//rr1=rr-n*tes*2.0 
			rr1=rr-n*tes*T(2);
			
			//call wekt(n,p0,dr3,d)       ! векторное произведение двух векторов вх-n,p0; вых-dr3 и d.
			dr3 = VecType(n).CrossProduct(p0);
			T d = r0p.Length();
				
			//tes=0
			//do I6=1,3
			//	tes=tes+rr(i6)*z(3,i6)
			//end do
			tes = (rr*VecType(tr.p3)).ManhattanLength();

			//et=dr3
			et = vector_cast<T, ComplexType>(dr3);

			//cn=CMPLX(n)
			cn=CVecType();
			
			//call wektk(cn,pr1,cdr3)         !векторного произведения двух комплексных векторов вх-cn,pr1; вых-cdr3.
			cdr3 = cn.CrossProduct(pr1);

			//et=(et+cdr3)
			et=(et+cdr3);
			
			//dr3=0.
			//tes=0
			//do I6=1,3
			//	tes=tes+rr(i6)*n(i6)
			//end do
			// dr3=p0*tes
			dr3 = p0 * tes;
			tes = (rr*n).ManhattanLength();

			//tes=0
			//do I6=1,3
			//	tes=tes+p0(i6)*n(i6)
			//end do
			
			tes = (rr*n).ManhattanLength();

			//ht=dr3-rr*tes
			ht = vector_cast<T, ComplexType>(dr3) - vector_cast<T, ComplexType>(rr) * tes;

			//! Вычисление интеграла по трем точкам
			//dr3=0.
			//tes=0
			//do I6=1,3
			//	tes=tes+rr1(i6)*n(i6)
			//end do
			//cdr3=pr1*tes
			tes = (rr1*n).ManhattanLength();
			cdr3 = pr1 * tes;
			
			//tes5=0
			//do I6=1,3
			//	tes5=tes5+pr1(i6)*n(i6)
			//end do
			//cdr3=cdr3-rr1*tes5
			CVecType tmpN = vector_cast<T, ComplexType>(VecType(n));
			tes5 = (pr1*tmpN).ManhattanLength();
			cdr3 = cdr3 - vector_cast<T, ComplexType>(rr1) * tes5;

			//tes=0
			//do I6=1,3
			//	tes=tes+rr(i6)*z(3,i6)
			//end do
			tes = (rr*VecType(tr.p3)).ManhattanLength();

			//ht=(ht+cdr3)
			ht = (ht + cdr3);
			//--------------------------------------
			//dr3=0.
			//call wekt(p1,r1,dr3,d)      ! векторное произведение двух векторов вх-n,p0; вых-dr3 и d.
			//a=(0,0)
			//do i6=1,3
			//	a=a+et(i6)*dr3(i6)
			//end do
			dr3 = VecType(p1).CrossProduct(r1);			
			d = dr3.ManhattanLength();
			CVecType tmpDR3 = vector_cast<T, ComplexType>(dr3);
			a = (et*tmpDR3).ManhattanLength();

			//tes=0
			//dr3=0.
			//do I6=1,3
			//		tes=tes+p1(i6)*r1(i6)
			//end do
			tes = (p1*r1).ManhattanLength();

			//dr3=r1*tes
			//dr3=p1-dr3
			//do i6=1,3
			//	a=a+ht(i6)*dr3(i6)
			//end do
			dr3 = r1*tes;
			dr3 = p1 - dr3;
			tmpDR3 = vector_cast<T, ComplexType>(dr3);
			a = (ht*tmpDR3).ManhattanLength();
		}
		//f(i2)=a
		//end do 
		f[i] = a;
	}


	//---------------------------------------------------- 
	//if (abs(a1).lt.ad.and.abs(a2).lt.ad) then
	if (abs(a1) < ad && abs(a2) < ad)
	{
		//bi(1)=1./6.
		//bi(2)=1./6.
		//bi(3)=1./2.		
		bi[0] = T(1)/6;
		bi[1] = T(1)/6;
		bi[2] = T(1)/6;
	}
	else
	{
		//if (abs(a1).ge.ad.and.abs(a2).ge.ad.and.abs(a1-a2).lt.ad) then
		if (abs(a1) >= ad && abs(a2) >= ad && abs(a1-a2) < ad)
		{
			//bi(1)=-1./(a1*a1)-(1./a1-ji)/a1+(1/(a1*a1)+(1./a1-ji)*(1./a1-ji))*(cdexp(ji*a1)-1)/(ji*a1)
			//bi(1)=-bi(1)/2.
			//bi(2)=bi(1)
			//bi(3)=(cdexp(ji*a1)*(1.-ji*a1)-1.)/(a1*a1)
			bi[0] = -T(1)/(a1*a1)-(T(1)/a1-ji)/a1+(T(1)/(a1*a1)+(T(1)/a1-ji)*(T(1)/a1-ji))*(exp(ji*a1)-T(1))/(ji*a1);
			bi[0] = -bi[0]/ComplexType(2);
			bi[1] = bi[0];
			bi[2]=(exp(ji*a1)*(T(1)-ji*a1)-T(1))/(a1*a1);
		}
		//else if (abs(a1).lt.ad.and.abs(a2).ge.ad)
		else if (abs(a1) < ad && abs(a2)>=ad)
		{
			//bi(1)=-((cdexp(ji*a2)-1)/(ji*a2)-1-a2*ji/2)/(a2*a2)
			//bi(2)=-(1-(cdexp(ji*a2)-1)/(ji*a2)+(ji*cdexp(ji*a2)*a2-cdexp(ji*a2)+1)/(ji*a2))/(a2*a2)
			//bi(3)=((cdexp(ji*a2)-1)/(ji*a2)-1)/(ji*a2)
			bi[0] =-((exp(ji*a2)-T(1))/(ji*a2)-T(1)-a2*ji/ComplexType(2))/(a2*a2);
			bi[1] =-(T(1)-(exp(ji*a2)-T(1))/(ji*a2)+(ji*exp(ji*a2)*a2-exp(ji*a2)+T(1))/(ji*a2))/(a2*a2);
			bi[2] = ((exp(ji*a2)-T(1))/(ji*a2)-T(1))/(ji*a2);
		}
		//else if (abs(a2).lt.ad.and.abs(a1).ge.ad)
		else if (abs(a2) < ad && abs(a1)>=ad)
		{
			//bi(1)=-(1-(cdexp(ji*a1)-1)/(ji*a1)+(ji*cdexp(ji*a1)*a1-cdexp(ji*a1)+1)/(ji*a1))/(a1*a1)
			//bi(2)=-((cdexp(ji*a1)-1)/(ji*a1)-1-a1*ji/2)/(a1*a1)
			//bi(3)=((cdexp(ji*a1)-1)/(ji*a1)-1)/(ji*a1)
			bi[0]=-(T(1)-(exp(ji*a1)-T(1))/(ji*a1)+(ji*exp(ji*a1)*a1-exp(ji*a1)+T(1))/(ji*a1))/(a1*a1);
			bi[1]=-((exp(ji*a1)-T(1))/(ji*a1)-T(1)-a1*ji/T(2))/(a1*a1);
			bi[2]=((exp(ji*a1)-T(1))/(ji*a1)-T(1))/(ji*a1);
		}
		else
		{
			//bi(1)=-((cdexp(ji*a2)-1)/(ji*a2)-(cdexp(ji*a1)-1)/(ji*a1)-(a2-a1)*(ji*cdexp(ji*a1)*a1-cdexp(ji*a1)+1)/(ji*a1*a1))/((a2-a1)*(a2-a1))
			//bi(2)=-((cdexp(ji*a1)-1)/(ji*a1)-(cdexp(ji*a2)-1)/(ji*a2)-(a1-a2)*(ji*cdexp(ji*a2)*a2-cdexp(ji*a2)+1)/(ji*a2*a2))/((a2-a1)*(a2-a1))
			//bi(3)=((cdexp(ji*a2)-1)/(ji*a2)-(cdexp(ji*a1)-1)/(ji*a1))/(ji*(a2-a1))
			bi[0]=-((exp(ji*a2)-T(1))/(ji*a2)-(exp(ji*a1)-T(1))/(ji*a1)-(a2-a1)*(ji*exp(ji*a1)*a1-exp(ji*a1)+T(1))/(ji*a1*a1))/((a2-a1)*(a2-a1));
			bi[1]=-((exp(ji*a1)-T(1))/(ji*a1)-(exp(ji*a2)-T(1))/(ji*a2)-(a1-a2)*(ji*exp(ji*a2)*a2-exp(ji*a2)+T(1))/(ji*a2*a2))/((a2-a1)*(a2-a1));
			bi[2]= ((exp(ji*a2)-T(1))/(ji*a2)-(exp(ji*a1)-T(1))/(ji*a1))/(ji*(a2-a1));
		}
	}
	//------------------------------------------------
	//tes = 0.0
	//do i6=1,3
	//	tes=tes+(r1(i6)+rr(i6))*z(3,i6)
	//end do
	tes = ((r1+rr)*VecType(tr.p3)).ManhattanLength();

	ep = exp(-ji*k0*tes)*dd*((f[0]-f[2])*bi[1]+(f[1]-f[2])*bi[0]+f[2]*bi[2]);

	//return
	//end
	return;
}