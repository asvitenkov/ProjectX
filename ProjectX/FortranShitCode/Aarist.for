!   /*Программа расчета поля на элементарном треугольнике*/
!z  - координаты точек 
!k0 - волновое число
!n0 - нормали к вершинам треугольника
!ep - расчитанное поле
!dd - площадь треугольника
!e1,m1 - относительные проницаемости
!tol - признак (0-металл, -1 - диэлектрик)
!P0,P1 - поляризации передатчика и приемника
subroutine aarist (z,k0,n0,ep,dd,e1,m1,tol,P0,P1)
  complex ji
  complex e1,m1,cn(3),cdr3(3),pr1(3)
  complex*16 et(3),ht(3),a,bi(3),c,c0,c1,p1t(3),c2,ep
  complex*16 tes5
  complex*16 f(3)

  real tol,z(3,3),r0(3),r1(3),rr(3),p0(3),p1(3),n0(3,3)
  real*8 da,db,dg,a1,a2,na,nb,nc,nd,h
  real*8 dr0(3),dr1(3),dr2(3)
  real*8 p0t(3),rr1(3),r0t(3)
  real  r0p(3),dr3(3),d,n(3)
  real ad,k0,dd,tes
	
  COMMON /S/ R0,R1			!S=VM  (VM=M/TAU;   tau=1.5) / R0(1)=CP, R0(2)=-SP*SF, R0(3)=-SP*CF / R1(1)=-CP, R1(2)=SP*SFG, R1(3)=SP*CFG
	common /s1/ RR				!s1=P1 (P1 - поляризации ) /RR(1)=-cp, RR(2)=SP*SF, RR(3)=SP*CF
  
  ad=1e-3
  ji=(0,1)

  ! Площадь
  da=(z(2,2)-z(3,2))*(z(1,3)-z(3,3))-(z(1,2)-z(3,2))*(z(2,3)-z(3,3)) ! da=z(,)*z(,)-z(,)*z(,)
  db=(z(2,3)-z(3,3))*(z(1,1)-z(3,1))-(z(1,3)-z(3,3))*(z(2,1)-z(3,1)) ! db=z(,)*z(,)-z(,)*z(,) 
  dg=(z(2,1)-z(3,1))*(z(1,2)-z(3,2))-(z(2,2)-z(3,2))*(z(1,1)-z(3,1)) ! dg=z(,)*z(,)-z(,)*z(,)
  dd=dsqrt(da*da+db*db+dg*dg)                                        ! dd=0.5*dsqrt(da*da+db*db+dg*dg)
  
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

! Вычисление полного рассеянного поля в вершинах треугольника 
  do i2=1,3 ! НОМЕР ТОЧКИ ТРЕУГОЛЬНИКА
    do i6=1,3! КООРДИНАТА ТРЕУГОЛЬНИКА(X=1),(Y=2), (Z=3)  
      n(i6)=n0(i2,i6)
    end do
    tes=0

    do i6=1,3 ! НОМЕР ТОЧКИ ТРЕУГОЛЬНИКА
      tes=tes+p0(i6)*n(i6)
    end do
    p0t=p0-n*tes
    tes=0

    do i6=1,3 ! НОМЕР ТОЧКИ ТРЕУГОЛЬНИКА
      tes=tes+rr(i6)*n(i6)
    end do
    a=0
    if (abs(tes).gt.0.001) then
      tes2=tes
      c0=csqrt(1.-(1-tes2*tes2)/(e1*m1))
      c2=c0
      ! Расчет поля для металлического треугольника
      if (tol.ge.0.0) then 
        c1=k0*csqrt(e1*m1)*tol*c0
        c=csqrt(m1/e1)*c0*cdsin(c1)/cdcos(c1)
        c0=c*tes2
        c1=(ji*c0-1)/(ji*c0+1)
        call wekt(n,rr,r0p,amod)   ! векторное произведение двух векторов вх-n,rr; вых-r0p и amod
        r0t=rr-n*tes2
        p1t=c1*p0t+2*ji*c/(ji*c0+1)*(r0t*(r0t(1)*p0(1)+r0t(2)*p0(2)
        #     +r0t(3)*p0(3))/(ji*c+tes2)+r0p*(r0p(1)*p0(1)+r0p(2)*p0(2)
        #     +r0p(3)*p0(3))/e1/m1/(ji*c+c2*c2/tes2))
      else
        !Расчет поля для диэлектрического треугольника
        c=csqrt(m1/e1)*c0
        c0=c*tes2
        c1=(-c0-1)/(-c0+1)
        call wekt(n,rr,r0p,amod)    ! векторное произведение двух векторов вх-n,rr; вых-r0p и amod
        r0t=rr-n*tes2
        p1t=c1*p0t+2*c/(-c0+1)*(r0t*(r0t(1)*p0(1)+r0t(2)*p0(2)
        #     +r0t(3)*p0(3))/(c-tes2)+r0p*(r0p(1)*p0(1)+r0p(2)*p0(2)
        #     +r0p(3)*p0(3))/e1/m1/(c-c2*c2/tes2))
      end if
     
      tes5=0
      do i6=1,3
        tes5=tes5+rr(i6)*p1t(i6)
      end do

      pr1=p1t+n*tes5/tes
      rr1=rr-n*tes*2.0 

      call wekt(n,p0,dr3,d)       ! векторное произведение двух векторов вх-n,p0; вых-dr3 и d.
      tes=0

      do I6=1,3
        tes=tes+rr(i6)*z(3,i6)
      end do

      et=dr3
      cn=CMPLX(n)
      call wektk(cn,pr1,cdr3)         !векторного произведения двух комплексных векторов вх-cn,pr1; вых-cdr3.
  	  et=(et+cdr3)
      dr3=0.
      tes=0

      do I6=1,3
        tes=tes+rr(i6)*n(i6)
      end do
      dr3=p0*tes
      tes=0
      do I6=1,3
        tes=tes+p0(i6)*n(i6)
      end do
      ht=dr3-rr*tes

      ! Вычисление интеграла по трем точкам
      dr3=0.
      tes=0
      do I6=1,3
        tes=tes+rr1(i6)*n(i6)
      end do
      
      cdr3=pr1*tes
      tes5=0
      do I6=1,3
        tes5=tes5+pr1(i6)*n(i6)
      end do
      
      cdr3=cdr3-rr1*tes5
      tes=0
      do I6=1,3
        tes=tes+rr(i6)*z(3,i6)
      end do
      
      ht=(ht+cdr3)
      c--------------------------------------
      dr3=0.
      call wekt(p1,r1,dr3,d)      ! векторное произведение двух векторов вх-n,p0; вых-dr3 и d.
      a=(0,0)
      do i6=1,3
        a=a+et(i6)*dr3(i6)
      end do
      
      tes=0
      dr3=0.
      do I6=1,3
      tes=tes+p1(i6)*r1(i6)
      end do
      dr3=r1*tes
      dr3=p1-dr3
      do i6=1,3
        a=a+ht(i6)*dr3(i6)
      end do
    end if
  f(i2)=a
  end do 

c---------------------------------------------------- 
if (abs(a1).lt.ad.and.abs(a2).lt.ad) then
  bi(1)=1./6.
  bi(2)=1./6.
  bi(3)=1./2.
  else 
    if (abs(a1).ge.ad.and.abs(a2).ge.ad.and.abs(a1-a2).lt.ad) then
      bi(1)=-1./(a1*a1)-(1./a1-ji)/a1+(1/(a1*a1)+(1./a1-ji)*(1./a1-ji))*(cdexp(ji*a1)-1)/(ji*a1)
      bi(1)=-bi(1)/2.
      bi(2)=bi(1)
      bi(3)=(cdexp(ji*a1)*(1.-ji*a1)-1.)/(a1*a1)
    else 
      if (abs(a1).lt.ad.and.abs(a2).ge.ad) then
        bi(1)=-((cdexp(ji*a2)-1)/(ji*a2)-1-a2*ji/2)/(a2*a2)
        bi(2)=-(1-(cdexp(ji*a2)-1)/(ji*a2)+(ji*cdexp(ji*a2)*a2-cdexp(ji*a2)+1)/(ji*a2))/(a2*a2)
        bi(3)=((cdexp(ji*a2)-1)/(ji*a2)-1)/(ji*a2)
      else 
        if (abs(a2).lt.ad.and.abs(a1).ge.ad) then
          bi(1)=-(1-(cdexp(ji*a1)-1)/(ji*a1)+(ji*cdexp(ji*a1)*a1-cdexp(ji*a1)+1)/(ji*a1))/(a1*a1)
          bi(2)=-((cdexp(ji*a1)-1)/(ji*a1)-1-a1*ji/2)/(a1*a1)
          bi(3)=((cdexp(ji*a1)-1)/(ji*a1)-1)/(ji*a1)
        else
          bi(1)=-((cdexp(ji*a2)-1)/(ji*a2)-(cdexp(ji*a1)-1)/(ji*a1)-(a2-a1)*(ji*cdexp(ji*a1)*a1-cdexp(ji*a1)+1)/(ji*a1*a1))/((a2-a1)*(a2-a1))
          bi(2)=-((cdexp(ji*a1)-1)/(ji*a1)-(cdexp(ji*a2)-1)/(ji*a2)-(a1-a2)*(ji*cdexp(ji*a2)*a2-cdexp(ji*a2)+1)/(ji*a2*a2))/((a2-a1)*(a2-a1))
          bi(3)=((cdexp(ji*a2)-1)/(ji*a2)-(cdexp(ji*a1)-1)/(ji*a1))/(ji*(a2-a1))
        end if
!------------------------------------------------
      tes=0.0
      do i6=1,3
        tes=tes+(r1(i6)+rr(i6))*z(3,i6)
      end do

      ep=cexp(-ji*k0*tes)*dd*((f(1)-f(3))*bi(2)+(f(2)-f(3))*bi(1)+f(3)*bi(3))
      return
      end