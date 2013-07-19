!Подпрограмма расчета поля, рассеянного кромкой
!P0, P1 - поляризации передатчика и приемника
!lym - длина волны
!ev - поле, рассеянное кромкой
!epr1 - ЭПР кромки
!ei0,mi0 - относительные проницаемости материала покрытия кромки
!tol0 - радиус покрытия на кромках (здесь равен 0)
!tet00 - угол поворота кромки (здесь pi/2)
!tau0 - параметр раствора кромки (здесь 3/2)
!acr0 - радиус цилиндра
!dleen0 - длина цилиндра
	Subroutine PointBle(P0,P1,lym,ev,epr1,ei0,mi0,tol0,tet00,tau0,acr0,dleen0)		 
	use MSFLIB
        Real,dimension(3) :: N,Ne,nen,K,SmCentre,SmCentre1,SmCentre3,GlavNorP1,ObGlavNorP1,SechKr
	Real,dimension(3) :: BKrPoint1,BKrPoint2,bkrpoint11,bkrpoint22,Point1,Point2
	Real,dimension(3) :: BKrPoint1cos,BKrPoint2cos,SechKrcos,SechOb
	Real,dimension(3) :: Kon,Nach,NachKr
	Real,dimension(2) :: Fi,Fii
	Real d,d1,d2,par,l,kap1,kap2,Rp,Rm,k0,pi,d7,d8,ffi
	Real d4,d5,Proiz1,Proiz2,prom22,a,b,Znak,prom221
	complex,dimension(3) :: IntMl1,IntMl2,ObMLO,IO1,IO2,IO,ISK,IMLL,MLL,obMll
	complex:: jjj=(0,1)
        real d21,d22,k20,ad21,l21,s20
        complex ev,prom11,d3,temp,Epsil,Mu,promev
        complex ei0,mi0
        real tol0
	Integer Kol,i,Kolichestvo,NNN,index,param,n5,k5,number7,njfi0
	Real	Epr1,lyambda,tau,teto,r11,r2,tet00,tau0,acr0,dleen0
	Real,dimension(3) ::R0,Rn,VektRaz,R1,PPer,PP,akk,bkk
	real rr0(3),rr1(3),p0(3),p1(3),lym
	Real,dimension(3) ::RnKr,R0Kr,VektRazKr,PPerKr,PPKr,OrtVektRaz,OrtK
        Real temp111,ddc
	real pola, polb

      real ad(3),ad0(3),af(3),pp0(3),spp,ad0n(3),tt0,aKon(3),aNach(3),an(3),ane(3),ateto
      integer zpp0
      REAL E00(3,3),E11(3,3),emm(3,3)
      REAL A00(3,3),B00(3,3),B11(3,3),AA0(3)


	  COMMON /S/ Rr0,Rr1
	  common /indd/ index
	  common /qa/ Pola,Polb

      character (12) ob [allocatable] (:), ok [allocatable] (:)
      real       ee [allocatable] (:,:), ff [allocatable] (:,:)
      REAL   DD0   [ALLOCATABLE] (:,:), e1 [ALLOCATABLE] (:,:,:)
      integer kolm [allocatable] (:), nallm [allocatable] (:)

      REAL AA    [ALLOCATABLE] (:,:),   PPP    [ALLOCATABLE] (:,:,:)
      REAL TF    [ALLOCATABLE] (:,:),   DD    [ALLOCATABLE] (:,:)
      REAL MPP   [ALLOCATABLE] (:,:),   BB     [ALLOCATABLE] (:,:,:)
      INTEGER  NPP   [ALLOCATABLE] (:), ZPP   [ALLOCATABLE] (:,:)
      real ceps1 [allocatable] (:), cmi1 [allocatable] (:)
      real ceps2 [allocatable] (:), cmi2 [allocatable] (:)     
      real ctol [allocatable] (:), Epr11 [allocatable] (:)
	  complex evv [allocatable] (:) 

      real  aNachm [ALLOCATABLE] (:,:), aKonm [ALLOCATABLE] (:,:)
      real  atetom [ALLOCATABLE] (:), taum [ALLOCATABLE] (:)
      real  rpm [ALLOCATABLE] (:), rmm [ALLOCATABLE] (:)
      real  aNm [ALLOCATABLE] (:,:), aNem [ALLOCATABLE] (:,:)    
      complex Epsilm [ALLOCATABLE] (:), Mum [ALLOCATABLE] (:)

 
	  lyambda=lym; rn=rr1; r1=rr1; r0=rr0; 
	  pper=p0;
	  pp=p1;
      PI=3.141592653589; k0=2*Pi; promev=0
	  SmCentre(1)=0.0; SmCentre(2)=0.0; SmCentre(3)=0

      knal=1
      nAL=3
      kol=2
      allocate (ob(knal),ok(knal),ee(3,knal),ff(3,knal),kolm(knal),nallm(knal))
      ALLOCATE (AA(3,NAL),TF(3,NAL),DD(3,NAL),ZPP(10,NAL),MPP(10,NAL))
      ALLOCATE (NPP(NAL),BB(3,3,NAL),PPP(10,NAL,3))
      ALLOCATE (ceps1(nal),cmi1(nal),ctol(nal))
      allocate (ceps2(nal),cmi2(nal))

      allocate (DD0(3,nal),e1(3,3,knal))  

!Определение раазмеров гладкой части для алгоритма видимости
      call fileindat (knal,nal,nallm,ee,ff,AA,TF,DD,ceps1,ceps2,cmi1,cmi2,ctol,npp,ppp,zpp,lyambda,acr0,dleen0)

      PPp( 1,  1,3)=  -dleen0/2
      PPp( 2,  1,3)=   dleen0/2
      AA(1,  2)=   acr0
      AA(2,  2)=   acr0
      AA(1,  3)=   acr0
      AA(2,  3)=   acr0


      allocate  (aNachm(3,kol), aKonm(3,kol))
      allocate (atetom(kol), taum(kol), rpm(kol), rmm(kol))
      allocate (aNm(3,kol),aNem(3,kol))    
      allocate (Epsilm(kol),Mum(kol))


      nom=1
      do n1=1,knal
      nall=nallm(n1)
      nomen=nom+nall-1

      CFF=COS(ee(1,n1)*PI/180)
      SFF=SIN(ee(1,n1)*PI/180)
      CTT=COS(ee(2,n1)*PI/180)
      STT=SIN(ee(2,n1)*PI/180)
      CKK=COS(ee(3,n1)*PI/180)
      SKK=SIN(ee(3,n1)*PI/180)


      e11(1,1)=CTT*CKK
      e11(1,2)=-CTT*SKK
      e11(1,3)=STT
      e11(2,1)=STT*SFF*CKK+CFF*SKK
      e11(2,2)=CFF*CKK-SFF*STT*SKK
      e11(2,3)=-CTT*SFF
      e11(3,1)=-STT*CFF*CKK+SFF*SKK 
      e11(3,2)=SFF*CKK+CFF*STT*SKK
      e11(3,3)=CTT*CFF
      
      e1(1,1,n1)=e11(1,1)
      e1(1,2,n1)=e11(1,2)
      e1(1,3,n1)=e11(1,3)
      e1(2,1,n1)=e11(2,1)
      e1(2,2,n1)=e11(2,2)
      e1(2,3,n1)=e11(2,3)
      e1(3,1,n1)=e11(3,1)
      e1(3,2,n1)=e11(3,2)
      e1(3,3,n1)=e11(3,3)
      
      
      do nn=1,3
      af(nn)=ff(nn,n1)
      end do

      do Kk=nom,nomen


      
      CFF=COS(TF(1,Kk)*PI/180)
      SFF=SIN(TF(1,Kk)*PI/180)
      CTT=COS(TF(2,Kk)*PI/180)
      STT=SIN(TF(2,Kk)*PI/180)
      CKK=COS(TF(3,Kk)*PI/180)
      SKK=SIN(TF(3,Kk)*PI/180)

      emm(1,1)=CTT*CKK
      emm(1,2)=-CTT*SKK
      emm(1,3)=STT
      emm(2,1)=STT*SFF*CKK+CFF*SKK
      emm(2,2)=CFF*CKK-SFF*STT*SKK
      emm(2,3)=-CTT*SFF
      emm(3,1)=-STT*CFF*CKK+SFF*SKK 
      emm(3,2)=SFF*CKK+CFF*STT*SKK
      emm(3,3)=CTT*CFF
      do nn=1,3
      ad0(nn)=dd(nn,kk)
      end do

      CALL GMPRD (e11,emm,b11,3,3,3)
      CALL GMPRD (e11,ad0,ad,3,3,1)
      ad=ad+af
      call linrg (3,b11,3,b00,3)
      do nn=1,3
      dd(nn,kk)=ad(nn)
      do nm=1,3
      bb(nn,nm,kk)=b00(nn,nm)
      end do
	end do


      IF (NPP(kk).GT.0) THEN 
      DO I7=1,NPP(kk)
      do nn=1,3
      pp0(nn)=PPP(I7,kk,Nn)
      end do
      zpp0=ZPP(I7,kk)
 

     CALL GMPRD (e11,pp0,ad0,3,3,1)
      
      ad=ad0+af
	tt0=(ad0(1)*ad(1)+ad0(2)*ad(2)+ad0(3)*ad(3))/(ad0(1)**2+ad0(2)**2+ad0(3)**2)
	ad0=ad0*tt0
      if (tt0.lt.0) ZPP0=-1*ZPP0

       spp=0.
       DO Nn=1,3 
        spp=spp+ad0(Nn)**2
        END DO

       DO Nn=1,3 
       PPP(I7,kk,Nn)=ad0(nn)
       end do
       zpp(i7,kk)=zpp0
       MPP(i7,kk)=spp

      END DO
      END IF

      end do
      nom=nom+nall
      end do



!**********************************************************************
	!Вычисляется и нормируется вектор VektRaz (R0-Rn)
	VektRaz=R0-Rn; call ort(VektRaz,OrtVektRaz)



	allocate (evv(Kol*2))
	allocate (epr11(Kol*2))
	evv=0; epr11=0; ev=0; epr1=0
	prom22=0; prom221=0; prom11=0; promev=0
      nom=1
      do n1=1,knal
       nall = kolm(n1)
       nomen=nom+nall-1

      e11(1,1)=e1(1,1,n1)
      e11(1,2)=e1(1,2,n1)
      e11(1,3)=e1(1,3,n1)
      e11(2,1)=e1(2,1,n1)
      e11(2,2)=e1(2,2,n1)
      e11(2,3)=e1(2,3,n1)
      e11(3,1)=e1(3,1,n1)
      e11(3,2)=e1(3,2,n1)
      e11(3,3)=e1(3,3,n1)

      do nn=1,3
      af(nn)=ff(nn,n1)
      end do
	



!***************begin of circle part**********************
! начало вычисления кромок

do ik=1,2

  if (ik.eq.1) then
		teto=tet00; tau=tau0
		SmCentre1(1)=0; SmCentre1(2)=0; SmCentre1(3)=-dleen0/2
  else
		teto=-tet00; tau=tau0
				SmCentre1(1)=0; SmCentre1(2)=0; SmCentre1(3)=dleen0/2
  end if
        rm=0.5
        rp=tol0/lym
        if(ei0.eq.1.0.and.mi0.eq.1.0) then
        rp=0.3
        end if   

        N(1)=1.0; N(2)=0.0; N(3)=0.0
        Ne(1)=0.0; Ne(2)=0.0; Ne(3)=1.0

        Epsil=ei0
        Mu=mi0
        a=acr0+1e-3; b=acr0+1e-3
        NNN=0;

        Pola=a
        Polb=b 
		a=a/lym
		b=b/lym
	call WEKT1(Ne,N,K)
!************************************************************
!Блок проверки и вычисления свечения всей кромки 
if (r0(3)>=cosd(15.0).and.ik==1.or.r0(3)<cosd(165.0).and.ik==2) then
ak=0   ! Начальный предел интегрирования 
bk=2*pi  !Конечный предел интегрирования
njfi=20 ! Количество интервалов на которых применяется ф-ла Гаусса
SmCentre3=SmCentre1/lym

call all_edge(ak,bk,njfi,a,b,R0,PP,RN,VektRaz,Epsil,Mu,Rp,Rm,teto,tau,IO2,SmCentre3)
	prom11=dot_product(PPer,IO2)
	promev=prom11*120*pi*lyambda**2
	prom22=CABS(prom11)
	prom22=pi*(120*Pi*prom22*lyambda)**2	
	prom221=0
        ev=ev+promev             
        EPR1=EPR1+prom22 
else if (abs(r0(3))<cosd(3.0)) then 
!****************************************************************

	call KoordObeKrom (Ne,N,K,OrtVektRaz,SmCentre,VektRazKr)

!l=(-1)*b*VektRazKr(1)/(a*VektRazKr(2))
l=b*OrtVektRaz(1)/(a*OrtVektRaz(2))
	Fi(1)=ATanD(l)
	if (Fi(1).LT.0)	Then
					  Fi(1)=180+Fi(1); Fi(2)=180+Fi(1)
					Else
					  Fi(2)=180+Fi(1)
	end if

BKrPoint1(1)=b*sinD(Fi(1)); BKrPoint1(2)=a*cosD(Fi(1)); BKrPoint1(3)=0
BKrPoint2(1)=b*sinD(Fi(2)); BKrPoint2(2)=a*cosD(Fi(2)); BKrPoint2(3)=0

d1=a*b
d2=sqrt((a**2)*(sind(fi(1))**2)+(b**2)*(cosd(fi(1))**2))
d2=d2**3; kap1=d1/d2; kap1=abs(kap1)
d1=a*b
d2=sqrt((a**2)*(sind(fi(2))**2)+(b**2)*(cosd(fi(2))**2))
d2=d2**3; kap2=d1/d2; kap2=abs(kap2)
! блок определения видимости блестящих точек на кромках
do i=1,nnn
	Read(51,*)SechOb(1),SechOb(2),SechOb(3),par
	call KoordObeKrom (Ne,N,K,SechOb,SmCentre,SechKr)
	sechkr=sechkr/lyambda
	call Ort(BKrPoint1,BKrPoint1cos)
	call Ort(BKrPoint2,BKrPoint2cos)
	call Ort(SechKr,SechKrcos)
	Fii(1)=dot_product(SechKrcos,BKrPoint1cos)
	Fii(2)=dot_product(SechKrcos,BKrPoint2cos)
	BKrPoint1cos=BKrPoint1*Fii(1)
	BKrPoint2cos=BKrPoint2*Fii(2)
	d1=Sqrt(BKrPoint1cos(1)**2+BKrPoint1cos(2)**2+BKrPoint1cos(3)**2)
	d2=Sqrt(BKrPoint2cos(1)**2+BKrPoint2cos(2)**2+BKrPoint2cos(3)**2)
	d=Sqrt(SechKr(1)**2+SechKr(2)**2+SechKr(3)**2)
	if ((d.LT.d1).AND.(par.EQ.1).and.(fii(1).gt.0)) Then 
			BKrPoint1=555;
	end if
	if ((d.GT.d1).AND.(par.EQ.-1).and.(fii(1).gt.0)) Then 
			BKrPoint1=555;
	end if
	if ((d.LT.d2).AND.(par.EQ.1).and.(fii(2).gt.0)) Then 
			BKrPoint2=555;
	end if
	if ((d.GT.d2).AND.(par.EQ.-1).and.(fii(2).gt.0)) Then 
		  	BKrPoint2=555;
	end if
	if ((fii(1).LT.0).and.(par.EQ.-1)) Then
			BKrPoint1=555;
	end if
	if ((fii(2).LT.0).and.(par.EQ.-1)) Then
			BKrPoint2=555;
	end if
end do


Point1=BKrPoint1*lym+SmCentre1
Point2=BKrPoint2*lym+SmCentre1


if (Point1(1).LE.100) Then
   call vidim (Point1,nal,dd,bb,aa,ppp,npp,zpp,mpp,m2)
   if (m2.EQ.1) Then 
		  Point1=555
   end if
end if
if (Point2(1).LE.100) Then
   call vidim (Point2,nal,dd,bb,aa,ppp,npp,zpp,mpp,m2)
   if (m2.EQ.1) Then 
		  Point2=555
   end if
end if

!*************
!Проверка для Teto 
	call ort(VektRaz,OrtVektRaz)
	call ort(Ne,Ne)
	d7=dot_product(OrtVektRaz,Ne)
!Окончание проверки и пересчета Teto
!************
!  Расчет рассеяния от блестящих точек на кромках
call KoordObekrom(Ne,N,K,R0,SmCentre,R0Kr)
call KoordObekrom(Ne,N,K,Rn,SmCentre,RnKr)
call KoordObekrom(Ne,N,K,PP,SmCentre,PPKr)
fi(1)=fi(1)*pi/180; fi(2)=fi(2)*pi/180
prom22=0
prom221=0
if (point1(1).LT.50) Then
	call ml(fi(1),R0,PP,Rn,Epsil,Mu,Rp,Rm,teto,tau,IntMl1)
   ObMLO=IntMl1
	GlavNorP1(1)=-b*sin(Fi(1)); GlavNorP1(2)=-a*cos(Fi(1)); GlavNorP1(3)=0
	call RKoordKromObe (Ne,N,K,GlavNorP1,SmCentre,ObGlavNorP1)
	call RKoordKromObe (Ne,N,K,bkrpoint1,SmCentre1,bkrpoint11)	
	call Ort(ObGlavNorP1,ObGlavNorP1)
	Proiz1=dot_product(GlavNorP1,VektRaz)
	Proiz2=dot_product((bkrPoint1+SmCentre1/lym),VektRaz)
	if (Proiz1.LT.0) Then
						Znak=-1.0
					 Else
						Znak=1.0
	end if
	d3=cmplx(0,1)*k0*Proiz2+cmplx(0,1)*Znak*pi/4
	d3=cexp(d3); d5=ABS(Proiz1)
	d4=SQRT((2*Pi*a**2)/(k0*d5))
	temp=d3*d4
	IO1=d3*d4*ObMlO
	prom11=dot_product(PPer,IO1)
	promev=prom11*120*pi*lyambda**2
 	prom22=CABS(prom11)
	prom22=pi*(120*Pi*prom22*lyambda)**2	
 	prom221=prom22
	ev=ev+promev
        EPR1=EPR1+prom221
end if
prom22=0
if (point2(1).LT.50) Then
		call ml(fi(2),R0,PP,Rn,Epsil,Mu,Rp,Rm,teto,tau,IntMl2)
ObMLO=IntMl2
	GlavNorP1(1)=-b*sin(Fi(2)); GlavNorP1(2)=-a*cos(Fi(2)); GlavNorP1(3)=0
	call RKoordKromObe (Ne,N,K,GlavNorP1,SmCentre,ObGlavNorP1)
	call RKoordKromObe (Ne,N,K,bkrpoint2,SmCentre1,bkrpoint22)	
    call Ort(ObGlavNorP1,ObGlavNorP1)
		Proiz1=dot_product(GlavNorP1,VektRaz)
	Proiz2=dot_product((bkrPoint2+SmCentre1/lym),VektRaz)
	if (Proiz1.LT.0) Then
						Znak=-1
					 Else
						Znak=1
	end if
	d3=cexp(jjj*k0*Proiz2+jjj*Znak*pi/4); d5=ABS(Proiz1)
		d4=SQRT((2*Pi*a**2)/(k0*d5))
temp=d3*d4
    IO2=d3*d4*ObMlO
	prom11=dot_product(PPer,IO2)
	promev=prom11*120*pi*lyambda**2
        ev=ev+promev
	prom22=CABS(prom11)
	prom22=pi*(120*Pi*prom22*lyambda)**2	
        EPR1=EPR1+prom22
end if


 end if
	

end do



!***************end of cirle part*************************

        nom=nom+nall
        end do

!EPR1 = 23

return
end

!____________________***___________________________________________________________________

	Subroutine PointBle_comm(P0,P1,lym,ev,epr1,ei0,mi0,tol0,tet00,tau0,acr0,dleen0,S,indd,qa)		 
	use MSFLIB
        Real,dimension(3) :: N,Ne,nen,K,SmCentre,SmCentre1,SmCentre3,GlavNorP1,ObGlavNorP1,SechKr
	Real,dimension(3) :: BKrPoint1,BKrPoint2,bkrpoint11,bkrpoint22,Point1,Point2
	Real,dimension(3) :: BKrPoint1cos,BKrPoint2cos,SechKrcos,SechOb
	Real,dimension(3) :: Kon,Nach,NachKr
	Real,dimension(2) :: Fi,Fii
	Real d,d1,d2,par,l,kap1,kap2,Rp,Rm,k0,pi,d7,d8,ffi
	Real d4,d5,Proiz1,Proiz2,prom22,a,b,Znak,prom221
	complex,dimension(3) :: IntMl1,IntMl2,ObMLO,IO1,IO2,IO,ISK,IMLL,MLL,obMll
	complex:: jjj=(0,1)
        real d21,d22,k20,ad21,l21,s20
        complex ev,prom11,d3,temp,Epsil,Mu,promev
        complex ei0,mi0
        real tol0
	Integer Kol,i,Kolichestvo,NNN,index,param,n5,k5,number7,njfi0
	Real	Epr1,lyambda,tau,teto,r11,r2,tet00,tau0,acr0,dleen0
	Real,dimension(3) ::R0,Rn,VektRaz,R1,PPer,PP,akk,bkk
	real rr0(3),rr1(3),p0(3),p1(3),lym
	Real,dimension(3) ::RnKr,R0Kr,VektRazKr,PPerKr,PPKr,OrtVektRaz,OrtK
        Real temp111,ddc
	real pola, polb

      real ad(3),ad0(3),af(3),pp0(3),spp,ad0n(3),tt0,aKon(3),aNach(3),an(3),ane(3),ateto
      integer zpp0
      REAL E00(3,3),E11(3,3),emm(3,3)
      REAL A00(3,3),B00(3,3),B11(3,3),AA0(3)


	  COMMON /S/ Rr0,Rr1
	  common /indd/ index
	  common /qa/ Pola,Polb

      character (12) ob [allocatable] (:), ok [allocatable] (:)
      real       ee [allocatable] (:,:), ff [allocatable] (:,:)
      REAL   DD0   [ALLOCATABLE] (:,:), e1 [ALLOCATABLE] (:,:,:)
      integer kolm [allocatable] (:), nallm [allocatable] (:)

      REAL AA    [ALLOCATABLE] (:,:),   PPP    [ALLOCATABLE] (:,:,:)
      REAL TF    [ALLOCATABLE] (:,:),   DD    [ALLOCATABLE] (:,:)
      REAL MPP   [ALLOCATABLE] (:,:),   BB     [ALLOCATABLE] (:,:,:)
      INTEGER  NPP   [ALLOCATABLE] (:), ZPP   [ALLOCATABLE] (:,:)
      real ceps1 [allocatable] (:), cmi1 [allocatable] (:)
      real ceps2 [allocatable] (:), cmi2 [allocatable] (:)     
      real ctol [allocatable] (:), Epr11 [allocatable] (:)
	  complex evv [allocatable] (:) 

      real  aNachm [ALLOCATABLE] (:,:), aKonm [ALLOCATABLE] (:,:)
      real  atetom [ALLOCATABLE] (:), taum [ALLOCATABLE] (:)
      real  rpm [ALLOCATABLE] (:), rmm [ALLOCATABLE] (:)
      real  aNm [ALLOCATABLE] (:,:), aNem [ALLOCATABLE] (:,:)    
      complex Epsilm [ALLOCATABLE] (:), Mum [ALLOCATABLE] (:)

 
	  lyambda=lym; rn=rr1; r1=rr1; r0=rr0; 
	  pper=p0;
	  pp=p1;
      PI=3.141592653589; k0=2*Pi; promev=0
	  SmCentre(1)=0.0; SmCentre(2)=0.0; SmCentre(3)=0

      knal=1
      nAL=3
      kol=2
      allocate (ob(knal),ok(knal),ee(3,knal),ff(3,knal),kolm(knal),nallm(knal))
      ALLOCATE (AA(3,NAL),TF(3,NAL),DD(3,NAL),ZPP(10,NAL),MPP(10,NAL))
      ALLOCATE (NPP(NAL),BB(3,3,NAL),PPP(10,NAL,3))
      ALLOCATE (ceps1(nal),cmi1(nal),ctol(nal))
      allocate (ceps2(nal),cmi2(nal))

      allocate (DD0(3,nal),e1(3,3,knal))  

!Определение раазмеров гладкой части для алгоритма видимости
      call fileindat (knal,nal,nallm,ee,ff,AA,TF,DD,ceps1,ceps2,cmi1,cmi2,ctol,npp,ppp,zpp,lyambda,acr0,dleen0)

      PPp( 1,  1,3)=  -dleen0/2
      PPp( 2,  1,3)=   dleen0/2
      AA(1,  2)=   acr0
      AA(2,  2)=   acr0
      AA(1,  3)=   acr0
      AA(2,  3)=   acr0


      allocate  (aNachm(3,kol), aKonm(3,kol))
      allocate (atetom(kol), taum(kol), rpm(kol), rmm(kol))
      allocate (aNm(3,kol),aNem(3,kol))    
      allocate (Epsilm(kol),Mum(kol))


      nom=1
      do n1=1,knal
      nall=nallm(n1)
      nomen=nom+nall-1

      CFF=COS(ee(1,n1)*PI/180)
      SFF=SIN(ee(1,n1)*PI/180)
      CTT=COS(ee(2,n1)*PI/180)
      STT=SIN(ee(2,n1)*PI/180)
      CKK=COS(ee(3,n1)*PI/180)
      SKK=SIN(ee(3,n1)*PI/180)


      e11(1,1)=CTT*CKK
      e11(1,2)=-CTT*SKK
      e11(1,3)=STT
      e11(2,1)=STT*SFF*CKK+CFF*SKK
      e11(2,2)=CFF*CKK-SFF*STT*SKK
      e11(2,3)=-CTT*SFF
      e11(3,1)=-STT*CFF*CKK+SFF*SKK 
      e11(3,2)=SFF*CKK+CFF*STT*SKK
      e11(3,3)=CTT*CFF
      
      e1(1,1,n1)=e11(1,1)
      e1(1,2,n1)=e11(1,2)
      e1(1,3,n1)=e11(1,3)
      e1(2,1,n1)=e11(2,1)
      e1(2,2,n1)=e11(2,2)
      e1(2,3,n1)=e11(2,3)
      e1(3,1,n1)=e11(3,1)
      e1(3,2,n1)=e11(3,2)
      e1(3,3,n1)=e11(3,3)
      
      
      do nn=1,3
      af(nn)=ff(nn,n1)
      end do

      do Kk=nom,nomen


      
      CFF=COS(TF(1,Kk)*PI/180)
      SFF=SIN(TF(1,Kk)*PI/180)
      CTT=COS(TF(2,Kk)*PI/180)
      STT=SIN(TF(2,Kk)*PI/180)
      CKK=COS(TF(3,Kk)*PI/180)
      SKK=SIN(TF(3,Kk)*PI/180)

      emm(1,1)=CTT*CKK
      emm(1,2)=-CTT*SKK
      emm(1,3)=STT
      emm(2,1)=STT*SFF*CKK+CFF*SKK
      emm(2,2)=CFF*CKK-SFF*STT*SKK
      emm(2,3)=-CTT*SFF
      emm(3,1)=-STT*CFF*CKK+SFF*SKK 
      emm(3,2)=SFF*CKK+CFF*STT*SKK
      emm(3,3)=CTT*CFF
      do nn=1,3
      ad0(nn)=dd(nn,kk)
      end do

      CALL GMPRD (e11,emm,b11,3,3,3)
      CALL GMPRD (e11,ad0,ad,3,3,1)
      ad=ad+af
      call linrg (3,b11,3,b00,3)
      do nn=1,3
      dd(nn,kk)=ad(nn)
      do nm=1,3
      bb(nn,nm,kk)=b00(nn,nm)
      end do
	end do


      IF (NPP(kk).GT.0) THEN 
      DO I7=1,NPP(kk)
      do nn=1,3
      pp0(nn)=PPP(I7,kk,Nn)
      end do
      zpp0=ZPP(I7,kk)
 

     CALL GMPRD (e11,pp0,ad0,3,3,1)
      
      ad=ad0+af
	tt0=(ad0(1)*ad(1)+ad0(2)*ad(2)+ad0(3)*ad(3))/(ad0(1)**2+ad0(2)**2+ad0(3)**2)
	ad0=ad0*tt0
      if (tt0.lt.0) ZPP0=-1*ZPP0

       spp=0.
       DO Nn=1,3 
        spp=spp+ad0(Nn)**2
        END DO

       DO Nn=1,3 
       PPP(I7,kk,Nn)=ad0(nn)
       end do
       zpp(i7,kk)=zpp0
       MPP(i7,kk)=spp

      END DO
      END IF

      end do
      nom=nom+nall
      end do



!**********************************************************************
	!Вычисляется и нормируется вектор VektRaz (R0-Rn)
	VektRaz=R0-Rn; call ort(VektRaz,OrtVektRaz)



	allocate (evv(Kol*2))
	allocate (epr11(Kol*2))
	evv=0; epr11=0; ev=0; epr1=0
	prom22=0; prom221=0; prom11=0; promev=0
      nom=1
      do n1=1,knal
       nall = kolm(n1)
       nomen=nom+nall-1

      e11(1,1)=e1(1,1,n1)
      e11(1,2)=e1(1,2,n1)
      e11(1,3)=e1(1,3,n1)
      e11(2,1)=e1(2,1,n1)
      e11(2,2)=e1(2,2,n1)
      e11(2,3)=e1(2,3,n1)
      e11(3,1)=e1(3,1,n1)
      e11(3,2)=e1(3,2,n1)
      e11(3,3)=e1(3,3,n1)

      do nn=1,3
      af(nn)=ff(nn,n1)
      end do
	



!***************begin of circle part**********************
! начало вычисления кромок

do ik=1,2

  if (ik.eq.1) then
		teto=tet00; tau=tau0
		SmCentre1(1)=0; SmCentre1(2)=0; SmCentre1(3)=-dleen0/2
  else
		teto=-tet00; tau=tau0
				SmCentre1(1)=0; SmCentre1(2)=0; SmCentre1(3)=dleen0/2
  end if
        rm=0.5
        rp=tol0/lym
        if(ei0.eq.1.0.and.mi0.eq.1.0) then
        rp=0.3
        end if   

        N(1)=1.0; N(2)=0.0; N(3)=0.0
        Ne(1)=0.0; Ne(2)=0.0; Ne(3)=1.0

        Epsil=ei0
        Mu=mi0
        a=acr0+1e-3; b=acr0+1e-3
        NNN=0;

        Pola=a
        Polb=b 
		a=a/lym
		b=b/lym
	call WEKT1(Ne,N,K)
!************************************************************
!Блок проверки и вычисления свечения всей кромки 
if (r0(3)>=cosd(15.0).and.ik==1.or.r0(3)<cosd(165.0).and.ik==2) then
ak=0   ! Начальный предел интегрирования 
bk=2*pi  !Конечный предел интегрирования
njfi=20 ! Количество интервалов на которых применяется ф-ла Гаусса
SmCentre3=SmCentre1/lym

call all_edge(ak,bk,njfi,a,b,R0,PP,RN,VektRaz,Epsil,Mu,Rp,Rm,teto,tau,IO2,SmCentre3)
	prom11=dot_product(PPer,IO2)
	promev=prom11*120*pi*lyambda**2
	prom22=CABS(prom11)
	prom22=pi*(120*Pi*prom22*lyambda)**2	
	prom221=0
        ev=ev+promev             
        EPR1=EPR1+prom22 
else if (abs(r0(3))<cosd(3.0)) then 
!****************************************************************

	call KoordObeKrom (Ne,N,K,OrtVektRaz,SmCentre,VektRazKr)

!l=(-1)*b*VektRazKr(1)/(a*VektRazKr(2))
l=b*OrtVektRaz(1)/(a*OrtVektRaz(2))
	Fi(1)=ATanD(l)
	if (Fi(1).LT.0)	Then
					  Fi(1)=180+Fi(1); Fi(2)=180+Fi(1)
					Else
					  Fi(2)=180+Fi(1)
	end if

BKrPoint1(1)=b*sinD(Fi(1)); BKrPoint1(2)=a*cosD(Fi(1)); BKrPoint1(3)=0
BKrPoint2(1)=b*sinD(Fi(2)); BKrPoint2(2)=a*cosD(Fi(2)); BKrPoint2(3)=0

d1=a*b
d2=sqrt((a**2)*(sind(fi(1))**2)+(b**2)*(cosd(fi(1))**2))
d2=d2**3; kap1=d1/d2; kap1=abs(kap1)
d1=a*b
d2=sqrt((a**2)*(sind(fi(2))**2)+(b**2)*(cosd(fi(2))**2))
d2=d2**3; kap2=d1/d2; kap2=abs(kap2)
! блок определения видимости блестящих точек на кромках
do i=1,nnn
	Read(51,*)SechOb(1),SechOb(2),SechOb(3),par
	call KoordObeKrom (Ne,N,K,SechOb,SmCentre,SechKr)
	sechkr=sechkr/lyambda
	call Ort(BKrPoint1,BKrPoint1cos)
	call Ort(BKrPoint2,BKrPoint2cos)
	call Ort(SechKr,SechKrcos)
	Fii(1)=dot_product(SechKrcos,BKrPoint1cos)
	Fii(2)=dot_product(SechKrcos,BKrPoint2cos)
	BKrPoint1cos=BKrPoint1*Fii(1)
	BKrPoint2cos=BKrPoint2*Fii(2)
	d1=Sqrt(BKrPoint1cos(1)**2+BKrPoint1cos(2)**2+BKrPoint1cos(3)**2)
	d2=Sqrt(BKrPoint2cos(1)**2+BKrPoint2cos(2)**2+BKrPoint2cos(3)**2)
	d=Sqrt(SechKr(1)**2+SechKr(2)**2+SechKr(3)**2)
	if ((d.LT.d1).AND.(par.EQ.1).and.(fii(1).gt.0)) Then 
			BKrPoint1=555;
	end if
	if ((d.GT.d1).AND.(par.EQ.-1).and.(fii(1).gt.0)) Then 
			BKrPoint1=555;
	end if
	if ((d.LT.d2).AND.(par.EQ.1).and.(fii(2).gt.0)) Then 
			BKrPoint2=555;
	end if
	if ((d.GT.d2).AND.(par.EQ.-1).and.(fii(2).gt.0)) Then 
		  	BKrPoint2=555;
	end if
	if ((fii(1).LT.0).and.(par.EQ.-1)) Then
			BKrPoint1=555;
	end if
	if ((fii(2).LT.0).and.(par.EQ.-1)) Then
			BKrPoint2=555;
	end if
end do


Point1=BKrPoint1*lym+SmCentre1
Point2=BKrPoint2*lym+SmCentre1


if (Point1(1).LE.100) Then
   call vidim (Point1,nal,dd,bb,aa,ppp,npp,zpp,mpp,m2)
   if (m2.EQ.1) Then 
		  Point1=555
   end if
end if
if (Point2(1).LE.100) Then
   call vidim (Point2,nal,dd,bb,aa,ppp,npp,zpp,mpp,m2)
   if (m2.EQ.1) Then 
		  Point2=555
   end if
end if

!*************
!Проверка для Teto 
	call ort(VektRaz,OrtVektRaz)
	call ort(Ne,Ne)
	d7=dot_product(OrtVektRaz,Ne)
!Окончание проверки и пересчета Teto
!************
!  Расчет рассеяния от блестящих точек на кромках
call KoordObekrom(Ne,N,K,R0,SmCentre,R0Kr)
call KoordObekrom(Ne,N,K,Rn,SmCentre,RnKr)
call KoordObekrom(Ne,N,K,PP,SmCentre,PPKr)
fi(1)=fi(1)*pi/180; fi(2)=fi(2)*pi/180
prom22=0
prom221=0
if (point1(1).LT.50) Then
	call ml(fi(1),R0,PP,Rn,Epsil,Mu,Rp,Rm,teto,tau,IntMl1)
   ObMLO=IntMl1
	GlavNorP1(1)=-b*sin(Fi(1)); GlavNorP1(2)=-a*cos(Fi(1)); GlavNorP1(3)=0
	call RKoordKromObe (Ne,N,K,GlavNorP1,SmCentre,ObGlavNorP1)
	call RKoordKromObe (Ne,N,K,bkrpoint1,SmCentre1,bkrpoint11)	
	call Ort(ObGlavNorP1,ObGlavNorP1)
	Proiz1=dot_product(GlavNorP1,VektRaz)
	Proiz2=dot_product((bkrPoint1+SmCentre1/lym),VektRaz)
	if (Proiz1.LT.0) Then
						Znak=-1.0
					 Else
						Znak=1.0
	end if
	d3=cmplx(0,1)*k0*Proiz2+cmplx(0,1)*Znak*pi/4
	d3=cexp(d3); d5=ABS(Proiz1)
	d4=SQRT((2*Pi*a**2)/(k0*d5))
	temp=d3*d4
	IO1=d3*d4*ObMlO
	prom11=dot_product(PPer,IO1)
	promev=prom11*120*pi*lyambda**2
 	prom22=CABS(prom11)
	prom22=pi*(120*Pi*prom22*lyambda)**2	
 	prom221=prom22
	ev=ev+promev
        EPR1=EPR1+prom221
end if
prom22=0
if (point2(1).LT.50) Then
		call ml(fi(2),R0,PP,Rn,Epsil,Mu,Rp,Rm,teto,tau,IntMl2)
ObMLO=IntMl2
	GlavNorP1(1)=-b*sin(Fi(2)); GlavNorP1(2)=-a*cos(Fi(2)); GlavNorP1(3)=0
	call RKoordKromObe (Ne,N,K,GlavNorP1,SmCentre,ObGlavNorP1)
	call RKoordKromObe (Ne,N,K,bkrpoint2,SmCentre1,bkrpoint22)	
    call Ort(ObGlavNorP1,ObGlavNorP1)
		Proiz1=dot_product(GlavNorP1,VektRaz)
	Proiz2=dot_product((bkrPoint2+SmCentre1/lym),VektRaz)
	if (Proiz1.LT.0) Then
						Znak=-1
					 Else
						Znak=1
	end if
	d3=cexp(jjj*k0*Proiz2+jjj*Znak*pi/4); d5=ABS(Proiz1)
		d4=SQRT((2*Pi*a**2)/(k0*d5))
temp=d3*d4
    IO2=d3*d4*ObMlO
	prom11=dot_product(PPer,IO2)
	promev=prom11*120*pi*lyambda**2
        ev=ev+promev
	prom22=CABS(prom11)
	prom22=pi*(120*Pi*prom22*lyambda)**2	
        EPR1=EPR1+prom22
end if


 end if
	

end do



!***************end of cirle part*************************

        nom=nom+nall
        end do

!Rr0 = 31; Rr1 = 32
!index = 33
!Pola = 34; Polb = 35
!EPR1 = 55

return
end

! Пдпрограмма вычисления отражения от всей кромки (для фронтальных ракурсов)
subroutine all_edge(ak,bk,njfi,aaa,bbb,RN,PX,R,DR,EPS1,MI1,RP,RM,TETO,TAUO,ISK,Centre)
real ak,bk,MFIK(njfi,5),anj(5),r(3),RN(3),px(3),DR(3)
real LD,RP,RM,teto,tauo,d,Centre(3)
complex EPS1,MI1,MLK(3),ISK(3),dmk,SGF1,SGF2,SGF3
integer ii,i,njfi
real pi

pi=3.141592653589
CALL UZL(ak,bk,njfi,MFIK)
FI=0.0
ISK(1)=0
ISK(2)=0
ISK(3)=0
ANJ(1)=0.2369267
ANJ(2)=0.4786285
ANJ(3)=0.5688887
ANJ(4)=ANJ(2)
ANJ(5)=ANJ(1)
DO I=1,5
 SGF1=0
 SGF2=0
 SGF3=0
 DO ii=1,NJFI
  call ML(mfik(ii,i),RN,PX,R,EPS1,MI1,RP,RM,TETO,TAUO,MLK)
  Ld=sqrt((aaa**2)*(sin(mfik(ii,i))**2)+&
  (bbb**2)*(cos(mfik(ii,i))**2))
  dmk=CEXP(cmplx(0,1)*(2*PI*(DR(1)*(LD*sin(mfik(ii,i)+Centre(1)))+&
  DR(2)*(LD*cos(mfik(ii,i)+Centre(2)))+DR(3)*Centre(3) )))

  SGF1=SGF1+MLK(1)*DMK
  SGF2=SGF2+MLK(2)*DMK
  SGF3=SGF3+MLK(3)*DMK
  end do

  ISK(1)=ISK(1)+ANJ(I)*SGF1
      ISK(2)=ISK(2)+ANJ(I)*SGF2
      ISK(3)=ISK(3)+ANJ(I)*SGF3
  	end do
      ISK(1)=ISK(1)*(BK-AK)/(2*NJFI)
      ISK(2)=ISK(2)*(BK-AK)/(2*NJFI)
      ISK(3)=ISK(3)*(BK-AK)/(2*NJFI)
	  ISK=ISK*LD !"*" на LD так как интегрирование по дуге
return
end

! Пересчет в систему координат кромки
      SUBROUTINE KOORD(USF,TET,B,BM)                    
      REAL TAU1,TAU2,TAU3,Q(3),W1,W2,W3,USF,TET            
      REAL BM(3),B(3)
	  if (tet>0) then
      BM(1)=B(3)
      BM(2)=B(1)*sin(usf)+B(2)*cos(usf)                          
      BM(3)=B(1)*cos(usf)-B(2)*sin(usf)
	  else 
      BM(1)=-B(3)
      BM(2)=B(1)*sin(usf)+B(2)*cos(usf)                          
      BM(3)=-B(1)*cos(usf)+B(2)*sin(usf)
	  end if
	                                 
      RETURN
      END                                 
! Обратный пересчет из системы координат кромки

      SUBROUTINE KOORD1(USF,TET,B,BM)
      REAL TAU1,TAU2,TAU3,Q1,Q2,Q3,W1,W2,W3,USF,TET
      COMPLEX  BM(3),B(3)
	  if (tet>0) then
      BM(1)=B(3)*cos(usf)+B(2)*sin(usf)
      BM(2)=-B(3)*sin(usf)+B(2)*cos(usf)                          
      BM(3)=B(1)
	  else 
      BM(1)=-B(3)*cos(usf)+B(2)*sin(usf)
      BM(2)=B(3)*sin(usf)+B(2)*cos(usf)                          
      BM(3)=-B(1)
	  end if

      RETURN
      END                      


subroutine WEKT1(A,B,C)
!Вычисление вектора С перпендикулярного векторам A,B, т.е. векторное произведение      
      real,dimension(3) :: a,b,c
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)
      call Ort(C,C)
      return
      END 

! Вычисление узлов интегрирования для составной пятиточечной формулы Гаусса
	  subroutine UZL(A,B,NJ,mMI)
      integer I,J,nj,L
      real P1,P2,A,B
      real MMI(nj,5),x1(5)
      X1(1)=-0.9061797
      X1(2)=-0.5384691
      X1(3)=0
      X1(4)=-X1(2)
      X1(5)=-X1(1)
      DO l=1,NJ
       j=l-1
       P1=J*(B-A)/real(NJ,8)+A
       P2=(J+1)*(B-A)/real(NJ,8)+A
       DO I=1,5
        MMI(l,I)=(P1+P2)/2+X1(I)*(B-A)/(2*real(NJ,8))
       end do
      end do
      return
      END

! Обратный Пересчет центра координат кромки
  subroutine RKoordKromObe (A,B,C,D,F,E)
  real,dimension(3) :: a,b,c,f
  real, dimension(3) :: d,e


  E(1)=D(1)*A(1)+D(2)*B(1)+D(3)*C(1)
  E(2)=D(1)*A(2)+D(2)*B(2)+D(3)*C(2)
  E(3)=D(1)*A(3)+D(2)*B(3)+D(3)*C(3)
  return
  End
! Вычисление длины вектора
subroutine Ort(A,B)
  real,dimension(3) :: a,b
  real d
  d=sqrt(A(1)*A(1)+A(2)*A(2)+A(3)*A(3))
  B=A/d
  return
  End

! Прямой Пересчет центра координат кромки

  subroutine KoordObeKrom (A,B,C,D,F,E)
  real,dimension(3) :: a,b,c,d,f,e

  E(1)=(D(1))*A(1)+(D(2))*A(2)+(D(3))*A(3)
  E(2)=(D(1))*B(1)+(D(2))*B(2)+(D(3))*B(3)
  E(3)=(D(1))*C(1)+(D(2))*C(2)+(D(3))*C(3)
  return
  End


 
! Вычисления приведенных коєффициентов отраженного поля в базисных поляризациях
      subroutine UV(BR,BF,U,V,UR,UF,VR,VF,tau,kapo,mmax,r,kap)
      integer in,mmax,m
      complex cm1(41),cm2(41),am1(41),am2(41),bm1(41),bm2(41)
      common /ghjk/ cm1,cm2,am1,am2,bm1,bm2,in
      complex U,V,SU,SV,kap,z,js,ys,hs
      complex sur,SVR,SUF,SVF,UR,VR,UF,VF,HMP,HM1,IMP,IM1,HM
      real BR,BF,IM2,IM11,SIM,COM,kapo,r
	  real s,vm,tau
      U=0
      V=0
      UR=0
      UF=0
      VR=0
      VF=0
      IF (BR-R) 1,2,2
    2 DO m=0,MMAX
      IN=M+1
      VM=M/TAU
      SIM=SIN(VM*BF)
      COM=COS(VM*BF)
      S=VM
      Z=KAPO*BR
      CALL COMBES(S,Z,JS,YS,HS)
      HM=HS
      IM2=JS
      S=VM+1
      Z=KAPO*BR
      CALL COMBES(S,Z,JS,YS,HS)
      HMP=HS
      HM1=(VM/Z)*HM-HMP
      IM11=REAL(HM1,8)
      SU=(BM1(IN)*IM2+CM1(IN)*HM)*SIM
      SV=(BM2(IN)*IM2+CM2(IN)*HM)*COM
      SUR=KAPO*(BM1(IN)*IM11+CM1(IN)*HM1)*SIM
      SVR=KAPO*(BM2(IN)*IM11+CM2(IN)*HM1)*COM
      SUF=VM*(BM1(IN)*IM2+CM1(IN)*HM)*COM
      SVF=-VM*(BM2(IN)*IM2+CM2(IN)*HM)*SIM
      UR=UR+SUR
      VR=VR+SVR
      UF=UF+SUF
      VF=VF+SVF
      U=U+SU
      V=V+SV
     end do
      GO TO 50
    1 DO m=0,MMAX
      IN=M+1
      VM=M/TAU
      SIM=SIN(VM*BF)
      COM=COS(VM*BF)
      S=VM
      Z=KAP*BR
      CALL COMBES(S,Z,JS,YS,HS)
      IM=JS
      S=VM+1
      Z=KAP*BR
      CALL COMBES(S,Z,JS,YS,HS)
      IMP=JS
      IM1=(VM/Z)*IM-IMP
      SUR=KAP*AM1(IN)*IM1*SIM
      SVR=KAP*AM2(IN)*IM1*COM
      SUF=VM*AM1(IN)*IM*COM
      SVF=-VM*AM2(IN)*IM*SIM
      UR=UR+SUR
      VR=VR+SVR
      UF=UF+SUF
      VF=VF+SVF
      SU=AM1(in)+im*SIM
      SV=AM2(IN)*IM*COM
      U=U+SU
      V=V+SV
	end do
   50 continue
      return
      END
! Программа интегрирования по тору кромки
     subroutine ML(USF,RN,PX,R,EPS1,MI1,RP,RM,TETO,TAUO,MLO)
      real ww
      REAL MFI(5,5)
      integer I,ii,njf
      real PI
      real USF,TET,RN(3),PX(3),R(3),RNM(3),AF,BF,RP,MODGF
      real TETO  
      real TAU,TAUO
	  real RM,FI,PM(3),RPM(3)
      complex EPS1,MI1,IF1(3)
      complex sgf1,SGF2,SGF3,ww1,dmn
      complex MLO(3)
      REAL  ANJ(5)
      complex GF11(3),GF21(3)
      COMPLEX GF1(5,5),GF2(5,5),GF3(5,5)
      external IRSF
      PI=3.141592653589
	  TAU=TAUO
      TET=TETO
      AF=0
      BF=TAU*PI
      NJF=5
      CALL KOORD(USF,TET,RN,RNM)
      CALL KOORD(USF,TET,PX,PM)
      CALL KOORD(USF,TET,R,RPM)
      CALL IRSF(AF,BF,NJF,RNM(1),RNM(2),RNM(3),PM(1),PM(2),PM(3),&
	  EPS1,MI1,RP,RM,tau,RPM(1),RPM(2),RPM(3),GF1,GF2,GF3,mfi)
      ANJ(1)=0.2369267
      ANJ(2)=0.4786285
      ANJ(3)=0.5688887
      ANJ(4)=ANJ(2)
      ANJ(5)=ANJ(1)
      IF1=0
      IF2=0
      IF3=0
      DO I=1,5
      SGF1=0
      SGF2=0
      SGF3=0
      DO ii=1,NJF
      FI=MFI(ii,I)
      ww=-2*PI*RM*(RPM(1)*COS(FI)+RPM(2)*SIN(FI))
      ww1=cmplx(0.,ww)
      dmn=cexp(ww1)
      GF21(1)=GF1(ii,I)
      GF21(2)=GF2(ii,I)
      GF21(3)=GF3(ii,I)
      CALL KOORD1(USF,TET,GF21,GF11)
      SGF1=SGF1+GF11(1)*DMN
      SGF2=SGF2+GF11(2)*DMN
      SGF3=SGF3+GF11(3)*DMN
     end do
      IF1(1)=IF1(1)+ANJ(I)*SGF1
      IF1(2)=IF1(2)+ANJ(I)*SGF2
      IF1(3)=IF1(3)+ANJ(I)*SGF3
  	end do
      IF1(1)=IF1(1)*(BF-AF)/(2*NJF)
      IF1(2)=IF1(2)*(BF-AF)/(2*NJF)
      IF1(3)=IF1(3)*(BF-AF)/(2*NJF)
  	  MLO=RM*IF1
      return
      end
	 

      !Подпрограмма вычисления значений функции G
	  !tau-параметр раскрыва клина (0<tau<2)
	  !RO1M,RO2M,RO3M-волновой вектор падающей волны
	  !RN1M,RN2M,RN3M-волновой вектор отраженной волны
	  !M1,M2,M3- вектор поляризации падающей волны
	  !A,B- пределы интегрирования по нашлепке
	  !NJ-количество интервалов для вычисления интеграла по нашлёпке
      !EPS1,MI1- относительные электрическая и магнитная проницаемости
	  !среды (комплексные)
	  !MMAX-количество членов в разложении
	  !mfi-массив узлов для интегрирования
      subroutine irsf(A,B,NJ,RO1M,RO2M,RO3M,m1,m2,m3,EPS1,MI1,R,RM,&
      TAU,RN1M,RN2M,RN3M,GF1,GF2,GF3,mfi)
	  Use IMSL
      real RN1M,RN2M,RN3M,RO1M,RO2M,RO3M
      real mfi(nj,5)
      complex gf1(nj,5),gf2(nj,5),gf3(nj,5)
      real A,B,EEE
      integer NJ,II,M,MMAX,IN,i
      complex UR,UF,VR,VF
      complex JS(1),YS(1),HS
      complex KAP1,PM,QM,ME,KK
      complex IM,IMP,IM1,HM,HMP,HM1,ALF,EM,CM,DM,QM1,QM2,&
      U1B,U2B,U3B,V1B,V2B,V3B,Et1,Et2,Et3,H1,H2,H3
      complex E11,E12,E13,H11,H12,H13
      real  SMN1,tau
      real FIO,FI,RM,R,VM,KAPO,M1,M2,M3,KO,C1,C2,C3,S3,s
      real PI
      complex KAP,EPS1,MI1,KAP2,EPS2,MI2,ALFO
      real vv,vv1
      complex ww,ww1,z
      real A0,C0,S0
      real IM2,IM11,M12,ME0,SI0,CO0
      complex cm1(41),am1(41),cm2(41),am2(41),bm1(41),bm2(41)
      common /ghjk/ cm1,cm2,am1,am2,bm1,bm2,in
      data bm1 /41*0./
      data bm2 /41*0./
      data cm1 /41*0./
      data cm2 /41*0./
      data am1 /41*0./
      data am2 /41*0./
      PI=3.141592653589
      ME0=120*PI
      ME=ME0*CSQRT(MI1/EPS1)
      KO=2*PI
      C1=RO1M
      C2=RO2M
      C3=RO3M
      if (RO3M.gt.1.0) then
	  RO3M=1.0
	  endif
	  if (RO3M.lt.-1.0) then
	  RO3M=-1.0
	  endif   

      S3=SIN(ACOS(RO3M))
      KK=KO*CSQRT(EPS1*MI1)
      C0=-C1/S3		
      S0=-C2/S3			
      IF (ABS(C0).gt.1) C0=C0/abs(c0)
      A0=ASIN(C0)
      IF (S0.ge.0) THEN
      FIO=PI/2-A0
      ELSE
      FIO=3*PI/2+A0
      end if
      KAPO=KO*S3
      KAP=KO*CSQRT(EPS1*MI1-C3**2)
     M12=m2*C1-M1*C2
  KAP1=1/KAPO**2-1/KAP**2
      MMAX=TAU*(ABS(KAPO*RM)+4)
      DO m=0,mmax
      VM=M/TAU
      SI0=SIN(VM*FIO)
      CO0=COS(VM*FIO)
      IF (M.eq.0) THEN
      SMN1=2
      ELSE
      SMN1=1
      end if
      vv=-PI*vm*0.5
      ww=cmplx(0.,vv)
      ALFO=4/TAU/SMN1*cEXP(ww)
      IN=M+1
      S=VM
      Z=(KAP*R)
	  if (s==0.and.cabs(z)==0) then 
	  JS(1)=cmplx(1,0)
	  else
      if (s>0.and.cabs(z)==0) then
      JS(1)=cmplx(0,0)
	  else
	    CALL CBJS(S,Z,1,JS)
      end if
	  end if
	  IM=JS(1)
      S=VM+1
      Z=KAP*R
     	  if (s==0.and.cabs(z)==0) then 
	  JS(1)=cmplx(1,0)
	  else
       if (s>0.and.cabs(z)==0) then
      JS(1)=cmplx(0,0)
	  else
	    CALL CBJS(S,Z,1,JS)
      end if
	  end if
      IMP=JS(1)
      IM1=(VM/Z)*IM-IMP
      S=VM
      Z=KAPO*R
     	  if (s==0.and.cabs(z)==0) then 
	  JS(1)=cmplx(1,0)
	  else
       if (s>0.and.cabs(z)==0) then
      JS(1)=cmplx(0,0)
	  else
	    CALL CBJS(S,Z,1,JS)
      end if
	  end if
       CALL CBYS(S,Z,1,YS)
      HM=JS(1)+cmplx(0,1)*YS(1)
      IM2=JS(1)
      S=VM+1
      Z=KAPO*R
    	  if (s==0.and.cabs(z)==0) then 
	  JS(1)=cmplx(1,0)
	  else
       if (s>0.and.cabs(z)==0) then
      JS(1)=cmplx(0,0)
	  else
	    CALL CBJS(S,Z,1,JS)
      end if
	  end if
      CALL CBYS(S,Z,1,YS)
      HMP=JS(1)+cmplx(0,1)*YS(1)
      HM1=(VM/Z)*HM-HMP
      IM11=REAL(HM1,8)
      ALF=C3*VM/R*IM2*KAP1
      PM=-ME0/KAPO*IM11+KK*ME*IM1*IM2/(KAP*KO*IM)
      QM=1/(ME0*KAPO)*IM11-KK/ME*IM1*IM2/(KAP*KO*IM)
      EM=-C3*VM/R*HM*KAP1
      DM=KK/(KAP*KO*ME)*IM1*HM/IM-1/(KAPO*ME0)*HM1
      CM=-ME*KK/(KAP*KO)*IM1*HM/IM+ME0/KAPO*HM1
      QM1=ALFO*(ALF*M3*SI0+PM*M12/ME0*CO0)
      QM2=ALFO*(-ALF*M12/ME0*CO0+QM*M3*SI0)
      CM1(IN)=(QM2+EM*QM1/CM)/(EM/CM*EM+DM)
      CM2(IN)=(QM1-EM*QM2/DM)/(EM/DM*EM+CM)
      BM1(IN)=ALFO*M3*SI0
      BM2(IN)=ALFO*M12/ME0*CO0
      AM1(IN)=BM1(IN)*IM2/IM+CM1(IN)*HM/IM
      AM2(IN)=BM2(IN)*IM2/IM+CM2(IN)*HM/IM
      end do 
      CALL UZL(A,B,NJ,MFI)
	  DO ii=1,nj
      DO I=1,5
      FI=mfi(ii,i)
      CALL UV(RM,FI,U3B,V3B,UR,UF,VR,VF,tau,kapo,mmax,r,kap)
	  U1B=UR*COS(FI)-UF*SIN(FI)/RM
      U2B=UR*SIN(FI)+UF*COS(FI)/RM
      V1B=VR*COS(FI)-VF*SIN(FI)/RM
      V2B=VR*SIN(FI)+VF*COS(FI)/RM
      IF (RM.ge.R) THEN
      EPS2=1
      MI2=1
      KAP2=KAPO
      ELSE
      EPS2=EPS1
      MI2=MI1
      KAP2=KAP
      END if
      vv=KO*C3/(KAP2**2)
      ww=cmplx(0.,vv)
      vv1=120*pi*MI2/(KO*(EPS2*MI2-C3**2))
      ww1=cmplx(0.,vv1)
      Et1=ww*u1b+ww1*v2b
      Et2=ww*u2b-ww1*v1b
      Et3=U3B
      EEE=SQRT(cABS(Et1)**2+cABS(Et2)**2+cABS(Et3)**2)
      E11=-Et3*SIN(FI)
      E12=Et3*COS(FI)
      E13=-Et2*COS(FI)+Et1*SIN(FI)
      vv=-(1/(120*PI))*EPS2/(KO*(EPS2*MI2-C3**2))
      ww=cmplx(0.,vv)
      vv1=KO*C3/(KAP2**2)
      ww1=cmplx(0.,vv1)
      h1=ww*u2b+ww1*v1b
      ww=cmplx(0.,-vv)
      H2=ww*u1b+ww1*v2b
      H3=V3B
      H11=-H3*SIN(FI)
      H12=H3*COS(FI) 
	  H13=-H2*COS(FI)+H1*SIN(FI)
      GF1(ii,I)=(H11-(1/(120*PI))*(E12*RN3M-E13*RN2M))
      GF2(ii,I)=(H12-(1/(120*PI))*(E13*RN1M-E11*RN3M))
      GF3(ii,I)=(H13-(1/(120*PI))*(E11*RN2M-E12*RN1M))
     end do
	end do
      return
      end
!Вычисление гамма функции	 
      SUBROUTINE gmma(XX,GX,IER)
      INTEGER ier
	  REAL xx,GX,ERR,X,Y,GY
	  IF(XX-57.)6,6,4
    4 IER=2
      GX=1.E38
      RETURN
    6 X=XX
      ERR=1.0E-6
      IER=0
      GX=1.0
      IF(X-2.0)50,50,15
   10 IF(X-2.0)110,110,15
   15 X=X-1.0
      GX=GX*X
      GO TO 10
   50 IF(X-1.0)60,120,110
   60 IF(X-ERR)62,62,80
   62 Y=FLOAT(INT(X))-X
      IF(ABS(Y)-ERR)130,130,64
   64 IF(1.0-Y-ERR)130,130,70
   70 IF(X-1.0)80,80,110
   80 GX=GX/X
      X=X+1.0
      GO TO 70
  110 Y=X-1.0
      GY=1.0+Y*(-0.5771017+Y*(+0.9858540+Y*(-0.8764218+Y*(+0.8328212+ &
      Y*(-0.5684729+Y*(+0.2548205+Y*(-0.05149930)))))))
      GX=GX*GY
  120 RETURN
  130 IER=1
      RETURN
      END
!Вычисление функций бесселя нецелого индекса от комплексного аргумента  
      subroutine combes (S,Z,JS,YS,HS)
      real a,b
      external gmma
      integer K,L,ii,IER,JJ
      real q,u,v,w,v1,w1,u2,u1,s
      complex  Z,ZQ,ZQ1,P1,P2,P11,P12,S1,S2,S3,S4,ww
      complex s11,s12,cdmlp
      complex JQ,YQ,YQ1,JS,JQ1,YS,HS,JN,JQK,YQK,FAK
      real pi
      pi=3.141592653589
      Q=S-int(S)
      IF (cABS(Z)-8) 10,11,11
   10 call gmma(s+1,U,ier)
      a=1/u
      P1=cmplx(a)
      S1=p1
      DO K=1,80
       P1=-P1*(Z*z)/(4*K)/(S+K)
       S1=S1+P1
      end do
	  JS=cdmlp(Z/2,s)*S1
      IF (Q.ge.1E-2) goto 100
      S2=0
      DO L=1,S
       P1=1/real(L,8)
       S2=S2+P1
      end do
      FAK=1
      DO jj=1,s
       l=jj-1
       FAK=FAK/real(S-L,8)
      end do
      S3=FAK*S2
      DO K=1,80
       S1=0
       DO L=1,K
        P1=1/real(L,8)
        S1=S1+P1
       end do
       S2=0
       DO L=1,K+S
        P1=1/real(L,8)
        S2=S2+P1
       end do
       FAK=(-FAK/((S+K)*K))*cdmlp(Z*0.5,2.0)
       S3=S3+FAK*(S1+S2)
     end do
      FAK=1
      ii=int(s-1)
      DO L=1,ii
       FAK=FAK*(S-L)
      end do
      IF (S.eq.0) THEN
      S4=0
      ELSE
      S4=FAK
      end if
      ii=s-1
      DO K=1,ii
	   FAK=FAK/((S-K)*K)*cdmlp(Z/2,2.0)
       S4=S4+FAK
      end do
      YS=(2/PI)*JS*cLOG(1.781072418*Z/2)-  &
      (1/PI)*cdmlp(Z/2,s)*S3-(1/PI)*cdmlp(Z/2,-s)*S4
      GO TO 101
  100 S=-S
      call gmma(s+1,U,ier)
      P1=1/U
      S1=1/U
      DO K=1,80
       P1=-P1*cdmlp(Z,2.0)/(4E0*K)/(S+K)
       S1=S1+P1
      end do
      JN=cdmlp(Z/2,S)*S1
      YS=(1/SIN((-S)*PI))*(JS*COS((-S)*PI)-JN)
      GO TO 101
   11 call gmma(q+.5,U,ier)
      V=U
      call gmma(q-.5,w,ier)
      ZQ=Z-((2*Q+1)/4E0)*PI
      call gmma(q+1.5e0,u1,ier)
      V1=U1
      call gmma(q+.5,w1,ier)
      call gmma(q+2.5,u2,ier)
      ZQ1=Z-((2*(Q+1)+1)/4E0)*PI
      P1=(1,0)
      IF (abs(w)-1e30) 30,30,31
  31  P2=(0.,0.)
      goto 32
  30  P2=U1/(W*2*Z)
  32  P11=(1,0)
      P12=U2/(W1*2*Z)
      S1=P1
      S2=P2
      S11=P11
      S12=P12
      DO K=1,4
      P1=-P1*(Q+.5+2*(K-1))*(Q+.5+2*K-1)*&
      (Q+.5-2*K)*(Q+.5-2*K+1)/(4*(Z**2)*(2*K-1)*2*K)
      P2=-P2*(Q+1.5+2*(K-1))*(Q+1.5+2*K-1)*	&
      (Q-.5-2*K)*(Q-.5-2*K+1)/(4*(Z**2)*(2*K+1)*2*K)
      P11=-P11*(Q+1.5+2*(K-1))*(Q+1.5+2*K-1)*  &
      (Q+1.5-2*K)*(Q+1.5-2*K+1)/(4*(Z**2)*(2*K-1)*2*K)
      P12=-P12*(Q+2.5+2*(K-1))*(Q+2.5+2*K-1)* &
      (Q+.5-2*K)*(Q+.5-2*K+1)/(4*(Z**2)*(2*K+1)*2*K)
      S1=S1+P1
      S2=S2+P2
      S11=S11+P11
       S12=S12+P12
      end do
      JQ=CSQRT(2/(PI*Z))*(S1*CCOS(ZQ)-S2*CSIN(ZQ))
      JS=JQ
      YQ=CSQRT(2/(PI*Z))*(S1*CSIN(ZQ)+S2*CCOS(ZQ))
      YS=YQ
      IF (int(S).eq.0)  GO TO 101
      JQ1=CSQRT(2/(PI*Z))*(S11*CCOS(ZQ1)-S12*CSIN(ZQ1))
      JS=JQ1
      YQ1=CSQRT(2/(PI*Z))*(S11*CSIN(ZQ1)+S12*CCOS(ZQ1))
      YS=YQ1
      IF (int(S).eq.1) GO TO 101
      jj=int(s)
      DO K=2,jj
        JQK=(2*(Q+K-1)/Z)*JQ1-JQ
        JQ=JQ1
       JQ1=JQK
       JS=JQK
       YQK=(2*(Q+K-1)/Z)*YQ1-YQ
       YQ=YQ1
       YQ1=YQK
       YS=YQK
      end do
  101 a=real(ys)
      b=imag(ys)
      ww=cmplx(-b,a)
      HS=JS+ww
      return 
	end
   function cdmlp(z,s)
  complex cdmlp,z
  real s
  if  (real(z,8)==0.and.imag(z)==0) then
  cdmlp=z
  else
  cdmlp=z**s
  end if
  if  (real(z,8)==0.and.imag(z)==0 .and. s<=0) then
  write(*,*) 'ERROR! Dividing by zero!!'
  stop
  end if
  end function

!  Перемножение матриц
      SUBROUTINE GMPRD1(A,B,R,N,M,L)
      real :: A(*),B(*),R(*)
      INTEGER IR,IK,N,M,J,K,JI,IB,I                                   
      IR=0                                                             
      IK=-M                                                            
      DO 100 K=1,L                                                     
      IK=IK+M                                                          
      DO 100 J=1,N                                                     
      IR=IR+1                                                          
      JI=J-N                                                           
      IB=IK                                                            
      R(IR)=0                                                          
      DO 100 I=1,M                                                     
      JI=JI+N                                                          
      IB=IB+1                                                          
  100 R(IR)=R(IR)+A(JI)*B(IB)                                          
      RETURN                                                           
      END  
! Подпрогрмма видимости блестящих точек кромки и их затенения гладкими элементами поверхности
         subroutine vidim (zjk,nal,dd,b,aa,pp,npp,zpp,mpp,m2)
        integer nal,m1,m2
        integer npp(nal),zpp(10,nal)
        real zjk(3),x1(3),yjk(3),kjk(3),mpp(10,nal),tdl
        real dd(3,nal),b11(3,3),aa(3,nal),PP(10,NAL,3)
        REAL B(3,3,NAL),SL(3),PL(3)
        REAL S0,B0,C0,D0,T11,T22
        real r0(3),r1(3),rr(3)
		real rr0(3),rr1(3)
        COMMON /S/ Rr0,Rr1
		r1=rr1
		r0=rr0
         do i6=1,3
         rr(i6)=-r1(i6)
         end do
            M1=0
            M2=0 
            DO K=1,NAL
			  DO I6=1,3
             DO I7=1,3
               B11(I6,I7)=B(I6,I7,K)
               END DO
             END DO
             CALL GMPRD1(B11,rr,SL,3,3,1) 
             CALL GMPRD1(B11,r0,PL,3,3,1)
              DO I6=1,3
                X1(I6)=ZJK(I6)-DD(I6,K)
                END DO
            CALL GMPRD1(B11,X1,YJK,3,3,1)
             S0=0.
             B0=0.
             C0=0.
              DO I6=1,3
            S0=S0+(PL(I6)/AA(I6,K))**2
            B0=B0+YJK(I6)*PL(I6)/(AA(I6,K)**2)
            C0=C0+(YJK(I6)/AA(I6,K))**2
              END DO
            D0=B0**2-S0*(C0-1)
            IF (D0.LT.0.) THEN 
               M1=0
              ELSE 
               T11=(-B0+SQRT(D0))/S0
               T22=(-B0-SQRT(D0))/S0
                IF (ABS(C0-1).LT.1.E-4) THEN 
                  IF (ABS(T11).LT.ABS(T22)) THEN
                    T11=0. 
                   ELSE
                    T22=0.
                   END IF
                END IF
                IF (ABS(T11).LT.1.E-2.AND.ABS(T22).LT.1.E-2) THEN 
                 T11=0
                 T22=0
                END IF
                IF (T11.LT.0) THEN 
                    DO I6=1,3
                  KJK(I6)=ZJK(I6)+R0(I6)*T11
                    END DO                     
                  IF (NPP(K).GT.0) THEN 
                    DO I8=1,NPP(K)
                      TDL=0.
                       DO I6=1,3
                         TDL=TDL+KJK(I6)*PP(I8,K,I6)
                         END DO
                    IF ((ZPP(I8,K)*TDL).GT.(ZPP(I8,K)*MPP(I8,K))) THEN 
                        M1=0
                        GOTO 71
                      ELSE 
                        M1=1
                      END IF
                    END DO
                  ELSE 
                    M1=1
                  END IF
                END IF
                IF (M1.EQ.1) THEN 
                   M2=1
                   return
                END IF
   71           IF (T22.LT.0.) THEN 
                  DO I6=1,3 
                   KJK(I6)=ZJK(I6)+R0(I6)*T22
                   END DO
                  IF (NPP(K).GT.0) THEN
                    DO I8=1,NPP(K)
                      TDL=0.
                       DO I6=1,3
                         TDL=TDL+KJK(I6)*PP(I8,K,I6)
                         END DO
                 IF ((ZPP(I8,K)*TDL).GT.(ZPP(I8,K)*MPP(I8,K))) THEN 
                         M1=0
                         GOTO 72
                       ELSE 
                         M1=1
                       END IF
                     END DO
                   ELSE
                     M1=1
                   END IF 
                 END IF
             END IF
             IF (M1.EQ.1) THEN 
               M2=1
               return
            END IF
    72    CONTINUE        
             S0=0.
             B0=0.
             C0=0.
              DO I6=1,3
            S0=S0+(SL(I6)/AA(I6,K))**2
            B0=B0+YJK(I6)*SL(I6)/(AA(I6,K)**2)
            C0=C0+(YJK(I6)/AA(I6,K))**2
              END DO
            D0=B0**2-S0*(C0-1)
            IF (D0.LT.0.) THEN 
               M1=0
              ELSE 
               T11=(-B0+SQRT(D0))/S0
               T22=(-B0-SQRT(D0))/S0
                IF (ABS(C0-1).LT.1.E-4) THEN 
                  IF (ABS(T11).LT.ABS(T22)) THEN
                    T11=0. 
                   ELSE
                    T22=0.
                   END IF
                END IF
                IF (ABS(T11).LT.1.E-2.AND.ABS(T22).LT.1.E-2) THEN 
                 T11=0
                 T22=0
                END IF
                IF (T11.LT.0) THEN 
                    DO I6=1,3
                  KJK(I6)=ZJK(I6)+Rr(I6)*T11
                    END DO                     
                  IF (NPP(K).GT.0) THEN 
                    DO I8=1,NPP(K)
                      TDL=0.
                       DO I6=1,3
                         TDL=TDL+KJK(I6)*PP(I8,K,I6)
                         END DO
                    IF ((ZPP(I8,K)*TDL).GT.(ZPP(I8,K)*MPP(I8,K))) THEN 
                        M1=0
                        GOTO 711
                      ELSE 
                        M1=1
                      END IF
                    END DO
                  ELSE 
                    M1=1
                  END IF
                END IF
                IF (M1.EQ.1) THEN 
                   M2=1
                   return
                END IF
  711           IF (T22.LT.0.) THEN 
                  DO I6=1,3 
                   KJK(I6)=ZJK(I6)+R0(I6)*T22
                   END DO
                  IF (NPP(K).GT.0) THEN
                    DO I8=1,NPP(K)
                      TDL=0.
                       DO I6=1,3
                         TDL=TDL+KJK(I6)*PP(I8,K,I6)
                         END DO
                 IF ((ZPP(I8,K)*TDL).GT.(ZPP(I8,K)*MPP(I8,K))) THEN 
                         M1=0
                         GOTO 722 
                       ELSE 
                         M1=1
                       END IF
                     END DO
                   ELSE
                     M1=1
                   END IF 
                 END IF
             END IF
             IF (M1.EQ.1) THEN 
               M2=1
               return
            END IF
          

  722     END DO
       return
       end   


