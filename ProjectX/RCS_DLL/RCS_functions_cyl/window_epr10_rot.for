    !Подпрограмма расчета рассеянного цилиндром поля
    !FFIGn - азимут в градусах
    !PSn1 - угол места в градусах
    !kon0 - частота в ГГц
    !ei0 - относительная диэлектрическая проницаемость цилиндра
    !mi0 - относительная магнитная проницаемость цилиндра
    !tol0 - признак материала (0- металл, -1 - диэлектрик)
    !acr - радиус цилиндра в м
    !dleen - длина в м
    !epr,field - выходные данные (ЭПР в м2, поле )
      subroutine sux (FFIGn,PSn1,kon0,ei0,mi0,tol0,acr,dleen,
     +epr,field)
      use imsl
	real epr(3)
	complex field(3)
	real ukon,nkon(3),dbbi,bbi,bbi4

      INTEGER NAL,ARI
      INTEGER N,I6,I7,I8,J4,OO,OL,NN,I2,I3
      INTEGER I,ROG,J,L,K,MM,NI,NJ,M1,M2,II,KK,LL,NK
      INTEGER nid,nin,kni,njd,njn,knj,ci1,cj
      INTEGER TMRI0,istat,ip,cn,cm,kk1
      real tet0,tau, acr,dleen
      REAL*8 PI
      REAL DT,TI,TAUI,DFI,FIJ,TDL,kof,kon,kok,dk,kon1,kon0
      REAL CP,SP,FFI,CF,SF,SFG,CFG
      REAL ZZ(3),H(3),IS(3),X(3),X1(3),LOO(3),ZJK(3)
      REAL TMIN,TMAX,FI0,L01,GNAT,FFIG,PS,GG,GG1
      REAL R0(3),R1(3),rr(3),YJK(3),KJK(3)
      real zjk4(2,2,3),zjkk(3,3),n04(2,2,3),n0(3,3)
      real rt(3),rt1(3),rf(3),rf1(3)
      real dtau,taui4,fij4,ksi,dksi,ksi4,lyambda
      integer m41,m42
      REAL A00(3,3),B00(3,3),B11(3,3),AA0(3)
      REAL S0,B0,C0,D0,T11,T22
      integer nglob
      real PSn1,FFIGn,ffigk,dffig
      REAL P0(3),P1(3),PK0(3),PK1(3)
      REAL EKR(4),DKO(2)
      REAL LZ,DZ,RZ
      REAL CTT,STT,CFF,SFF,SKK,CKK
      REAL TIRI0,DTRI0
      REAL I10,FN0,ADF0,DHR0,DLR0,tol,e1r,e1m,m1r,m1m,tol0
      
      real ad1(3),ad0(3),af(3),pp0(3),spp,ad0n(3),tt0
      integer zpp0
      REAL E00(3,3),E11(3,3),emm(3,3)
      Real temp


      COMMON /S/ R0,R1
      common /s1/ RR                                            
      common /ngl/ nglob
      character rmr*8
      complex ji,d1,d2,cv,ad,bd,c1,c2,ei1,mi1,ei0,mi0
      complex ev
	complex*16 ep
      real en

      character (12) ob [allocatable] (:), ok [allocatable] (:)  
      real       ee [allocatable] (:,:), ff [allocatable] (:,:)
      integer nallm [allocatable] (:)
      REAL   DD0   [ALLOCATABLE] (:,:)



      complex ev0[allocatable](:,:)
      real en0[allocatable](:,:)
      REAL IH0	 [ALLOCATABLE] (:),	rldp  [ALLOCATABLE] (:)
      complex cak[ALLOCATABLE] (:),	akf   [ALLOCATABLE] (:)
      REAL MAS   [ALLOCATABLE] (:,:),   MZ    [ALLOCATABLE] (:,:)
      REAL COL   [ALLOCATABLE] (:,:),   UOL   [ALLOCATABLE] (:,:)
      REAL VOL   [ALLOCATABLE] (:,:),   LOES  [ALLOCATABLE] (:,:)
      REAL POL   [ALLOCATABLE] (:,:),   AA    [ALLOCATABLE] (:,:)
      REAL TF    [ALLOCATABLE] (:,:),   DD    [ALLOCATABLE] (:,:)
      REAL IL    [ALLOCATABLE] (:,:)  
      REAL MPP   [ALLOCATABLE] (:,:),   ALFA  [ALLOCATABLE] (:)
      REAL ALF   [ALLOCATABLE] (:),     LOS   [ALLOCATABLE] (:)
      REAL T1    [ALLOCATABLE] (:),     T2    [ALLOCATABLE] (:)
      REAL DL    [ALLOCATABLE] (:) 
      REAL A0    [ALLOCATABLE] (:,:,:), B     [ALLOCATABLE] (:,:,:)
      REAL B1    [ALLOCATABLE] (:,:,:), PP    [ALLOCATABLE] (:,:,:)
      REAL H1    [ALLOCATABLE] (:),     H2    [ALLOCATABLE] (:)
      INTEGER  NPP   [ALLOCATABLE] (:), ZPP   [ALLOCATABLE] (:,:)
      real ceps1 [allocatable] (:), cmi1 [allocatable] (:)
      real ceps2 [allocatable] (:), cmi2 [allocatable] (:)     
      real ctol [allocatable] (:)
      complex ep0[allocatable](:,:,:), ep1[allocatable](:,:,:,:)
      complex epk0[allocatable](:,:,:), epk1[allocatable](:,:,:,:)


      ji=(0,1)
	tet0=1.570796
      tau=1.5
  760 FORMAT(' Bл ўли«Ё §  Ја ­Ёжл ¤Ё Ї®§®­ ')
      pi=3.14159265
	kon1=30/kon0/100
      kon=2.*pi/kon1
      kok=kon  
      dk=1
      kk=1 
	kk1=1
      NID=1
      KNI=1
	! Количество разбиений на поверхности каждого элемента
      nin=400
      njn=401
      NJD=1
      KNJ=1
      psn = psn1+90
      psk=psn
      dps=1 

       ffigk=FFIGn 
       dffig=1 

       ggn=0
       ggk=0
       dgg=1

      knal=  1
      nal= 3

      allocate (ob(knal),ok(knal),ee(3,knal),ff(3,knal),nallm(knal))
 

      ALLOCATE (COL(NAL,3),UOL(NAL,3),VOL(NAL,3),LOES(NAL,3),POL(NAL,3),
     +AA(3,NAL),TF(3,NAL),DD(3,NAL),IL(NAL,3),ZPP(10,NAL),MPP(10,NAL),
     +ALFA(NAL),ALF(NAL),LOS(NAL),T1(NAL),T2(NAL),DL(NAL), 
     +NPP(NAL),A0(3,3,NAL),B(3,3,NAL),B1(3,3,NAL),PP(10,NAL,3))
      ALLOCATE (ceps1(nal),cmi1(nal),ctol(nal))
      allocate (ceps2(nal),cmi2(nal),ep1(kk,kni,knj,nal))
      allocate (EP0(kk,kni,knj),ev0(kk,3),en0(kk,3))
      allocate (EPk0(kk,kni,knj),epk1(kk,kni,knj,nal))
      
      allocate (DD0(3,nal))
! Подпрограмма вычисления параметров элементов 
      call fileindat (knal,nal,nallm,ee,ff,AA,TF,DD,ceps1,ceps2,
     +cmi1,cmi2,ctol,npp,pp,zpp,kon1,acr,dleen)
    ! Изменение размеров элементов в случае диэлектрического цилиндра
      if (tol0.ne.0.0) then
      PP( 1,  1,3)=  -dleen/2
      PP( 2,  1,3)=   dleen/2
      AA(1,  2)=   acr
      AA(2,  2)=   acr
      AA(1,  3)=   acr
      AA(2,  3)=   acr
	end if

      nom=1
      do n1=1,knal
      nall = nallm(n1)
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
      do n=1,3
      af(n)=ff(n,n1)
      end do
      do K=nom,nomen
      
      CFF=COS(TF(1,K)*PI/180)
      SFF=SIN(TF(1,K)*PI/180)
      CTT=COS(TF(2,K)*PI/180)
      STT=SIN(TF(2,K)*PI/180)
      CKK=COS(TF(3,K)*PI/180)
      SKK=SIN(TF(3,K)*PI/180)

      emm(1,1)=CTT*CKK
      emm(1,2)=-CTT*SKK
      emm(1,3)=STT
      emm(2,1)=STT*SFF*CKK+CFF*SKK
      emm(2,2)=CFF*CKK-SFF*STT*SKK
      emm(2,3)=-CTT*SFF
      emm(3,1)=-STT*CFF*CKK+SFF*SKK 
      emm(3,2)=SFF*CKK+CFF*STT*SKK
      emm(3,3)=CTT*CFF
      do n=1,3
      ad0(n)=dd(n,k)
      end do
      CALL GMPRD (e11,emm,b11,3,3,3)
      CALL GMPRD (e11,ad0,ad1,3,3,1)
      ad1=ad1+af
      call linrg (3,b11,3,b00,3)
      do n=1,3
      dd(n,k)=ad1(n)
      do nm=1,3
      b1(n,nm,k)=b11(n,nm)
      b(n,nm,k)=b00(n,nm)
      end do
	end do


      IF (NPP(k).GT.0) THEN 
      DO I7=1,NPP(k)
      do n=1,3
      pp0(n)=PP(I7,k,N)
      end do
      zpp0=ZPP(I7,k)

      CALL GMPRD (e11,pp0,ad0,3,3,1)
      
      ad1=ad0+af
	tt0=(ad0(1)*ad1(1)+ad0(2)*ad1(2)+ad0(3)*ad1(3))/
     1(ad0(1)**2+ad0(2)**2+ad0(3)**2)
	ad0=ad0*tt0
      if (tt0.lt.0) ZPP0=-1*ZPP0

       spp=0.
       DO N=1,3 
        spp=spp+ad0(N)**2
        END DO

       DO N=1,3 
       PP(I7,k,N)=ad0(n)
       end do
       zpp(i7,k)=zpp0
       MPP(i7,k)=spp



      END DO
      END IF




      end do
      nom=nom+nall
      end do
  


      DO N=1,3
       DO K=1,3
        DO I=1,NAL
          A0(N,K,I)=0.
        END DO
       END DO
      END DO
      DO 27 K=1,NAL
      A0(1,1,K)=AA(1,K)
      A0(2,2,K)=AA(2,K)
      A0(3,3,K)=AA(3,K)
   27 CONTINUE


      II=1
! Основной цикл по азимуту и углу места
      DO FFIG=FFIGn,FFIGk,dffig
      DO PS1=PSn,PSk,dps
        psr = ps1 
      PS=psr*PI/180. 
! Определение направлений облучения и приема и поляризаций
      CP=COS(PS)
      SP=SIN(PS)
      FFI=FFIG*PI/180.
      CF=COS(FFI)
      SF=SIN(FFI)
      R0(1)=CP
      R0(2)=-SP*SF
      R0(3)=-SP*CF
      rr(1)=-cp
      Rr(2)=SP*SF
      Rr(3)=SP*CF
      R1(1)=-CP
      DO GG=GGn,GGk,dgg
      GG1=GG*PI/180.
      SFG=SIN(FFI+GG1)
      CFG=COS(FFI+GG1)
      R1(2)=SP*SFG
      R1(3)=SP*CFG
      P1(1)=0.
      P1(2)=-CFG
      P1(3)=SFG
      PK1(1)=-SP
      PK1(2)=-CP*SFG
      PK1(3)=-CP*CFG
       DO II=1,1,1
      P0(1)=0
      P0(2)=CF
      P0(3)=-SF
      PK0(1)=SP
      PK0(2)=CP*SF
      PK0(3)=CP*CF
      do cj=1,knj
      do ci1=1,kni
      do i9=1,kk1
       ep0(i9,ci1,cj)=0
       epk0(i9,ci1,cj)=0
 
	do i=1,nal
       epk1(i9,ci1,cj,i)=0
	 ep1(i9,ci1,cj,i)=0
	end do
      end do
      end do
      end do
      CALL BPEMJ(B,A0,DD,VOL,UOL,POL,LOES,COL,ALFA,ALF,LOS,T1,T2,NAL)
      TMIN=T1(1)
      TMAX=T2(1)
      IF (NAL.GT.1) THEN 
      DO I=2,NAL
      IF (TMIN.GT.T1(I)) THEN
         TMIN=T1(I)
        END IF
      END DO
      DO I=2,NAL
      IF (TMAX.LT.T2(I)) THEN
         TMAX=T2(I)
       END IF
      END DO
      END IF 
      
  321 FORMAT(' TMIN=',F6.2,'TMAX=',F6.2)
      KK=1
      LL=1
      J4=1
      ep=0.
      

C Основной цикл по элементам гладкой поверхности (l=1 - боковая поверхность, l=2,3 - торцы цилиндра)
      DO L=1,NAL

         ei1=ei0
         mi1=mi0
         tol=tol0
         if (ei1.eq.1.0.and.mi1.eq.1.0) then
          tol=0.0
         end if
 
         ci1=0
         do while (ci1.lt.kni)
         ni=nin+nid*ci1
         ci1=ci1+1
         cj=0
         do while (cj.lt.knj)
         nj=njn+njd*cj
	 cj=cj+1
         dksi=pi/float(ni)
	if (L==1) then
      ukon=0
	dbbi=(PP( 2,  1,3)-PP( 1,  1,3))/ni
	end if

         do i=1,ni
         
	 ksi=-pi/2.0+dksi*(i-1)
	if (L==1) then
	bbi=PP( 2,  1,3)-dbbi*(i-1)
	end if
! Определение углов расчета для торцов
         taui=sin(ksi)
         if (taui .ge. 1.0) then 
          taui=1.0
         end if
	 if (taui.le.(alf(l)).or.L==1) then
          DO I6=1,3 
           ZZ(I6)=LOES(L,I6)
           END DO
          IF (TAUI.GE.-1.0.AND.TAUI.LE.-ALF(L)) THEN                                         
          FI0=-PI*0.5
	  ELSE
	   afio=TAUI/SQRT(1.-TAUI**2)/ALFA(L)
	bfio=TAUI/SQRT(1.-TAUI**2)/ALFA(L)
	IF (ABS(bfio).gt.1) bfio=bfio/abs(bfio)
           FI0=ASIN(bfio)      
          END IF
! Определение углов расчета для цилиндра
	if (L==1) then
		FI0=-PI*0.5
      end if

         DFI=(PI-2.*FI0)/NJ
         m41=0
         m42=0
         DO J=1,NJ
	  FIJ=FI0+DFI*(J-1)
           iin1=0
	     iin2=0
! Цикл по квадрату гладкой поверхности
            do i4=1,2
             do j4=1,2
               do i6=1,3
                zjk4(i4,j4,i6)=0
               end do
       
               ksi4=ksi+dksi*(i4-1)
	if (L==1) then
	bbi4=bbi-dbbi*(i4-1)
	end if

               taui4=sin(ksi4)
               fij4=fij+dfi*(j4-1)
            
          DO I6=1,3
            X(I6)=-TAUI4*ZZ(I6)+SQRT(1.-TAUI4**2)*(UOL(L,I6)*COS(FIJ4)+
     +            VOL(L,I6)*SIN(FIJ4))
            END DO
          K=3
          M2=0
          M1=0
c  Определение координат и нормалей точек на торцах         
           DO I6=1,3
             DO I7=1,3
               A00(I6,I7)=A0(I6,I7,L)
               B00(I6,I7)=B1(I6,I7,L)
               B11(I6,I7)=B(I6,I7,L)
               END DO
             END DO 
          CALL GMPRD (A00,X,X1,3,3,1)
          CALL GMPRD (B00,X1,LOO,3,3,1)
          if (abs(taui4).ne.1.0) then
                    DO I6=1,3
            rt(i6)=-zz(i6)-(UOL(L,I6)*COS(FIJ4)+VOL(L,I6)*SIN(FIJ4))*
     +            taui4/SQRT(1.-TAUI4**2)
            rf(i6)=SQRT(1.-TAUI4**2)*(-UOL(L,I6)*SIN(FIJ4)+
     +            VOL(L,I6)*COS(FIJ4))
            END DO
          CALL GMPRD (A00,rt,rt1,3,3,1)
          CALL GMPRD (B00,rt1,rt,3,3,1)
          CALL GMPRD (A00,rf,rf1,3,3,1)
          CALL GMPRD (B00,rf1,rf,3,3,1)
          call wekt(rt,rf,rt1,d)
          rt1=rt1/d

            DO I6=1,3 
              n04(i4,j4,i6)=-rt1(i6)
              END DO
          else
             n04(i4,j4,1)=loes(l,1)
             n04(i4,j4,2)=loes(l,2)
             n04(i4,j4,3)=loes(l,3)
          end if
             
            DO I6=1,3 
              ZJK(I6)=LOO(I6)+DD(I6,L)
              zjk4(i4,j4,i6)=zjk(i6)
              END DO
! Определение координат и нормалей на боковой поверхности
	if (L==1) then
	       n04(i4,j4,3)=0
             n04(i4,j4,1)=-cos(fij4)
             n04(i4,j4,2)=-sin(fij4)
             nkon(3)=0
             nkon(1)=-cos(fij4)
             nkon(2)=-sin(fij4)
		
         zjk4(i4,j4,3)=bbi4
         zjk4(i4,j4,2)=-(AA(1,  1))*sin(fij4)
         zjk4(i4,j4,1)=-(AA(1,  1))*cos(fij4)
         zjk(3)=bbi4
         zjk(2)=-(AA(1,  1))*sin(fij4)
         zjk(1)=-(AA(1,  1))*cos(fij4)

	end if

! Определение видимости точек по отсечениям элементов
            IF (NPP(L).GT.0) THEN 
               DO I8=1,NPP(L)
                 TDL=0.
                 DO I6=1,3
                  TDL=TDL+ZJK(I6)*PP(I8,L,I6)
                  END DO 
               IF ((ZPP(I8,L)*TDL).GT.(ZPP(I8,L)*MPP(I8,L))) THEN
                  m2=1
                 GOTO 60
                END IF
               END DO
            END IF
           M2=0
           MM=3
! Определение видимости боковой поверхности по ее освещенной части
	if (L==1) then
           vidkon=dot_product(nkon,r0)
           if (vidkon.GE.0.0) then
           m2=1
	     goto 60
           end if
	END IF
! Определение затенения точек другими элементами гладкой поверхности
            DO K=1,NAL
	if (k==L) cycle
              DO I6=1,3
                X1(I6)=ZJK(I6)-DD(I6,K)
                END DO
           DO I6=1,3
             DO I7=1,3
               B11(I6,I7)=B(I6,I7,k)
               END DO
             END DO 
            CALL GMPRD(B11,X1,YJK,3,3,1)
             S0=0.
             B0=0.
             C0=0.
              DO I6=1,3
            S0=S0+(POL(K,I6)/AA(I6,K))**2
            B0=B0+YJK(I6)*POL(K,I6)/(AA(I6,K)**2)
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
                   GOTO 60
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
               GOTO 60
             END IF
   72     END DO
             IF (M2.EQ.1) THEN 
                GOTO 60
               END IF
              DO I6=1,3
                AA0(I6)=AA(I6,L)
                END DO
   60 if (m2.EQ.1.and.i4.eq.1.and.j4.eq.1) then
         iin1=iin1+1
	   iin2=iin2+1
	end if
      if (m2.EQ.1.and.i4.eq.2.and.j4.eq.2) then
        iin1=iin1+1
	  iin2=iin2+1
      end if
      if (m2.EQ.1.and.i4.eq.1.and.j4.eq.2) then
        iin1=iin1+1
      end if
      if (m2.EQ.1.and.i4.eq.2.and.j4.eq.1) then
        iin2=iin2+1
      end if
c------------konets tsyklov vidimosti
              end do
            end do
! Вычисление поля в первом треугольнике квадрата
      if (iin1.lt.2.and.taui.ne.-1.0) then
        do i6=1,3
           zjkk(1,i6)=zjk4(1,1,i6)
           zjkk(2,i6)=zjk4(1,2,i6)
           zjkk(3,i6)=zjk4(2,2,i6)
           n0(1,i6)=n04(1,1,i6)
           n0(2,i6)=n04(1,2,i6)
           n0(3,i6)=n04(2,2,i6)
	end do
      do i9=1,kk1
      kof=kon+dk*(i9-1)
      call aarist (zjkk,kof,n0,ep,dob,ei1,mi1,tol,P0,P1)
      ep0(i9,ci1,cj)=ep0(i9,ci1,cj)+ep
	epk1(i9,ci1,cj,l)=epk1(i9,ci1,cj,l)+ep
      ep1(i9,ci1,cj,l)=ep1(i9,ci1,cj,l)+ep
! Вычисление поля для вертикальной поляризации диэлектрического цилиндра	     
      if(tol.ne.0.0) then 
      call aarist (zjkk,kof,n0,ep,dob,ei1,mi1,tol,Pk0,Pk1)
	end if
      epk0(i9,ci1,cj)=epk0(i9,ci1,cj)+ep     
      epk1(i9,ci1,cj,l)=epk1(i9,ci1,cj,l)+ep

      end do
      dob1=dob1+dob
      end if
! вычисление поля для второго треугольника в квадрате
      if (iin2.lt.2.and.taui.ne.1.0) then
        do i6=1,3
           zjkk(1,i6)=zjk4(1,1,i6)
           zjkk(2,i6)=zjk4(2,1,i6)
           zjkk(3,i6)=zjk4(2,2,i6)
           n0(1,i6)=n04(1,1,i6)
           n0(2,i6)=n04(2,1,i6)
           n0(3,i6)=n04(2,2,i6)
	end do
      do i9=1,kk1
      kof=kon+dk*(i9-1)
      call aarist (zjkk,kof,n0,ep,dob,ei1,mi1,tol,P0,P1)
      ep0(i9,ci1,cj)=ep0(i9,ci1,cj)+ep
      ep1(i9,ci1,cj,l)=ep1(i9,ci1,cj,l)+ep
! Вычисление поля для вертикальной поляризации диэлектрического цилиндра	     

      if(tol.ne.0.0) then 
      call aarist (zjkk,kof,n0,ep,dob,ei1,mi1,tol,Pk0,Pk1)
	end if
      epk0(i9,ci1,cj)=epk0(i9,ci1,cj)+ep     
      epk1(i9,ci1,cj,l)=epk1(i9,ci1,cj,l)+ep

      end do
      dob1=dob1+dob  
      end if



c Конец вычисления гладкой части



   61    END DO

	 end if
          
	 end do



      END DO
         end do
         end do
! Блок вычисления кромок
        ev0=0
	ev0=0
	en0=0
	if (tol0.eq.0.0) then
        do i9=1,kk1
        kof=kon+dk*(i9-1)
        lyambda=2*pi/kof
      call PointBle(P0,P1,lyambda,ev,en,ei0,mi0,tol0,tet0,tau,acr,dleen)
         ev0(i9,1)=ev
         en0(i9,1)=en
      call PointBle(PK0,PK1,lyambda,ev,en,ei0,mi0,tol0,tet0,tau,acr,
     +dleen)
         ev0(i9,2)=ev
         en0(i9,2)=en
      call PointBle(PK0,P1,lyambda,ev,en,ei0,mi0,tol0,tet0,tau,acr,
     +dleen)
         ev0(i9,3)=ev
         en0(i9,3)=en

        end do

      end if

	

! Блок сложения полей от гладкой и кромочной частей

         do ci1=1,kni
         do cj=1,knj
	 do i9=1,kk1
       	 kof=kon+dk*(i9-1)
	 ssert01=cabs(ep0(i9,ci1,cj)+ev0(i9,1))**2*kof*kof/4/pi
	 ssert02=cabs(epK0(i9,ci1,cj)+ev0(i9,2))**2*kof*kof/4/pi
	 ssert119=cabs(ev0(i9,3))*cabs(ev0(i9,3))*kof*kof/4/pi
	epr(1)=ssert01
	epr(2)=ssert02
	epr(3)=ssert119
	field(1)=(ep0(i9,ci1,cj)+ev0(i9,1))*kof
	field(2)=(epK0(i9,ci1,cj)+ev0(i9,2))*kof
	field(3)=(ev0(i9,3))*kof

          ep0(i9,ci1,cj)=0
         end do
         end do 
         end do
	 
	 END DO
	 END DO
	 END DO
	 END DO

	 return
	 END

C   /*Программа расчета поля на элементарном треугольнике*/
! z - координаты точек 
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
           real tol,z(3,3),r0(3),r1(3),rr(3),p0(3),p1(3),n0(3,3)
           complex*16 et(3),ht(3),a,bi(3),c,c0,c1,p1t(3),c2,ep
           complex*16 tes5
           real*8 da,db,dg,a1,a2,na,nb,nc,nd,h
           real*8 dr0(3),dr1(3),dr2(3)
	   real*8 p0t(3),rr1(3),r0t(3)
           real  r0p(3),dr3(3),d,n(3)
	   complex*16 f(3)
           real ad,k0,dd,tes
	COMMON /S/ R0,R1
	common /s1/ RR
	   ad=1e-3
           ji=(0,1)
           da=(z(2,2)-z(3,2))*(z(1,3)-z(3,3))-
     +      (z(1,2)-z(3,2))*(z(2,3)-z(3,3))
           db=(z(2,3)-z(3,3))*(z(1,1)-z(3,1))-
     +      (z(1,3)-z(3,3))*(z(2,1)-z(3,1))
           dg=(z(2,1)-z(3,1))*(z(1,2)-z(3,2))-
     +	   (z(2,2)-z(3,2))*(z(1,1)-z(3,1))
           dd=dsqrt(da*da+db*db+dg*dg)
           a1=0
           a2=0
           do i6=1,3
             dr0(i6)=rr(i6)+r1(i6)
             dr1(i6)=z(2,i6)-z(3,i6)
             dr2(i6)=z(1,i6)-z(3,i6)
             a1=a1+dr0(i6)*dr1(i6)
             a2=a2+dr0(i6)*dr2(i6)
           end do
           a1=-a1*k0
	   a2=-a2*k0
           do i2=1,3
           tes=0
           do I6=1,3
             tes=tes+rr(i6)*n0(i2,i6)
           end do
           if (tes.gt.0) then
           do i6=1,3
           n0(i2,i6)=-n0(i2,i6)
            end do
           end if
           end do
c Вычисление полного рассеянного поля в вершинах треугольника 
           do i2=1,3
            do i6=1,3
	     n(i6)=n0(i2,i6)
            end do
           tes=0
           do i6=1,3
            tes=tes+p0(i6)*n(i6)
           end do
           p0t=p0-n*tes
            tes=0
           do i6=1,3
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
           call wekt(n,rr,r0p,amod)
           r0t=rr-n*tes2
		 p1t=c1*p0t+2*ji*c/(ji*c0+1)*(r0t*(r0t(1)*p0(1)+r0t(2)*p0(2)
     #     +r0t(3)*p0(3))/(ji*c+tes2)+r0p*(r0p(1)*p0(1)+r0p(2)*p0(2)
     #     +r0p(3)*p0(3))/e1/m1/(ji*c+c2*c2/tes2))
	 else
!Расчет поля для диэлектрического треугольника
	     c=csqrt(m1/e1)*c0
           c0=c*tes2
           c1=(-c0-1)/(-c0+1)
           call wekt(n,rr,r0p,amod)
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

           call wekt(n,p0,dr3,d)
           tes=0
           do I6=1,3
             tes=tes+rr(i6)*z(3,i6)
           end do
           et=dr3
           cn=CMPLX(n)
	   call wektk(cn,pr1,cdr3)
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
           call wekt(p1,r1,dr3,d)
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
           else if (abs(a1).ge.ad.and.abs(a2).ge.ad.and.
     +	    abs(a1-a2).lt.ad) then
            bi(1)=-1./(a1*a1)-(1./a1-ji)/a1+(1/(a1*a1)+(1./a1-ji)*
     +        (1./a1-ji))*(cdexp(ji*a1)-1)/(ji*a1)
            bi(1)=-bi(1)/2.
            bi(2)=bi(1)
            bi(3)=(cdexp(ji*a1)*(1.-ji*a1)-1.)/(a1*a1)
           else if (abs(a1).lt.ad.and.abs(a2).ge.ad) then
            bi(1)=-((cdexp(ji*a2)-1)/(ji*a2)-1-a2*ji/2)/(a2*a2)
            bi(2)=-(1-(cdexp(ji*a2)-1)/(ji*a2)+(ji*cdexp(ji*a2)*a2-
     +        cdexp(ji*a2)+1)/(ji*a2))/(a2*a2)
            bi(3)=((cdexp(ji*a2)-1)/(ji*a2)-1)/(ji*a2)
           else if (abs(a2).lt.ad.and.abs(a1).ge.ad) then
            bi(1)=-(1-(cdexp(ji*a1)-1)/(ji*a1)+(ji*cdexp(ji*a1)*a1-
     +        cdexp(ji*a1)+1)/(ji*a1))/(a1*a1)
            bi(2)=-((cdexp(ji*a1)-1)/(ji*a1)-1-a1*ji/2)/(a1*a1)
            bi(3)=((cdexp(ji*a1)-1)/(ji*a1)-1)/(ji*a1)
           else
            bi(1)=-((cdexp(ji*a2)-1)/(ji*a2)-(cdexp(ji*a1)-1)/(ji*a1)-
     +       (a2-a1)*(ji*cdexp(ji*a1)*a1-cdexp(ji*a1)+1)/(ji*a1*a1))/
     +       ((a2-a1)*(a2-a1))
            bi(2)=-((cdexp(ji*a1)-1)/(ji*a1)-(cdexp(ji*a2)-1)/(ji*a2)-
     +       (a1-a2)*(ji*cdexp(ji*a2)*a2-cdexp(ji*a2)+1)/(ji*a2*a2))/
     +       ((a2-a1)*(a2-a1))
            bi(3)=((cdexp(ji*a2)-1)/(ji*a2)-(cdexp(ji*a1)-1)/(ji*a1))/
     +       (ji*(a2-a1))
	    end if
c------------------------------------------------
            tes=0.0
            do i6=1,3
             tes=tes+(r1(i6)+rr(i6))*z(3,i6)
	   end do
	    ep=cexp(-ji*k0*tes)*dd*((f(1)-f(3))*bi(2)+
     +        (f(2)-f(3))*bi(1)+f(3)*bi(3))
	    !ep = (9,10)
           return
           end
c____________________***___________________________

      subroutine aarist_comm (z,k0,n0,ep,dd,e1,m1,tol,P0,P1,S,s1)
           complex ji
           complex e1,m1,cn(3),cdr3(3),pr1(3)
           real tol,z(3,3),r0(3),r1(3),rr(3),p0(3),p1(3),n0(3,3)
           complex*16 et(3),ht(3),a,bi(3),c,c0,c1,p1t(3),c2,ep
           complex*16 tes5
           real*8 da,db,dg,a1,a2,na,nb,nc,nd,h
           real*8 dr0(3),dr1(3),dr2(3)
	   real*8 p0t(3),rr1(3),r0t(3)
           real  r0p(3),dr3(3),d,n(3)
	   complex*16 f(3)
           real ad,k0,dd,tes
	COMMON /S/ R0,R1
	common /s1/ RR
	   ad=1e-3
           ji=(0,1)
           da=(z(2,2)-z(3,2))*(z(1,3)-z(3,3))-
     +      (z(1,2)-z(3,2))*(z(2,3)-z(3,3))
           db=(z(2,3)-z(3,3))*(z(1,1)-z(3,1))-
     +      (z(1,3)-z(3,3))*(z(2,1)-z(3,1))
           dg=(z(2,1)-z(3,1))*(z(1,2)-z(3,2))-
     +	   (z(2,2)-z(3,2))*(z(1,1)-z(3,1))
           dd=dsqrt(da*da+db*db+dg*dg)
           a1=0
           a2=0
           do i6=1,3
             dr0(i6)=rr(i6)+r1(i6)
             dr1(i6)=z(2,i6)-z(3,i6)
             dr2(i6)=z(1,i6)-z(3,i6)
             a1=a1+dr0(i6)*dr1(i6)
             a2=a2+dr0(i6)*dr2(i6)
           end do
           a1=-a1*k0
	   a2=-a2*k0
           do i2=1,3
           tes=0
           do I6=1,3
             tes=tes+rr(i6)*n0(i2,i6)
           end do
           if (tes.gt.0) then
           do i6=1,3
           n0(i2,i6)=-n0(i2,i6)
            end do
           end if
           end do
c Вычисление полного рассеянного поля в вершинах треугольника 
           do i2=1,3
            do i6=1,3
	     n(i6)=n0(i2,i6)
            end do
           tes=0
           do i6=1,3
            tes=tes+p0(i6)*n(i6)
           end do
           p0t=p0-n*tes
            tes=0
           do i6=1,3
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
           call wekt(n,rr,r0p,amod)
           r0t=rr-n*tes2
		 p1t=c1*p0t+2*ji*c/(ji*c0+1)*(r0t*(r0t(1)*p0(1)+r0t(2)*p0(2)
     #     +r0t(3)*p0(3))/(ji*c+tes2)+r0p*(r0p(1)*p0(1)+r0p(2)*p0(2)
     #     +r0p(3)*p0(3))/e1/m1/(ji*c+c2*c2/tes2))
	 else
!Расчет поля для диэлектрического треугольника
	     c=csqrt(m1/e1)*c0
           c0=c*tes2
           c1=(-c0-1)/(-c0+1)
           call wekt(n,rr,r0p,amod)
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

           call wekt(n,p0,dr3,d)
           tes=0
           do I6=1,3
             tes=tes+rr(i6)*z(3,i6)
           end do
           et=dr3
           cn=CMPLX(n)
	   call wektk(cn,pr1,cdr3)
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
           call wekt(p1,r1,dr3,d)
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
           else if (abs(a1).ge.ad.and.abs(a2).ge.ad.and.
     +	    abs(a1-a2).lt.ad) then
            bi(1)=-1./(a1*a1)-(1./a1-ji)/a1+(1/(a1*a1)+(1./a1-ji)*
     +        (1./a1-ji))*(cdexp(ji*a1)-1)/(ji*a1)
            bi(1)=-bi(1)/2.
            bi(2)=bi(1)
            bi(3)=(cdexp(ji*a1)*(1.-ji*a1)-1.)/(a1*a1)
           else if (abs(a1).lt.ad.and.abs(a2).ge.ad) then
            bi(1)=-((cdexp(ji*a2)-1)/(ji*a2)-1-a2*ji/2)/(a2*a2)
            bi(2)=-(1-(cdexp(ji*a2)-1)/(ji*a2)+(ji*cdexp(ji*a2)*a2-
     +        cdexp(ji*a2)+1)/(ji*a2))/(a2*a2)
            bi(3)=((cdexp(ji*a2)-1)/(ji*a2)-1)/(ji*a2)
           else if (abs(a2).lt.ad.and.abs(a1).ge.ad) then
            bi(1)=-(1-(cdexp(ji*a1)-1)/(ji*a1)+(ji*cdexp(ji*a1)*a1-
     +        cdexp(ji*a1)+1)/(ji*a1))/(a1*a1)
            bi(2)=-((cdexp(ji*a1)-1)/(ji*a1)-1-a1*ji/2)/(a1*a1)
            bi(3)=((cdexp(ji*a1)-1)/(ji*a1)-1)/(ji*a1)
           else
            bi(1)=-((cdexp(ji*a2)-1)/(ji*a2)-(cdexp(ji*a1)-1)/(ji*a1)-
     +       (a2-a1)*(ji*cdexp(ji*a1)*a1-cdexp(ji*a1)+1)/(ji*a1*a1))/
     +       ((a2-a1)*(a2-a1))
            bi(2)=-((cdexp(ji*a1)-1)/(ji*a1)-(cdexp(ji*a2)-1)/(ji*a2)-
     +       (a1-a2)*(ji*cdexp(ji*a2)*a2-cdexp(ji*a2)+1)/(ji*a2*a2))/
     +       ((a2-a1)*(a2-a1))
            bi(3)=((cdexp(ji*a2)-1)/(ji*a2)-(cdexp(ji*a1)-1)/(ji*a1))/
     +       (ji*(a2-a1))
	    end if
c------------------------------------------------
            tes=0.0
            do i6=1,3
             tes=tes+(r1(i6)+rr(i6))*z(3,i6)
	   end do
	
	    !R0(1) = 25	
	    !R1 = 26
	    !RR = 27
	    ep=cexp(-ji*k0*tes)*dd*((f(1)-f(3))*bi(2)+
     +        (f(2)-f(3))*bi(1)+f(3)*bi(3))
	    !ep = (33,33)
           return
           end     


! Подпрограмма векторного произведения двух векторов
	   SUBROUTINE WEKT(A,B,C,D)
           REAL A(3),B(3),C(3),D                                       
           C(1)=A(2)*B(3)-A(3)*B(2)                                    
           C(2)=A(3)*B(1)-A(1)*B(3)                                    
           C(3)=A(1)*B(2)-A(2)*B(1)                                    
           D=SQRT(C(1)**2+C(2)**2+C(3)**2)                             
          RETURN                                                       
          END  
! Подпрограмма векторного произведения двух комплексных векторов
		                                                        
         SUBROUTINE WEKTK(RAK,L1K,RLK)
         COMPLEX RAK(3),L1K(3),RLK(3)
         RLK(1)=RAK(2)*L1K(3)-RAK(3)*L1K(2)
         RLK(2)=RAK(3)*L1K(1)-RAK(1)*L1K(3)
         RLK(3)=RAK(1)*L1K(2)-RAK(2)*L1K(1)
         RETURN
         END



! Программа расчета координатных данных для торцов цилиндра
         SUBROUTINE BPEMJ(B,A00,DD,VOL,UOL,POL,LOES,COL,ALFA,          
     +                    ALF,LOS,T1,T2,NL)                            
          INTEGER NL,JJ,I6,I7,I                                                
          REAL B(3,3,NL),A00(3,3,NL),DD(3,NL),T1(NL),T2(NL)                          
          REAL VOL(NL,3),UOL(NL,3),POL(NL,3),LOES(NL,3),COL(NL,3),     
     +    ALFA(NL),ALF(NL),LOS(NL),L01,LOO,VL(3),UL(3),TDL,ADL,BDL                 
          REAL BL(3,3),AL(3,3),A(3,3),EX,R0(3),     
     +    R1(3),R11(3),DI(3),PL(3),SL(3),CL(3),WL(3),LOE(3),L0(3),QL(3)
          COMMON /S/ R0,R1                                             
            DO JJ=1,NL                                                 
              DO I6=1,3                                                
                DO I7=1,3                                              
                  BL(I6,I7)=B(I6,I7,jj)                                
                  AL(I6,I7)=A00(I6,I7,jj)                              
                END DO                                                 
              END DO                                                   
             EX=1.0E-4                                                 
             CALL GMPRD(BL,R0,PL,3,3,1)                                
             CALL GMPRD(BL,R1,SL,3,3,1)                                
           CL=PL-SL                                                    
             CALL GMPRD(AL,CL,L0,3,3,1)                                
           LOO=SQRT(L0(1)**2+L0(2)**2+L0(3)**2)                        
           LOS(JJ)=LOO                                                 
           LOE=-L0/LOO 
           DO I6=1,3
             DO I7=1,3                                     
              A(I6,I7)=0.
             END DO
           END DO                                                       
            DO I=1,3                                                   
              A(I,I)=1.0/AL(I,I)                                       
            END DO                                                     
             CALL GMPRD(A,PL,QL,3,3,1)                                 
             CALL WEKT(LOE,QL,WL,L01)                                  
              IF (L01.LE.1E-4) THEN                                    
                R11=0.0                                                
                R11(1)=1.0                                             
                 DO I6=1,3                                             
                DI(I6)=AMIN1(ABS(LOE(I6)-R11(I6)),ABS(LOE(I6)+R11(I6))) 
                 END DO                                                
                 IF (DI(1).LE.EX.AND.DI(2).LE.EX.AND.DI(3).LE.EX) THEN 
                  R11=0.0                                              
                  R11(2)=1.0                                           
                 END IF                                                
              CALL WEKT(LOE,R11,UL,L01)                                
                UL=UL/L01                                              
                GOTO 91                                                
               END IF                                                  
             UL=-WL/L01                                                
   91         CALL WEKT(UL,LOE,VL,L01)                                 
                TDL=0.0                                                
                ADL=0.0                                                
                BDL=0.0                                                
              DO I6=1,3                                                
                TDL=TDL+QL(I6)*VL(i6)                                  
                ADL=ADL+QL(I6)*LOE(i6)                                 
                BDL=BDL+DD(I6,JJ)*(R0(I6)-R1(I6))                      
                UOL(JJ,I6)=UL(I6)                                      
                VOL(JJ,I6)=VL(I6)                                      
                COL(JJ,I6)=CL(I6)                                      
                POL(JJ,I6)=PL(I6)                                      
                LOES(JJ,I6)=LOE(I6)                                    
              END DO                                                   
             ALFA(JJ)=TDL/ADL                                          
             ALF(JJ)=ALFA(JJ)/SQRT(1.0+ALFA(JJ)**2)                    
             T1(JJ)=-LOO+BDL                                           
             T2(JJ)= LOO*ALF(JJ)+BDL                                   
            END DO                                                     
           RETURN                                                      
           END    
! произведение двух матриц		                                                      
      SUBROUTINE GMPRD(A,B,R,N,M,L)                                    
      DIMENSION A(*),B(*),R(*)
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
