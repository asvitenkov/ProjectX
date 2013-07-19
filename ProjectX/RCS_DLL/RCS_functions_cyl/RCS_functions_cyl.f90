!  RCS_functions_cyl.f90 
!
!  FUNCTIONS/SUBROUTINES exported from RCS_functions_cyl.dll:
!	RCS_functions_cyl      - subroutine 

!   Программа расчета поля на элементарном треугольнике
! z - координаты точек 
!k0 - волновое число
!n0 - нормали к вершинам треугольника
!ep - расчитанное поле
!dd - площадь треугольника
!e1,m1 - относительные проницаемости
!tol - признак (0-металл, -1 - диэлектрик)
!P0,P1 - поляризации передатчика и приемника
!
subroutine aarist_dll (i_z,i_k0,i_n0,i_dd,i_e1,i_m1,i_tol,i_P0,i_P1,i_ep)

! Expose subroutine RCS_functions_cyl to users of this DLL
  !
  !DEC$ ATTRIBUTES DLLEXPORT::aarist_dll

  ! Variables
  real z(3,3), k0, n0(3,3), dd, tol, P0(3), P1(3)
  complex*16 ep
  complex e1, m1

  real i_z(3,3), i_k0, i_n0(3,3), i_dd, i_tol, i_P0(3), i_P1(3)
  complex*16 i_ep
  complex i_e1, i_m1
   
  ! Body of aarist_dll
  z = i_z; k0 = i_k0; n0 = i_n0; dd = i_dd; tol = i_tol; P0 = i_P0; P1 = i_P1
  e1 = i_e1; m1 = i_m1

  call aarist(z,k0,n0,ep,dd,e1,m1,tol,P0,P1)
  i_ep = ep

end subroutine aarist_dll

!_____________________________________________________________________________________

subroutine aarist_dll_common (i_z,i_k0,i_n0,i_dd,i_e1,i_m1,i_tol,i_P0,i_P1,i_R0,i_R1,i_RR,i_ep)

! Expose subroutine RCS_functions_cyl to users of this DLL
  !
  !DEC$ ATTRIBUTES DLLEXPORT::aarist_dll_common

  ! Variables
  real z(3,3), k0, n0(3,3), dd, tol, P0(3), P1(3)
  real R0(3), R1(3), RR(3)
  complex*16 ep
  complex e1, m1

  real i_z(3,3), i_k0, i_n0(3,3), i_dd, i_tol, i_P0(3), i_P1(3)
  real i_R0(3), i_R1(3), i_RR(3)
  complex*16 i_ep
  complex i_e1, i_m1

  COMMON /S/ R0,R1
  common /s1/ RR
     
  ! Body of aarist_dll_common
  z = i_z; k0 = i_k0; n0 = i_n0; dd = i_dd; tol = i_tol; P0 = i_P0; P1 = i_P1
  e1 = i_e1; m1 = i_m1
  !r0 = i_r0; r1 = i_r1; rr = i_rr

  call aarist_comm(z,k0,n0,ep,dd,e1,m1,tol,P0,P1,S,s1)
  i_ep = ep
  i_R0 = R0
  i_R1 = R1
  i_RR = RR
  !i_ep = (10,10)

end subroutine aarist_dll_common

!_____________________________________________________________________________________

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

subroutine PointBle_dll(i_P0,i_P1,i_lym,i_ev,i_ei0,i_mi0,i_tol0,i_tet00,i_tau0,i_acr0,i_dleen0,i_epr1)

! Expose subroutine RCS_functions_cyl to users of this DLL
  !
  !DEC$ ATTRIBUTES DLLEXPORT::PointBle_dll

  ! Variables
        
  real P0(3), P1(3), lym, epr1, tol0, tet00, tau0, acr0, dleen0
  complex ev, ei0, mi0

  real i_P0(3), i_P1(3), i_lym, i_epr1, i_tol0, i_tet00, i_tau0, i_acr0, i_dleen0
  complex i_ev, i_ei0, i_mi0

     
  ! Body of PointBle_dll
  P0 = i_P0; P1 = i_P1; lym = i_lym; ev = i_ev; ei0 = i_ei0; mi0 = i_mi0; 
  tol0 = i_tol0; tet00 = i_tet00; tau0 = i_tau0; acr0 = i_acr0; dleen0 = i_dleen0

  call PointBle(P0,P1,lym,ev,epr1,ei0,mi0,tol0,tet00,tau0,acr0,dleen0)
  i_epr1 = epr1
  !i_epr1 = 1

end subroutine PointBle_dll

!_________________________________________________________________________________

subroutine PointBle_dll_common(i_P0,i_P1,i_lym,i_ev,i_ei0,i_mi0,i_tol0,i_tet00,i_tau0,i_acr0,i_dleen0,i_Rr0,i_Rr1,i_index,i_pola,i_polb,i_epr1)

! Expose subroutine RCS_functions_cyl to users of this DLL
  !
  !DEC$ ATTRIBUTES DLLEXPORT::PointBle_dll_common

  ! Variables
        
  real P0(3), P1(3), lym, epr1, tol0, tet00, tau0, acr0, dleen0
  real Rr0(3), Rr1(3), pola, polb
  integer index
  complex ev, ei0, mi0

  real i_P0(3), i_P1(3), i_lym, i_epr1, i_tol0, i_tet00, i_tau0, i_acr0, i_dleen0
  real i_Rr0(3), i_Rr1(3), i_pola, i_polb
  integer i_index
  complex i_ev, i_ei0, i_mi0

  
  COMMON /S/ Rr0,Rr1
  common /indd/ index
  common /qa/ Pola,Polb

     
  ! Body of PointBle_dll_common
  P0 = i_P0; P1 = i_P1; lym = i_lym; ev = i_ev; ei0 = i_ei0; mi0 = i_mi0; 
  tol0 = i_tol0; tet00 = i_tet00; tau0 = i_tau0; acr0 = i_acr0; dleen0 = i_dleen0

  call PointBle_comm(P0,P1,lym,ev,epr1,ei0,mi0,tol0,tet00,tau0,acr0,dleen0,S,indd,qa)
  i_epr1 = epr1
  !i_epr1 = 1
  i_Rr0 = Rr0; i_Rr1 = Rr1
  i_index = index
  i_Pola = Pola; i_Polb = Polb

end subroutine PointBle_dll_common