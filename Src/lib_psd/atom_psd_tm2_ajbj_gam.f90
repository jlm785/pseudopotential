!>  Calculates b(j) at fixed gamma and delta
!>  rpsi(k) = r(k)^(l+1) * exp(polyr)
!>  polyr = delta + gamma^2 + sum_j b(j)*r^(2*j+2)
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.013
!>  \date         1990s, 30 June 2021
!>  \copyright    GNU Public License v2

subroutine atom_psd_tm2_ajbj_gam(bjin, rcpn, delta, gama, aj, bj)

! extracted from pseudo2. 30 June 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)          ::  bjin(5)                          !<  right member of the equation
  real(REAL64), intent(in)          ::  rcpn(12)                         !<  r_cut^n
  real(REAL64), intent(in)          ::  delta                            !<  current value of delta
  real(REAL64), intent(in)          ::  gama                             !<  current value of gamma

! output

  real(REAL64), intent(out)         ::  aj(5,5), bj(5)                   !<  aj: auxiliary;  bj: polyr coefficients

! local variables

  integer          ::  info, ipvt(5)

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64


  bj(1) = bjin(1) - delta - gama*rcpn(2)
  bj(2) = bjin(2) - 2*gama*rcpn(1)
  bj(3) = bjin(3) - 2*gama
  bj(4) = bjin(4)
  bj(5) = bjin(5)

  aj(1,1) = rcpn(4)
  aj(1,2) = rcpn(6)
  aj(1,3) = rcpn(8)
  aj(1,4) = rcpn(10)
  aj(1,5) = rcpn(12)

  aj(2,1) = 4*rcpn(3)
  aj(2,2) = 6*rcpn(5)
  aj(2,3) = 8*rcpn(7)
  aj(2,4) = 10*rcpn(9)
  aj(2,5) = 12*rcpn(11)

  aj(3,1) = 4*3*rcpn(2)
  aj(3,2) = 6*5*rcpn(4)
  aj(3,3) = 8*7*rcpn(6)
  aj(3,4) = 10*9*rcpn(8)
  aj(3,5) = 12*11*rcpn(10)

  aj(4,1) = 4*3*2*rcpn(1)
  aj(4,2) = 6*5*4*rcpn(3)
  aj(4,3) = 8*7*6*rcpn(5)
  aj(4,4) = 10*9*8*rcpn(7)
  aj(4,5) = 12*11*10*rcpn(9)

  aj(5,1) = 4*3*2*ONE
  aj(5,2) = 6*5*4*3*rcpn(2)
  aj(5,3) = 8*7*6*5*rcpn(4)
  aj(5,4) = 10*9*8*7*rcpn(6)
  aj(5,5) = 12*11*10*9*rcpn(8)

  call dgefa(aj,5,5,ipvt,info)
  if(info /= 0) stop
  call dgesl(aj,5,5,ipvt,bj,0)

  return

end subroutine atom_psd_tm2_ajbj_gam
