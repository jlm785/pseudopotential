!>  Calculates b(j) at fixed delta and gamma = 0
!>  rpsi(k) = r(k)^(l+1) * exp(polyr)
!>  polyr = delta + sum_j b(j)*r^(2*j)
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.013
!>  \date         1990s, 30 June 2021
!>  \copyright    GNU Public License v2

subroutine atom_psd_tm2_ajbj(bjin, rcpn, delta, aj, bj)

! extracted from pseudo2. 30 June 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  real(REAL64), intent(in)         ::  bjin(5)
  real(REAL64), intent(in)         ::  rcpn(12)                         !<  r_cut^n
  real(REAL64), intent(in)         ::  delta


! output

  real(REAL64), intent(out)         ::  aj(5,5), bj(5)

! local variables

  integer          ::  info, ipvt(5)

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64

! counter

  integer      ::  j


  do j = 1,5
    bj(j) = bjin(j)
  enddo
  bj(1) = bjin(1) - delta

  aj(1,1) = rcpn(2)
  aj(1,2) = rcpn(4)
  aj(1,3) = rcpn(6)
  aj(1,4) = rcpn(8)
  aj(1,5) = rcpn(10)

  aj(2,1) = 2*rcpn(1)
  aj(2,2) = 4*rcpn(3)
  aj(2,3) = 6*rcpn(5)
  aj(2,4) = 8*rcpn(7)
  aj(2,5) = 10*rcpn(9)

  aj(3,1) = 2*ONE
  aj(3,2) = 4*3*rcpn(2)
  aj(3,3) = 6*5*rcpn(4)
  aj(3,4) = 8*7*rcpn(6)
  aj(3,5) = 10*9*rcpn(8)

  aj(4,1) = ZERO
  aj(4,2) = 4*3*2*rcpn(1)
  aj(4,3) = 6*5*4*rcpn(3)
  aj(4,4) = 8*7*6*rcpn(5)
  aj(4,5) = 10*9*8*rcpn(7)

  aj(5,1) = ZERO
  aj(5,2) = 4*3*2*1*ONE
  aj(5,3) = 6*5*4*3*rcpn(2)
  aj(5,4) = 8*7*6*5*rcpn(4)
  aj(5,5) = 10*9*8*7*rcpn(6)

  call dgefa(aj,5,5,ipvt,info)
  if(info /= 0) stop
  call dgesl(aj,5,5,ipvt,bj,0)

  return

end subroutine atom_psd_tm2_ajbj
