!>  Calculates the spherical bessel function of first kind j_n(x)
!>  Formulas 10.1.2, 10.1.19 and 10.1.59 of Abramowitz and Stegun
!>
!>  \author       SIESTA developers et al., Jose Luis Martins
!>  \version      5.805
!>  \date         19 April 2018
!>  \copyright    GNU Public License v2

! Copyright (C) 2019  Jose Luis Martins

! This file is free software: you can redistribute it and/or modify
! it under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This file is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.

! You should have received a copy of the GNU Lesser General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

function sbessj(n,x)

! Written 19 April 2018 from old atom code, Siesta code, and NumSBF code.
! Old code had absolute precision but not relative precision when
! compared with the results of Mathematica.

! Relative errors are only relevant for x ~ 0.75 n and n > 30
! If you want those values use backward recurrence.  
! See for example L.-Wu Cai, Comp. Phys. Comm. 182, 663 2011.

! Copyright José Luís Martins, INESC-MN

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  n                                !<  n >= 0 order of function
  real(REAL64), intent(in)          ::  x                                !<  argument

! output

  real(REAL64)                      ::  sbessj                           !< result

! local variables

  real(REAL64)               ::  by, bym, byp, ux                   !  recurrence variables
  real(REAL64)               ::  pref                               !  prefactor of series expansion
  real(REAL64)               ::  x2                                 !  0.5*x^2
  real(REAL64)               ::  sumx, fac                          !  series expansion
  logical                    ::  fail
  real(REAL64)               ::  expse                              !  x < expse uses series expansion

! counter

  integer                    ::  j

! parameters

  real(REAL64), parameter           ::  ZERO = 0.0_REAL64
  real(REAL64), parameter           ::  UM = 1.0_REAL64
  real(REAL64), parameter           ::  EPS = 1.0 / 10.0**15

! spherical bessel function of the first kind

  if(n < 0) then
    write(6,*) '  sbessj: negative order, n = ', n

    stop

! remove the following test to see the results in the region of low accuracy

  elseif(n > 50) then
    if(x < n-10 .and. x > 5.5*sqrt(UM*n)) then
      write(6,*) '  sbessj: precision would be low, n, x = ', n,x

      stop

    endif

  endif

! old code
! expse = 0.001
! SIESTA
! expse = max(1,2*n+1)
! just before maximum
! expse = (n + 0.5) + 0.8086165*(n+0.5)**0.33333333 - 3.1416/4.0

! from NumSBT, JD Talman, Comp. Phys. Comm., 180, 332.
! expse = 0.75 * n
  expse = 0.78 * n

  if(abs(x) > expse) then

!   recursion formula
    if(n == 0) then
      sbessj = sin(x) / x
    elseif(n == 1) then
      sbessj = (sin(x)/x - cos(x)) / x
    else
      ux = UM / x
      bym = sin(x) * ux
      by = (bym - cos(x)) * ux
      do j = 1,n-1
        byp = (2*j+1)*ux*by - bym
        bym = by
        by = byp
      enddo

      sbessj = by

    endif
        
  else

!   series expansion

    pref = UM
    if(n > 0) then
      do j = 1,n
        pref = pref*x/(2*j+1)
      enddo
    endif

!   maximum j can be estimated from Stirling formula

    x2 = x*x/2
    sumx = UM
    fac = UM
    fail = .TRUE.

    do j = 1,20+2*n
      fac = -fac*x2 / (j*(2*n+2*j+1))
      sumx = sumx + fac
      if(abs(fac) < EPS) then
        fail = .FALSE.             
        exit
      endif
    enddo

    if(fail) then
      write(6,*) '  sbessj:  did not converge, n, x = ', n, x

      stop

    endif

    sbessj = sumx*pref

  endif
  
  

  return

end function sbessj
