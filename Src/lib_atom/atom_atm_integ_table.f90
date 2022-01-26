!>  Calculates the integral f(x) dx between f_0 and f_j up to f_np
!>  using mostly Bode's rule, Abramowitz-Stegun 25.4.14, complemented by 25.4.15.
!>  For small number of intervals, uses, trapeze or Simpson (25.4.5, 25.4.13)
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.4
!>  \date         30 September 2021.
!>  \copyright    GNU Public License v2


subroutine atom_atm_integ_table(x, f, d, np)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  np                               !<  number of points
  real(REAL64), intent(in)          ::  f(0:np-1)                        !<  function to be integrated
  real(REAL64), intent(in)          ::  d(0:np-1)                        !<  step: d(x) / d (i)

! output

  real(REAL64), intent(out)         ::  x(0:np-1)                        !<  integral

! local variables

  integer                  ::  n                                         !  number of intervlas

  real(REAL64)             ::  x1
  integer                  ::  n4, nn

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64

! counter

  integer                  ::  i, j

  n = np - 1


  x(0) = ZERO
  if(n >= 1) then
    x(1) = (f(0)*d(0) + f(1)*d(1)) / 2
  endif
  if(n >= 2) then
    x(2) = (f(0)*d(0) + 4*f(1)*d(1) + f(2)*d(2)) / 3
  endif
  if(n >= 3) then
    x1 = (f(0)*d(0) + 3*f(1)*d(1) + 3*f(2)*d(2) + f(3)*d(3)) / 8
    x(3) = 3*x1
  endif
  if(n >= 4) then
    x1 = (7*f(0)*d(0) + 32*f(1)*d(1) + 12*f(2)*d(2) +                    &
            32*f(3)*d(3) + 7*f(4)*d(4)) / 90
    x(4) = 4*x1
  endif
  if(n >= 5) then
    x1 = (19*f(0)*d(0) + 75*f(1)*d(1) + 50*f(2)*d(2) +                   &
            50*f(3)*d(3) + 75*f(4)*d(4) + 19*f(5)*d(5)) / 288
    x(5) = 5*x1
  endif
  if(n >= 6) then
    x1 = (f(3)*d(3) + 3*f(4)*d(4) + 3*f(5)*d(5) + f(6)*d(6)) / 8
    x(6) = x(3) + 3*x1
  endif
  if(n >= 7) then
    x1 = (f(4)*d(4) + 3*f(5)*d(5) + 3*f(6)*d(6) + f(7)*d(7)) / 8
    x(7) = x(4) + 3*x1
  endif
  if(n >= 8) then
    x1 = (7*f(4)*d(4) + 32*f(5)*d(5) + 12*f(6)*d(6) +                    &
            32*f(7)*d(7) + 7*f(8)*d(8)) / 90
    x(8) = x(4) + 4*x1
  endif
  if(n >= 9) then
    x1 = (19*f(4)*d(4) + 75*f(5)*d(5) + 50*f(6)*d(6) +                   &
            50*f(7)*d(7) + 75*f(8)*d(8) + 19*f(9)*d(9)) / 288
    x(9) = x(4) + 5*x1
  endif
  if(n >= 10) then
    x1 = (19*f(5)*d(5) + 75*f(6)*d(6) + 50*f(7)*d(7) +                   &
            50*f(8)*d(8) + 75*f(9)*d(9) + 19*f(10)*d(10)) / 288
    x(10) = x(5) + 5*x1
  endif
  if(n >= 11) then
    x1 = (f(8)*d(8) + 3*f(9)*d(9) + 3*f(10)*d(10) + f(11)*d(11)) / 8
    x(11) = x(8) + 3*x1
  endif
  if(n >= 12) then
    n4 = n / 4
    do i = 3,n4
      nn = 4*i
      x1 = (7*f(nn-4)*d(nn-4) + 32*f(nn-3)*d(nn-3) + 12*f(nn-2)*d(nn-2)    &
             + 32*f(nn-1)*d(nn-1) + 7*f(nn)*d(nn)) / 90
      x(nn) = x(nn-4) + 4*x1
      do j = 1,3
        nn = 4*i+j
        if(nn > n) exit
        x1 = (19*f(nn-5)*d(nn-5) + 75*f(nn-4)*d(nn-4) + 50*f(nn-3)*d(nn-3)   &
           + 50*f(nn-2)*d(nn-2) + 75*f(nn-1)*d(nn-1) + 19*f(nn)*d(nn)) / 288
        x(nn) = x(nn-5) + 5*x1
      enddo
    enddo
  endif

  return

end subroutine atom_atm_integ_table
