!>  Calculates the integral f(x) dx between f_0 and f_np
!>  using mostly Bode's rule, Abramowitz-Stegun 25.4.14, complemented by 25.4.15.
!>  For small number of intervals, uses, trapeze or Simpson (25.4.5, 25.4.13)
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.4
!>  \date         29 September 2021.
!>  \copyright    GNU Public License v2


subroutine atom_atm_integ(x, f, d, np)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  np                               !<  number of points
  real(REAL64), intent(in)          ::  f(0:np-1)                        !<  function to be integrated
  real(REAL64), intent(in)          ::  d(0:np-1)                        !<  step: d(x) / d (i)

! output

  real(REAL64), intent(out)         ::  x                                !<  integral

! local variables

  integer                  ::  n                                         !  number of intervlas

  real(REAL64)             ::  x1, x2
  integer                  ::  n4, n5, nn

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64

! counter

  integer                  ::  i

  n = np - 1

  if(n < 1) then
    x = ZERO
  elseif(n == 1) then
    x = (f(0)*d(0) + f(1)*d(1)) / 2
  elseif(n == 2) then
    x = (f(0)*d(0) + 4*f(1)*d(1) + f(2)*d(2)) / 3
  elseif(n == 3) then
    x = (f(0)*d(0) + 3*f(1)*d(1) + 3*f(2)*d(2)+f(3)*d(3)) / 8
    x = 3*x
  elseif(n == 6) then
    x = (f(0)*d(0) + 3*f(1)*d(1) + 3*f(2)*d(2) + 2*f(3)*d(3) +           &
             3*f(4)*d(4) + 3*f(5)*d(5) + f(6)*d(6)) / 8
    x = 3*x
  elseif(n == 7) then
    x1 = (f(0)*d(0) + 3*f(1)*d(1) + 3*f(2)*d(2)+f(3)*d(3)) / 8
    x2 = (7*f(3)*d(3) + 32*f(4)*d(4) + 12*f(5)*d(5) + 32*f(6)*d(6) +     &
              7*f(7)*d(7)) / 90
    x = 3*x1 + 4*x2
  elseif(n == 11) then
    x1 = (f(0)*d(0) + 3*f(1)*d(1) + 3*f(2)*d(2)+f(3)*d(3)) / 8
    x2 = (7*f(3)*d(3) + 32*f(4)*d(4) + 12*f(5)*d(5) + 32*f(6)*d(6) +     &
              14*f(7)*d(7) + 32*f(8)*d(8) + 12*f(9)*d(9) +               &
              32*f(10)*d(10)+7*f(11)*d(11)) / 90
    x = 3*x1 + 4*x2
  else
    n5 = mod(n,4)
    n4 = (n - 5*n5) / 4
    x1 = ZERO
    x2 = ZERO
    if(n4 > 0) then
      do i = 1,n4
        nn = 4*i-4
        x1 = x1 + 32*f(nn+1)*d(nn+1) + 12*f(nn+2)*d(nn+2) +              &
                  32*f(nn+3)*d(nn+3) + 14*f(nn+4)*d(nn+4)
      enddo
      x1 = x1 + 7*(f(0)*d(0) - f(4*n4)*d(4*n4))
      x1 = (4*x1) / 90
    endif
    if(n5 > 0) then
      do i = 1,n5
        nn = 4*n4 + 5*i - 5
        x2 = x2 + 75*f(nn+1)*d(nn+1) + 50*f(nn+2)*d(nn+2) +             &
                  50*f(nn+3)*d(nn+3) + 75*f(nn+4)*d(nn+4) +             &
                  38*f(nn+5)*d(nn+5)
      enddo
      x2 = x2 + 19*(f(4*n4)*d(4*n4) - f(4*n4+5*n5)*d(4*n4+5*n5))
      x2 = (5*x2) / 288
    endif
    x = x1 + x2
  endif

  return

end subroutine atom_atm_integ
