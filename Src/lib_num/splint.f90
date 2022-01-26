!>  splint evaluates a cubic spline and its first and second
!>  derivatives at the abscissas in xi.  The spline (which
!>  is defined by x, y, and ypp) may have been determined by
!>  splift or any other fitting routine that
!>  provides second derivatives.
!>
!>  The code is based on, and compatible with, the 1978 subroutine
!>  written by Rondall E, Jones at Sandia.  The linear upward
!>  search was replaced by the hunt algorithm from Numerical Recipes.
!>
!>  \author       Jose Luis Martins
!>  \version      5.805
!>  \date         9 May 2018
!>  \copyright    GNU Public License v2

subroutine splint (x,y,ypp,n,xi,yi,ypi,yppi,ni,kerr)

! copyright José Luís Martins,  INESC-MN.


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter          :: iowrite = 6

! input

  integer, intent(in)               ::  n                                !<  number of data points that define the spline. n > 1
  real(REAL64), intent(in)          ::  x(n)                             !<  array of abscissas (in increasing order) that define the spline.
  real(REAL64), intent(in)          ::  y(n)                             !<  array of ordinates that define the spline.
  real(REAL64), intent(in)          ::  ypp(n)                           !<  array of second derivatives that define the spline.

  integer, intent(in)               ::  ni                               !<  the number of abscissas at which the spline is to be evaluated. ni > 0
  real(REAL64), intent(in)          ::  xi(ni)                           !<  array of abscissas (in arbitrary order) at which the spline is to be evaluated.

! output

  real(REAL64), intent(out)         ::  yi(ni)                           !<  array of values of the spline (ordinates) at xi
  real(REAL64), intent(out)         ::  ypi(ni)                          !<  array of values of the first derivative of spline at xi
  real(REAL64), intent(out)         ::  yppi(ni)                         !<  array of values of the second derivative of spline at xi

  integer, intent(out)              ::  kerr                             !<  status code. 1: all normal; 2: extrapolations occurred; 3: ni < 1, no interpolation.

! local variables

  real(REAL64)               ::  xx, xxold                               !  current and previous abscissa to be evaluated, xi(k), xi(k-1)
  integer                    ::  i                                       !  i is current index into x array.
  integer                    ::  k                                       !  k is index on value of xi being worked on.
  integer                    ::  il, ir                                  !  bracket interval for bissection
  integer                    ::  jh, inc                                 !  try index and increment for hunt phase

  real(REAL64)               ::  h, h2
  real(REAL64)               ::  xl, xl2, xl3
  real(REAL64)               ::  xr, xr2, xr3

! counter

  integer                    ::  j

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64

! check input

  if (ni <= 0) then
    write(iowrite,*) '   in splint,  the requested number of ',          &
              'interpolating points was not positive'
    kerr = 3
  else
    kerr = 1

    do k = 1,ni
      xx = xi(k)

!     checks for extrapolation

      if (xx < x(1)) then
        kerr = 2
        i = 1
      elseif (xx > x(n)) then
        kerr = 2
        i  = n-1
      else

!       use past information

        if(k == 1) then
          il = 1
          ir = n
        else
          if(xxold < xx) then
            il = i
            ir = n
          else
            il = 1
            ir = min(n,i+1)
          endif
        endif

!       hunt phase

        inc = 1
        do j = 1,n+10
          jh = il + inc
          if(jh > ir) then
            exit
          else
            if(xx > x(jh)) then
              il = jh
              inc = 2*inc
            else
              ir = jh
              exit
            endif
          endif
        enddo

!       bisection search

        do j = 1,n+10

          i  = (il+ir)/2

          if(i == il) exit

          if(xx == x(i)) exit

          if(xx < x(i)) then
            ir = i
          else
            il = i
          endif
        enddo

      endif

!     just being paranoid

      i = min(i,n-1)

      xxold = xx

!     interpolate

      h  = x(i+1) - x(i)
      if(h == ZERO) then
        write(iowrite,*) '  in splint   identical abscissas'
        stop
      endif
      h2 = h*h
      xr = (x(i+1)-xx)/h
      xr2= xr*xr
      xr3= xr*xr2
      xl = (xx-x(i))/h
      xl2= xl*xl
      xl3= xl*xl2
      yi(k) = y(i)*xr + y(i+1)*xl - h2*(ypp(i)*(xr-xr3) + ypp(i+1)*(xl-xl3))/6
      ypi(k) = (y(i+1)-y(i))/h + h*(ypp(i)*(UM-3*xr2) - ypp(i+1)*(UM-3*xl2))/6
      yppi(k) = ypp(i)*xr + ypp(i+1)*xl

    enddo
  endif

  return

end subroutine splint
