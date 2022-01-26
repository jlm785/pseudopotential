!>  splift fits an interpolating cubic spline to the n data point
!>  given in x and y and returns the first and second derivatives
!>  in yp and ypp.  The resulting spline (defined by x, y, and
!>  ypp) and its first and second derivatives may then be
!>  evaluated using splint.  The spline may be integrated using
!>  spliq.
!>  
!>  The code is heavily based on, and compatible with, the 1978 subroutine 
!>  written by Rondall E, Jones at Sandia. 
!>
!>  \author       Jose Luis Martins
!>  \version      5.805
!>  \date         9 May 2018
!>  \copyright    GNU Public License v2

subroutine splift (x,y,yp,ypp,n,w,ierr,isx,a1,b1,an,bn)

! copyright José Luís Martins,  INESC-MN.


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter          :: iowrite = 6

! input

  integer, intent(in)               ::  n                                !<  number of data points that define the spline. n > 1
  real(REAL64), intent(in)          ::  x(n)                             !<  array of abscissas (in increasing order) that define the spline.
  real(REAL64), intent(in)          ::  y(n)                             !<  array of ordinates that define the spline.

  integer, intent(in)               ::  isx                              !<  must be 0 on the initial call to splift.  If x and w, a_i, b_i did not change after previous call, may be set to 1.

  real(REAL64), intent(in)          ::  a1                               !<  specify the end conditions for the spline:
  real(REAL64), intent(in)          ::  b1                               !<  ypp(1) = a1*ypp(2) + b1, ypp(n) = an*ypp(n-1) + bn, with abs(a1) < 1  and  abs(an) < 1.
  real(REAL64), intent(in)          ::  an                               !<  The smoothest spline is obtained by a1=b1=an=bn=0.
  real(REAL64), intent(in)          ::  bn                               !<  If the data is to be extrapolated taking a1=an=0.5 and b1=bn=0 may yield better results.
                                                                         !<  a spline that has a given first derivative yp1 at x(1) and ypn at y(n) may be defined by using:
                                                                         !<  a1=-0.5, b1= 3.0*((y(2)-y(1))/(x(2)-x(1))-yp1)/(x(2)-x(1)),
                                                                         !<  an=-0.5, bn=-3.0*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)/(x(n)-x(n-1)

! output

  real(REAL64), intent(out)         ::  yp(n)                            !<  array of first derivatives of the spline (at the x(i)).
  real(REAL64), intent(out)         ::  ypp(n)                           !<  array of second derivatives that define the spline (at the x(i)).

  integer, intent(out)              ::  ierr                             !<  status code. 1: all normal; 2: n < 4; 3: abscissas are not increasing.

! input and output

  real(REAL64), intent(inout)       ::  w(n,3)                           !<  previous/next LU decomposition

! local variables

  real(REAL64)               ::  dold, dnew

! counter
  
  integer                    ::  i, j

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64

  ierr = 1
  if (n < 4) then
    ierr = 2
!   defaults to linear interpolation
    if(n > 0) then
      do i = 1,n
        yp(i) = ZERO
        ypp(i) = ZERO
      enddo
    endif
  else
    if (isx == 0) then

!     first time or new data

      do i = 2,n
        if (x(i)-x(i-1) <= ZERO) then
          ierr = 3
        endif
      enddo
      
      if(ierr /= 3) then

!       define the tridiagonal matrix

        w(1,3) = x(2)-x(1)
        do i = 2,n-1
          w(i,2) = w(i-1,3)
          w(i,3) = x(i+1)-x(i)
          w(i,1) = 2*(w(i,2)+w(i,3))
        enddo
        w(1,1) = 4*UM
        w(1,3) =-4*a1
        w(n,1) = 4*UM
        w(n,2) =-4*an

!       LU decomposition

        do i = 2,n
          w(i-1,3) = w(i-1,3)/w(i-1,1)
          w(i,1) = w(i,1) - w(i,2)*w(i-1,3)
        enddo
      endif

    endif
      
    if(ierr /= 3) then

!     define *constant* vector

      ypp(1) = 4*b1
      dold = (y(2)-y(1))/w(2,2)
      do i = 2,n-2
        dnew   = (y(i+1) - y(i))/w(i+1,2)
        ypp(i) = 6*(dnew - dold)
        yp(i)  = dold
        dold = dnew
      enddo

      dnew = (y(n)-y(n-1))/(x(n)-x(n-1))
      ypp(n-1) = 6*(dnew - dold)
      ypp(n) = 4*bn
      yp(n-1)= dold
      yp(n) = dnew

!     forward substitution

      ypp(1) = ypp(1)/w(1,1)
      do i = 2,n
        ypp(i) = (ypp(i) - w(i,2)*ypp(i-1))/w(i,1)
      enddo

!     backward substitution

      do j = 1,n-1
        i = n-j
        ypp(i) = ypp(i) - w(i,3)*ypp(i+1)
      enddo

!     compute first derivatives

      yp(1) = (y(2)-y(1))/(x(2)-x(1)) - (x(2)-x(1))*(2*ypp(1) + ypp(2))/6
      do i = 2,n-1
         yp(i) = yp(i) + w(i,2)*(ypp(i-1) + 2*ypp(i))/6
      enddo
      yp(n) = yp(n) + (x(n)-x(n-1))*(ypp(n-1) + 2*ypp(n))/6

      ierr = 1

    endif

  endif

  return

end subroutine splift
