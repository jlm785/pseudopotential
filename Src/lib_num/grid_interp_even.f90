!>  Interpolates an even function, f(-x) = f(x), from one grid to another
!>  grid using Lagrange interpolation (grid_interp) on the input grid
!>  extended with the symmetric points -x of its first points.
!>  Returns also the derivatives of the function on the new grid
!>  and an estimate of the respective errors.
!>
!>  \author       Jose Luis Martins, claude
!>  \version      6.1.0
!>  \date         8 November 2026.
!>  \copyright    GNU Public License v2


subroutine grid_interp_even(n, nd, nin, xin, fin, nout, xout, fout, dymax, ifail, mxdout)

! The input grid must be ordered with xin(1) >= 0.  A point at x = 0 is used
! if present, but it is not needed, so the value at the origin (often an
! undefined 0/0, for example rho = cdc/r**2) can be calculated from the
! other points and their symmetric points.

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdout                           !<  dimension of points in the output (new) grid

  integer, intent(in)               ::  n                                !<  order of Lagrange interpolation
  integer, intent(in)               ::  nd                               !<  order of derivatives

  integer, intent(in)               ::  nin                              !<  number of points in the input (old) grid
  real(REAL64), intent(in)          ::  xin(nin)                         !<  points on the input grid, ordered, xin(1) >= 0
  real(REAL64), intent(in)          ::  fin(nin)                         !<  f(x(i))

  integer, intent(in)               ::  nout                             !<  number of points in the output (new) grid
  real(REAL64), intent(in)          ::  xout(mxdout)                     !<  points on the output (new) grid

! output

  real(REAL64), intent(out)         ::  fout(mxdout,0:nd)                !<  fout(:,0): value of the function on the new grid; fout(:,n): n-th derivative
  real(REAL64), intent(out)         ::  dymax(0:nd)                      !<  estimate of the accuracy (from n-1 interpolation)
  integer, intent(out)              ::  ifail                            !<  ifail =/= 0 indicates failure of the process.

! allocatable arrays

  real(REAL64), allocatable         ::  xs(:)                            !  input grid with the symmetric points
  real(REAL64), allocatable         ::  fs(:)                            !  function on that grid

! local variables

  integer        ::  ip                                                  !  first point with x > 0
  integer        ::  nsym                                                !  number of symmetric points

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64

! counters

  integer     ::  i


  ifail = 0

  if(nin < 1) then
    ifail = 1

    return

  endif

  if(xin(1) < ZERO) then
    write(6,*) '  grid_interp_even,  negative xin(1) = ', xin(1)
    ifail = 6

    return

  endif

  ip = 1
  if(xin(1) == ZERO) ip = 2

! n+2 symmetric points cover the interpolation stencil of points near the origin

  nsym = min(n+2, nin-ip+1)

  allocate(xs(nin+nsym), fs(nin+nsym))

  do i = 1,nsym
    xs(i) = -xin(ip+nsym-i)
    fs(i) = fin(ip+nsym-i)
  enddo
  do i = 1,nin
    xs(nsym+i) = xin(i)
    fs(nsym+i) = fin(i)
  enddo

  call grid_interp(n, nd, nin+nsym, xs, fs, nout, xout, fout, dymax, ifail, mxdout)

  deallocate(xs, fs)

  return

end subroutine grid_interp_even
