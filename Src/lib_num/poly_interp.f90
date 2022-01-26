!>  Performs a Lagrange interpolation of a function and gives its
!>  Derivatives.
!>  Based on D. B. Hunter, The Computer Journal, 3, 270 (1961)
!>  with the small difference algorithm of Numerical Recipes Eq. 3.1.5
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.6
!>  \date         5 July 2021, 6 November 2021.
!>  \copyright    GNU Public License v3

subroutine poly_interp(y, dy, xin, yin, n, nd)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  n                                !<  order of interpolation
  integer, intent(in)               ::  nd                               !<  order of derivative
  real(REAL64), intent(in)          ::  xin(0:n)                         !<  grid points
  real(REAL64), intent(in)          ::  yin(0:n)                         !<  f(xin)

! output

  real(REAL64), intent(out)         ::  y(0:nd)                          !<  f(0) interpolated value at x=0
  real(REAL64), intent(out)         ::  dy(0:nd)                         !<  error estimate (last correction)

! work arrays

  real(REAL64), allocatable         ::  cmi(:,:), dmi(:,:)           !  difference arrays of Neville's algorithm
  real(REAL64), allocatable         ::  xnum(:)

! local variables

  integer         ::  ns
  real(REAL64)    ::  dp, dm

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64

!      counters

  integer     ::  i, m, k, kmax

  allocate(cmi(0:nd,0:n),dmi(0:nd,0:n))
  allocate(xnum(0:nd))

  ns = 0
  do i = 0,n
    if( abs(xin(i)) < abs(xin(ns)) ) ns = i
  enddo

  y(0) = yin(ns)
  ns = ns - 1

  do i = 0,n
    cmi(0,i) = yin(i)
    dmi(0,i) = yin(i)
  enddo

  if(nd > 0) then
    do k = 1,nd
      y(k) = ZERO
    enddo
    do i = 0,n
      do k = 1,nd
        cmi(k,i) = ZERO
        dmi(k,i) = ZERO
      enddo
    enddo
  endif

  do m = 1,n

    kmax = min(m,nd)
    do i = 0,n-m
      dm = xin(i)
      dp = xin(i+m)
      do k = 0,kmax
        xnum(k) = (dmi(k,i) - cmi(k,i+1)) / (dp - dm)
      enddo
      dmi(0,i) = dp*xnum(0)
      cmi(0,i) = dm*xnum(0)
      if(nd > 0) then
        do k = 1,kmax
          dmi(k,i) = dp*xnum(k) - k*xnum(k-1)
          cmi(k,i) = dm*xnum(k) - k*xnum(k-1)
        enddo
      endif
    enddo

    if (2*ns+1 < n - m) then
      do k = 0,kmax
        dy(k) = cmi(k,ns+1)
      enddo
    else
      do k = 0,kmax
        dy(k) = dmi(k,ns)
      enddo
      ns = ns - 1
    endif

    do k = 0,kmax
      y(k) = y(k) + dy(k)
    enddo

  enddo

  return

end subroutine poly_interp

