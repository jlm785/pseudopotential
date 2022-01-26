!>  Prints the Fourier transform of the potential.
!>
!>  \author       Sverre Froyen,  Norm Troullier, Jose Luis Martins
!>  \version      6.013
!>  \date         1980s and 1990s, 30 June 2021
!>  \copyright    GNU Public License v2

subroutine atom_psd_plot_ran(r, nr, vd, qmax, qstep, lp1, zion, itplot, mxdnr)

! ***********************************************************
! *                                                         *
! *    This is a plotting routine; the user should adjust   *
! *  for their own needs.  The potential is fitted with a   *
! *  second degree polynomial, which is muliplied with the  *
! *  appropriate functions and then integrated by parts     *
! *  to find the fourier transform.  The result is then     *
! *  printed to the current plot.dat file (unit=3) for      *
! *  later plotting.  A marker(marker fn#) is placed at     *
! *  the end of each set of data.                           *
! *                                                         *
! ***********************************************************

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension  of the number of radial points
  integer, intent(in)               ::  itplot                           !<  default plot file

  integer, intent(in)               ::  nr                               !<  number of radial points
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid
  real(REAL64), intent(in)          ::  vd(mxdnr)                        !<  potential

  real(REAL64), intent(in)          ::  qmax                             !<  maximum value of q in plot
  real(REAL64), intent(in)          ::  qstep                            !<  step for the plot

  integer, intent(in)               ::  lp1                              !<  l+1
  real(REAL64), intent(in)          ::  zion                             !<  effective Z

! local allocatable arrays

  real(REAL64), allocatable         ::  a(:), b(:), c(:)
  real(REAL64), allocatable         ::  qp(:), vql(:)

! local variables

  real(REAL64)    ::  rm, r0, rp
  real(REAL64)    ::  vm, v0, vp
  real(REAL64)    ::  d1, d2, d3
  real(REAL64)    ::  q, q2
  integer         ::  nq

! parameters

  real(REAL64), parameter       ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64

! counter

  integer             ::  j, k


!  The potential times r is fitted to the polynominal
!  a + bx + cx^2 at every other point.

!  Modified JLM, do 130 k = 2,nr-1,2 12 April 2018

!  Should use ft_local in some future.  If ain't broken...

  allocate(a(nr),b(nr),c(nr))

  rm = ZERO
  vm = 2*zion
  do k = 2,nr-1,2
    r0 = r(k)
    v0 = r0*vd(k) + 2*zion
    rp = r(k+1)
    vp = rp*vd(k+1) + 2*zion
    d1 = ONE/((rp-rm)*(r0-rm))
    d2 = ONE/((rp-r0)*(rm-r0))
    d3 = ONE/((r0-rp)*(rm-rp))
    a(k) = vm*d1 + v0*d2 + vp*d3
    b(k) =-vm*(r0+rp)*d1 - v0*(rm+rp)*d2 - vp*(rm+r0)*d3
    c(k) = vm*r0*rp*d1 + v0*rm*rp*d2 + vp*rm*r0*d3
    rm = rp
    vm = vp
  enddo

!  Find the fourier transform q^2/4pi/zion*vql. Everything is
!  rescaled  by zion.

  nq = nint(qmax/qstep)

  allocate(qp(nq),vql(nq))

  do j = 1,nq
    q = j*qstep
    qp(j) = q
    q2 = q*q
    vql(j) = zero
    rm = zero
    do k = 2,nr-1,2
      rp = r(k+1)
      vql(j) = vql(j) + (2*a(k)*rp+b(k))/q*sin(q*rp)                     &
                      - ((a(k)*rp+b(k))*rp+c(k)-2*a(k)/q2)*cos(q*rp)     &
                      - (2*a(k)*rm+b(k))/q*sin(q*rm)                     &
                      + ((a(k)*rm+b(k))*rm+c(k)-2*a(k)/q2)*cos(q*rm)
      rm = rp
    enddo
    vql(j) = vql(j)/2/zion - ONE
  enddo


  call atom_plot_one(1, nq, qp, vql, ZERO, lp1, 1, 'fn  ', '    ',       &
              .FALSE., itplot, mxdnr)


!  Print out the transforms( really q^2/(4pi*zion)*v(q) ) to
!  the current plot.dat file (unit=3) for latter plotting.
!
!   do j = 1,nq
!     write(itplot,'(1x,f12.5,3x,f15.8)') j*qstep, vql(j)
!   enddo
!
!   if(lp1 < 10) then
!     write(itplot,'(1x,"marker fn",i1)') lp1
!   elseif(lp1 < 100) then
!     write(itplot,'(1x,"marker fn",i2)') lp1
!   else
!     write(itplot,'(1x,"marker fnl")')
!   endif

  deallocate(a,b,c)
  deallocate(qp,vql)

  return

end subroutine atom_psd_plot_ran
