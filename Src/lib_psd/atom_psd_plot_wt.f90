!>  Prints the Fourier transform of the wave-functions
!>  Plot range (q=0-12) is hardcoded...
!>
!>  \author       Sverre Froyen,  Norm Troullier, Jose Luis Martins
!>  \version      6.0.6
!>  \date         1980s and 1990s, 30 June 2021, 3 November 2021.
!>  \copyright    GNU Public License v2

subroutine atom_psd_plot_wt(nr, r, drdi, vd, qmax, qstep, l, ist, itplot, mxdnr)

! ***********************************************************
! *                                                         *
! *    This is a plotting routine; the user should adjust   *
! *  for their own needs.  The result                       *
! *  is then printed to the current plot.dat file (unit=3)  *
! *  for later plotting of the data.  A marker (marker fw#) *
! *  is placed at the end of each set of data.              *
! *                                                         *
! ***********************************************************

!jlm    check dimensions 21 January 2007
! do not plot zeroes decision in calling subroutine. 3 November 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension  of the number of radial points
  integer, intent(in)               ::  itplot                           !<  default plot file

  integer, intent(in)               ::  nr                               !<  number of radial points
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r / d i
  real(REAL64), intent(in)          ::  vd(mxdnr)                        !<  potential

  real(REAL64), intent(in)          ::  qmax                             !<  maximum value of q in plot
  real(REAL64), intent(in)          ::  qstep                            !<  step for the plot

  integer, intent(in)               ::  l                                !<  angular momentum
  integer, intent(in)               ::  ist                              !<  sign/scaling

! local allocatable arrays

  real(REAL64), allocatable         ::  qp(:), vql(:)

! local variables

  integer         ::  nq

! parameters

  real(REAL64), parameter       ::  ZERO = 0.0_REAL64

! counter

  integer             ::  j


  nq = nint(qmax/qstep)

  allocate(qp(nq),vql(nq))

  do j = 1,nq
    vql(j) = ZERO
    qp(j) = qstep*j
  enddo

  call ft_radial(nr-1, nq-1, r, drdi, qp, vd, vql, l)

  call atom_plot_one(1, nq, qp, vql, ZERO, l, ist, 'fw  ', '    ',       &
              .FALSE., itplot, mxdnr)

  deallocate(qp,vql)

  return

end subroutine atom_psd_plot_wt
