!>  Retuns the values of delta, b(j) given a fixed value of gamma.
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.0.6
!>  \date         1980s and 1990s, 30 June 2021, 30 October 2021.
!>  \copyright    GNU Public License v2

subroutine atom_psd_tm2_gamfn(gama, v0pp, iflag,                         &
       jrc, r, drdi, bjin, rcpn, cdrc, lo, rpsi, delta, bj,              &
       mxdnr)

! *********************************************************
! *                                                       *
! *  njtj                                                 *
! *   Retuns the values of delta, alpha, alpha1, alpha2,  *
! *   alpha3, and alpha4 given a fixed value of gamma.    *
! *  njtj                                                 *
! *                                                       *
! *********************************************************

! use bjin, instead of recalculating. 1 July 2021. JLM
! lp -> lo, bkrk, 30 October 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension of radial points

  integer, intent(in)               ::  jrc                              !<  r(jrc) = r_cut

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  real(REAL64), intent(in)          ::  bjin(5)                          !<  right member of the equation
  real(REAL64), intent(in)          ::  rcpn(12)                         !<  r_cut^n

  real(REAL64), intent(in)          ::  cdrc                             !<  int_0^r_cut rho(r) dr

  integer, intent(in)               ::  lo                               !<  angular momentum l

  real(REAL64), intent(in)          ::  gama                             !<  value of gamma

! output

  real(REAL64), intent(out)         ::  v0pp                             !<  v''(0)

  real(REAL64), intent(out)         ::  delta                            !<  value of delta
  real(REAL64), intent(out)         ::  bj(5)                            !<  polyr coefficients

  integer, intent(out)              ::  iflag                            !<  iflag = 0: success

! input and output

  real(REAL64), intent(inout)       ::  rpsi(jrc)                        !<  r*psi, modified up to jrc.

! allocatable arrays

  real(REAL64), allocatable         ::  aj(:,:)

! local variables

  real(REAL64)        ::  cdps
  real(REAL64)        ::  fdold, fdnew, ddelta
  real(REAL64)        ::  bkrk(0:6)

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  ERRMIN = 1.0E-12_REAL64
  integer, parameter         ::  NTRY = 100

! counters

  integer  ::  j, k


  allocate(aj(5,5))

  delta = ZERO

  call atom_psd_tm2_ajbj_gam(bjin, rcpn, ZERO, gama, aj, bj)

! start iteration loop to find delta(with gama fixed)

  iflag = 1
  do j = 1,NTRY

!   generate pseudo wavefunction-note missing factor exp(delta)



    bkrk(0) = ZERO
    bkrk(1) = gama
    do k = 2,6
      bkrk(k) = bj(k-1)
    enddo

    call atom_psd_krk_psi(6, bkrk, .TRUE., jrc, r, lo, rpsi, mxdnr)

    call atom_psd_psi_charge_nrl(jrc, drdi, rpsi, cdps, jrc)

!   Calculate new delta(with gamma fixed), uses false position

    fdnew = log(cdrc/cdps) - 2*delta
    if (abs(fdnew) < errmin) then
      v0pp = 8*((2*lo+5)*bj(1) + gama*gama)
      iflag = 0

      exit

    endif
    if (j == 1) then
      ddelta = -ONE/2
    else
      ddelta = - fdnew*ddelta / (fdnew-fdold)
    endif
    delta = delta + ddelta

    call atom_psd_tm2_ajbj_gam(bjin, rcpn, delta, gama, aj, bj)

    fdold = fdnew
  enddo

  if (iflag /= 0) then
    write(6,*)
    write(6,*)
    write(6,'("  error in psd_tm2_gamfn - delta not found")')

  endif

  deallocate(aj)

  return

end subroutine atom_psd_tm2_gamfn
