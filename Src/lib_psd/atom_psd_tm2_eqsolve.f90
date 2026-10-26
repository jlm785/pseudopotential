!>  Solves the equations of the
!>  improved scheme of N. Troullier and J. L. Martins
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.0.6
!>  \date         April 1990, 29 October 2021, 26 October 2026.
!>  \copyright    GNU Public License v2

subroutine atom_psd_tm2_eqsolve(r, drdi, jrc, lo,                        &
     polydrc, cdrc, bkrk,                                                &
     mxdnr)

! extracted from the atom_psd_tm2 code. 29 October 2021. JLM
! iflag initialized, overflow and zero charge checks, hint on failure. 26 October 2026. claude
! roots inside the initial interval (atom_psd_tm2_inbrac). 26 October 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension  of the number of radial points

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  integer, intent(in)               ::  jrc                              !<  r(jrc) = r_cut

  real(REAL64), intent(in)          ::  cdrc                             !<  int_0^r_cut rho(r) dr
  real(REAL64), intent(in)          ::  polydrc(0:4)                     !<  n-th derivatives of Kerker polynomial ar rc

  integer, intent(in)               ::  lo                               !<  angular quantum number l

! output

  real(REAL64), intent(out)         ::  bkrk(0:6)                        !<  coefficient of polynomial

! local allocated variables

  real(REAL64), allocatable         ::  rpsi(:)                          !  wave-function up to jrc

! local variables

  real(REAL64)         ::  bjin(5)                                       !  derivatives of the Kerker polynomial
  real(REAL64)         ::  aj(5,5)                                       !  for solving linear system
  real(REAL64)         ::  bj(5)                                         !  polyr coefficients
  real(REAL64)         ::  gama                                          !  value of gamma
  real(REAL64)         ::  delta                                         !  value of delta
  real(REAL64)         ::  rcpn(12)                                      !  rc^n
  real(REAL64)         ::  cdps                                          !  int_0^r_cut rho_psd(r) dr

  real(REAL64)         ::  fdnew, fdold, ddelta
  integer              ::  iflag
  real(REAL64)         ::  x1, x2
  real(REAL64)         ::  polymax
  logical              ::  lfound

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  real(REAL64), parameter    ::  SMALL = 1.0E-12_REAL64
  real(REAL64), parameter    ::  PMAX = 100.0_REAL64                     !  maximum exponent of the trial wave-function (about 3 in solutions)

! counter

  integer      ::  j, k


! Set up matrix without the d2(potential(0)/dr2=0 condition
! to find an intial guess for gamma.

  rcpn(1) = r(jrc)
  do j = 2,12
    rcpn(j) = rcpn(j-1)*r(jrc)
  enddo

  do j = 1,5
    bjin(j) = polydrc(j-1)
  enddo

  delta = ZERO
  iflag = 1

! could allocate only up to jrc...

  allocate(rpsi(mxdnr))

  call atom_psd_tm2_ajbj(bjin, rcpn, ZERO, aj, bj)

! Start iteration loop to find delta, uses false postion.

  do j = 1,50

!   Generate pseudo wavefunction-note missing factor exp(delta).

    gama = bj(1)
    bj(1) = bj(2)
    bj(2) = bj(3)
    bj(3) = bj(4)
    bj(4) = bj(5)
    bj(5) = ZERO


    bkrk(0) = ZERO
    bkrk(1) = gama
    do k = 2,6
      bkrk(k) = bj(k-1)
    enddo

    call atom_psd_krk_polymax(6, bkrk, .TRUE., jrc, r, polymax, mxdnr)

    if(polymax > PMAX) exit

    call atom_psd_krk_psi(6, bkrk, .TRUE., jrc, r, lo, rpsi, mxdnr)

!   Integrate pseudo charge density from r = 0 to rc.

    call atom_psd_psi_charge_nrl(jrc, drdi, rpsi, cdps, mxdnr)

    if(cdps <= ZERO) exit

!   Calculate new delta

    fdnew = log(cdrc/cdps) - 2*delta
    if (abs(fdnew) < SMALL) then
      iflag = 0

      exit

    endif
    if (j == 1) then
      ddelta = -0.5
    else
      ddelta = - fdnew * ddelta / (fdnew-fdold)
    endif
    delta = delta + ddelta

    call atom_psd_tm2_ajbj(bjin, rcpn, delta, aj, bj)

    fdold = fdnew
  enddo

! End iteration loop for delta.

  if (iflag /= 0) then

    write(6,*)
    write(6,*)
    write(6,'(" error in atom_psd_tm2_eqsolve - nonconvergence")')
    write(6,'("  in finding starting delta for l = ",i3)') lo
    write(6,'("  try a different core radius for l = ",i3)') lo

    STOP

  endif

! Bracket the correct gamma, use gamma and -gamma
! from above as intial brackets, expands brackets
! until a root is found..
! If v''(0) has the same sign at both ends and the interval
! contains roots, uses the one closest to the initial gamma.

  x1 = gama
  x2 =-gama

  call atom_psd_tm2_inbrac(x1,x2, lfound,                                &
      jrc, r, drdi, bjin, rcpn, cdrc, lo, rpsi, delta, bj,               &
      mxdnr)

  if(lfound) then
    iflag = 0
  else
    call atom_psd_tm2_zrbac(x1,x2, iflag,                                &
        jrc, r, drdi, bjin, rcpn, cdrc, lo, rpsi, delta, bj,             &
        mxdnr)
  endif

  if(iflag /= 0) then

    write(6,'("   stopped in atom_psd_tm2_eqsolve. iflag = ",i3)') iflag
    write(6,'("   no Troullier-Martins solution found, try a different core radius for l = ",i3)') lo

    STOP

  endif

! Iteration loop to find correct gamma, uses
! bisection to find gamma.

  call atom_psd_tm2_rtbis(x1,x2, iflag,                                  &
   jrc, r, drdi, bjin, rcpn, cdrc, lo, rpsi, gama, delta, bj,            &
   mxdnr)

  if(iflag /= 0) then
    write(6,'("   stopped in atom_psd_tm2_eqsolve. iflag = ",i3)') iflag
    write(6,'("   no Troullier-Martins solution found, try a different core radius for l = ",i3)') lo

    STOP

  endif

  bkrk(0) = delta
  bkrk(1) = gama
  do j = 2,6
    bkrk(j) = bj(j-1)
  enddo

  deallocate(rpsi)

  return

end subroutine atom_psd_tm2_eqsolve
