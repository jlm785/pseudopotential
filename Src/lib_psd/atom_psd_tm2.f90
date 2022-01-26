!>  generates a pseudopotential using the
!>  improved scheme of N. Troullier and J. L. Martins
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.0.6
!>  \date         April 1990, 30 June 2021, 1 November 2021.
!>  \copyright    GNU Public License v2

subroutine atom_psd_tm2(ispp, nr, r, drdi, lrci,                         &
      no, lo, iso, rc,                                                   &
      vionic, vhxc, ev,                                                  &
      rpsi_ae, br_ae, rpsi_ps, vscr,                                     &
      iowrite, mxdnr)

! *************************************************************
! *                                                           *
! *     This routine was written by Norman J. Troullier Jr.   *
! *   April 1990, while at the U. of Minnesota.               *
! *                                                           *
! *   The general format of this routine is the same as the   *
! *   pseudo and pseudk routines.  Output/input is            *
! *   compatible.                                             *
! *                                                           *
! *************************************************************

! modified output for compatibility with other programs May 12, 2010. JLM
! modified loop to zero nops 12 April 2012. JLM
! major modifications, common code with other pseudos separated,  July 2021. JLM
! so->iso, vionic, vhxc, cdpsd, vscr. 15, 23 September 2021. JLM
! major rewriting. Only one  orbital at a time, 1 November 2021. JLM


 implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter          :: NTM2 = 6

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  character(len=1), intent(in)      ::  ispp                             !<  spin polarization  ' ', 's', 'r'

  integer, intent(in)               ::  nr                               !<  number of the number of radial points
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  logical, intent(in)               ::  lrci                             !<  uses rverse communication interface

  integer, intent(in)               ::  no                               !<  principal quantum number n
  integer, intent(in)               ::  lo                               !<  angular quantum number l
  integer, intent(in)               ::  iso                              !<  spin quantum number

  real(REAL64), intent(in)          ::  vionic(mxdnr)                    !<  r*ionic potential in Rydberg

  real(REAL64), intent(in)          ::  vhxc(mxdnr)                      !<  screening potential (Ry)

  real(REAL64), intent(in)          ::  rpsi_ae(mxdnr)                   !<  r*all-electron-wave-function (major component if relativistic)
  real(REAL64), intent(in)          ::  br_ae(mxdnr)                     !<  d rpsi_ae / dr or minor component

! output

  real(REAL64), intent(out)         ::  rpsi_ps(mxdnr)                   !<  r*pseudo-wave-function

  real(REAL64), intent(out)         ::  vscr(mxdnr)                      !<  screened-pseudo-potential in Rydberg

! input and output

  real(REAL64), intent(inout)       ::  ev                               !<  orbital energy

  real(REAL64), intent(inout)       ::  rc                               !<  core radius r_c(l)

! local allocatable arrays

  real(REAL64), allocatable         ::  ar(:)
  real(REAL64), allocatable         ::  br(:)

  real(REAL64), allocatable         ::  vj(:)

  real(REAL64), allocatable         ::  rz_all(:)
  real(REAL64), allocatable         ::  rx_all(:)
  real(REAL64), allocatable         ::  ax_all(:)
  real(REAL64), allocatable         ::  bx_all(:)

! local variables

  real(REAL64)         ::  bkrk(0:NTM2)                                !  coefficient of Kerker polynomial
  real(REAL64)         ::  polydrc(0:NTM2-2)                           !  n-th derivatives of Kerker polynomial ar rc

  real(REAL64)         ::  rzero, rextr
  integer              ::  nextr, nzero
  integer              ::  nn

  integer              ::  jrc, ist
  real(REAL64)         ::  cdrc

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  character(len=1), parameter   ::  IL(0:6) = (/'s','p','d','f','g','h','i'/)

! counters

  integer        ::  j


  allocate(ar(nr))
  allocate(br(nr))

  allocate(vj(nr))

  do j = 1,nr
    ar(j) = rpsi_ae(j)
    br(j) = br_ae(j)
  enddo

! Find last zero and extremum

  nn = no-lo + 2

  allocate(rz_all(nn))
  allocate(rx_all(nn))
  allocate(ax_all(nn))
  allocate(bx_all(nn))

  call atom_atm_zero_extr(ispp, ar, br, nr, r,                           &
      no, lo, iso ,vionic, vhxc, ev,                                     &
      nzero, nextr, rz_all, rx_all, ax_all, bx_all,                      &
      mxdnr)

  if(nzero == 0) then
    rzero = ZERO
  else
    rzero = rz_all(nzero)
  endif
  rextr = rx_all(nextr)

  deallocate(rz_all)
  deallocate(rx_all)
  deallocate(ax_all)
  deallocate(bx_all)

! Check rc if inside rzero,
! reset to .9 between rmax and rzero if inside
! if rc(lo(i)) is negative, rc(lo(i)) is percent of way
! betweeen rzero and rmax.

  if (rc > rzero) then
  elseif(rc >= ZERO) then
    rc = rzero + 0.9*(rextr-rzero)
  else
    rc = rzero - rc*(rextr-rzero)
  endif

! Find the index for odd grid point closest to rc.

  do j=1,nr
    jrc = j
    if (r(j) > rc) exit
  enddo

  jrc = jrc - 1
  rc = r(jrc)

  ist = 1
  if (ar(jrc) < ZERO) ist = -1

  do j = 1,nr
    ar(j) = ar(j)*ist
    br(j) = br(j)*ist
  enddo

! Find the integrated charge inside rc(1-charge outside).

  if (ispp == 'r') then
    call atom_psd_psi_charge_rel(jrc, drdi, ar, br, cdrc, nr)
  else
    call atom_psd_psi_charge_nrl(jrc, drdi, ar, cdrc, nr)
  endif

! Find the values for wave(arc), d(wave)/dr(arp), potential(vrc),
! d(potential)/dr(vrp), and d2(potential)/dr2(vrpp)
! and converts do d^n poly / d r^n (r_c)

  call atom_psd_tm2_atrc(ispp, jrc, r, ar(jrc), br(jrc),                 &
      lo, iso, ev, vionic, vhxc,                                         &
      polydrc,                                                           &
      mxdnr)

! Solves the set of equations for the Troullier-Martins pseudopotential

  if(lrci) then
    call atom_psd_tm2_eqsolve_rci(r, drdi, jrc, lo,                      &
      polydrc, cdrc, bkrk,                                               &
      mxdnr)
  else
    call atom_psd_tm2_eqsolve(r, drdi, jrc, lo,                          &
      polydrc, cdrc, bkrk,                                               &
      mxdnr)
  endif

! Construct pseudo-wave-function

  call atom_psd_krk_psi(NTM2, bkrk, .TRUE., jrc, r, lo, ar, mxdnr)

! Invert Schrödinger equation to find new potential.

  call atom_psd_krk_v(NTM2, bkrk, .TRUE., jrc, r, lo, ev, vj, mxdnr)

  do j = 1,jrc
    vscr(j) = vj(j)
  enddo
  if(jrc < nr) then
    do j = jrc+1,nr
      vscr(j) = vionic(j) / r(j) + vhxc(j)
    enddo
  endif

  do j = 1,nr
    rpsi_ps(j) = ar(j)
  enddo

  if(lo < 7) then
    write(iowrite,'(1x,i1,a1,f6.1,10f12.6)') lo+1, IL(lo), 0.5*ISO,      &
         ev, rc, cdrc, (bkrk(j),j=0,NTM2)
  else
    write(iowrite,'(1x,i1,i1,f6.1,10f12.6)') lo+1, lo, 0.5*ISO,          &
         ev, rc, cdrc, (bkrk(j),j=0,NTM2)
  endif

  deallocate(ar)
  deallocate(br)

  deallocate(vj)

  return

end subroutine atom_psd_tm2
