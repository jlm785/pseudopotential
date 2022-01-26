!>  generates a pseudopotential using the
!>  scheme of J. L. Martins, C. L. Reis and J. Pacheco
!>
!>  \author       Carlos Loia Reis, Jose Luis Martins
!>  \version      6.0.7
!>  \date         1 November 2021, 19 November 2021.
!>  \copyright    GNU Public License v2

subroutine atom_psd_mrpp(ispp, nr, r, drdi, d2rodr,                      &
      no, lo, iso, rc,                                                   &
      vionic, vhxc, ev1, ev2,                                            &
      rpsi_ae, br_ae, rpsi_ps, vscr,                                     &
      iowrite, mxdnr)


 implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter          :: MRPP = 8

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  character(len=1), intent(in)      ::  ispp                             !<  spin polarization  ' ', 's', 'r'

  integer, intent(in)               ::  nr                               !<  number of the number of radial points
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i
  real(REAL64), intent(in)          ::  d2rodr(mxdnr)                    !<  (d^2 r(i) / d i^2) / (d r / d i)

  integer, intent(in)               ::  no                               !<  principal quantum number n
  integer, intent(in)               ::  lo                               !<  angular quantum number l
  integer, intent(in)               ::  iso                              !<  spin quantum number

  real(REAL64), intent(in)          ::  vionic(mxdnr)                    !<  r*ionic potential in Rydberg

  real(REAL64), intent(in)          ::  vhxc(mxdnr)                      !<  screening potential (Ry)

  real(REAL64), intent(in)          ::  rpsi_ae(mxdnr,2)                 !<  r*all-electron-wave-function (major component if relativistic)
  real(REAL64), intent(in)          ::  br_ae(mxdnr,2)                   !<  d rpsi_ae / dr or minor component

! output

  real(REAL64), intent(out)         ::  rpsi_ps(mxdnr,2)                 !<  r*pseudo-wave-function

  real(REAL64), intent(out)         ::  vscr(mxdnr)                      !<  screened-pseudo-potential in Rydberg

! input and output

  real(REAL64), intent(inout)       ::  ev1                              !<  orbital energy of lower orbital
  real(REAL64), intent(inout)       ::  ev2                              !<  orbital energy of second orbital

  real(REAL64), intent(inout)       ::  rc                               !<  core radius r_c(l)

! local allocatable arrays

  real(REAL64), allocatable         ::  ar(:)
  real(REAL64), allocatable         ::  br(:)
  real(REAL64), allocatable         ::  ar2(:)
  real(REAL64), allocatable         ::  br2(:)

  real(REAL64), allocatable         ::  vj(:)

  real(REAL64), allocatable         ::  rz_all(:)
  real(REAL64), allocatable         ::  rx_all(:)
  real(REAL64), allocatable         ::  ax_all(:)
  real(REAL64), allocatable         ::  bx_all(:)

! local variables

  real(REAL64)         ::  bkrk(0:MRPP)                                  !  coefficient of Kerker polynomial
  real(REAL64)         ::  polydrc(0:MRPP-4)                             !  n-th derivatives of Kerker polynomial ar rc

  real(REAL64)         ::  rzero, rextr
  integer              ::  nextr, nzero
  integer              ::  nn

  integer              ::  jrc, ist
  real(REAL64)         ::  cdrc
  real(REAL64)         ::  cdrc2
  real(REAL64)         ::  ar2p                                          !  derivative of second function at rc

  integer              ::  iflag

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  TOL = 1.0E-12_REAL64
  character(len=1), parameter   ::  IL(0:6) = (/'s','p','d','f','g','h','i'/)

! counters

  integer        ::  j


  allocate(ar(nr))
  allocate(br(nr))
  allocate(ar2(nr))
  allocate(br2(nr))

  allocate(vj(nr))

  do j = 1,MRPP
    bkrk(j) = ZERO
  enddo
  do j = 1,MRPP-4
    polydrc(j) = ZERO
  enddo

  do j = 1,nr
    ar(j) = rpsi_ae(j,1)
    br(j) = br_ae(j,1)
    ar2(j) = rpsi_ae(j,2)
    br2(j) = br_ae(j,2)
  enddo

! Find last zero and extremum

  nn = no-lo + 2

  allocate(rz_all(nn))
  allocate(rx_all(nn))
  allocate(ax_all(nn))
  allocate(bx_all(nn))

  call atom_atm_zero_extr(ispp, ar, br, nr, r,                           &
      no, lo, iso ,vionic, vhxc, ev1,                                    &
      nzero, nextr, rz_all, rx_all, ax_all, bx_all,                      &
      mxdnr)

  rzero = rz_all(nzero)
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
    call atom_psd_psi_charge_rel(jrc, drdi, ar2, br2, cdrc2, nr)
  else
    call atom_psd_psi_charge_nrl(jrc, drdi, ar, cdrc, nr)
    call atom_psd_psi_charge_nrl(jrc, drdi, ar2, cdrc2, nr)
  endif

! Find the values for wave(arc), d(wave)/dr(arp), potential(vrc),
! d(potential)/dr(vrp), and d2(potential)/dr2(vrpp)
! and converts do d^n poly / d r^n (r_c)

  call atom_psd_mrpp_atrc(ispp, jrc, r, ar(jrc), br(jrc), ar2(jrc), br2(jrc),  &
      lo, iso, ev1, ev2, vionic, vhxc,                                         &
      polydrc, ar2p,                                                           &
      mxdnr)

! Solves the set of equations for the Troullier-Martins pseudopotential
! as the input of mrpp

    call atom_psd_tm2_eqsolve_rci(r, drdi, jrc, lo,                      &
      polydrc, cdrc, bkrk,                                               &
      mxdnr)

!   Solves the set of equations for MRPP

    call atom_psd_mrpp_eqsolve_rci(r, drdi, d2rodr, jrc, lo,             &
      polydrc, cdrc, cdrc2, ar2(jrc), ar2p, ev1, ev2, bkrk,              &
      iowrite, mxdnr)

! Construct pseudo-wave-function

  call atom_psd_krk_psi(MRPP, bkrk, .TRUE., jrc, r, lo, ar, mxdnr)

! Invert Schrödinger equation to find new potential.

  call atom_psd_krk_v(MRPP, bkrk, .TRUE., jrc, r, lo, ev1, vj, mxdnr)

  do j = 1,jrc
    vscr(j) = vj(j)
  enddo
  if(jrc < nr) then
    do j = jrc+1,nr
      vscr(j) = vionic(j) / r(j) + vhxc(j)
    enddo
  endif

  do j = 1,nr
    rpsi_ps(j,1) = ar(j)
  enddo

! finds the second wave-function

  do j = 2,nr
    vj(j) = vscr(j) + (lo*(lo+1))/(r(j)*r(j))
  enddo

  call atom_atm_difnrl_bound(nr, r, drdi, d2rodr, vj, ar, br,            &
      lo+2, lo, ev2, iflag, TOL,                                         &
      iowrite, mxdnr)

  do j = 1,nr
    rpsi_ps(j,2) = ar(j)
  enddo

  if(lo < 7) then
    write(iowrite,'(1x,i1,a1,f6.1,14f12.6)') lo+1, IL(lo), 0.5*ISO,      &
         ev1, rc, cdrc, (bkrk(j),j=0,MRPP)
  else
    write(iowrite,'(1x,i1,i1,f6.1,14f12.6)') lo+1, lo, 0.5*ISO,          &
         ev1, rc, cdrc, (bkrk(j),j=0,MRPP)
  endif

  deallocate(ar)
  deallocate(br)
  deallocate(ar2)
  deallocate(br2)

  deallocate(vj)

  return

end subroutine atom_psd_mrpp
