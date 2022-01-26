!>  Calculates the values of the derivatives of the Kerker
!>  polynomial at the core radius from the
!>  wave-function and screened all-electron potential
!>  at r(jrc) = r_cut, and their derivatives.
!>
!>  \author       Sverre Froyen,  Norm Troullier, Jose Luis Martins
!>  \version      6.0.6
!>  \date         1980s and 1990s, 30 June 2021, 30 October 2021
!>  \copyright    GNU Public License v2

subroutine atom_psd_mrpp_atrc(ispp, jrc, r, arjrc, brjrc, ar2jrc, br2jrc,   &
        lo, iso, ev1, ev2, vionic, vhxc,                                    &
        polydrc, ar2p,                                                      &
        mxdnr)

! extracted from pseudo2. 30 June 2021. JLM
! merged with atom_psd_tm2_bjin_rcpn. 30 October 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension  of the number of radial points

  character(len=1), intent(in)      ::  ispp                             !<  spin polarization  ' ', 's', 'r'

  integer, intent(in)               ::  jrc                              !<  r(jrc) = r_cut

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid

  real(REAL64), intent(in)          ::  arjrc, brjrc                     !<  wave-function ar(jrc), br(jrc)
  real(REAL64), intent(in)          ::  ar2jrc, br2jrc                   !<  second wave-function ar2(jrc), br2(jrc)

  integer, intent(in)               ::  lo                               !<  angular quantum number l
  integer, intent(in)               ::  iso                              !<  2*spin or 2*(j-l)

  real(REAL64), intent(in)          ::  ev1                              !<  orbital energy
  real(REAL64), intent(in)          ::  ev2                              !<  second orbital energy

  real(REAL64), intent(in)          ::  vionic(mxdnr)                    !<  r*ionic potential in Rydberg

  real(REAL64), intent(in)          ::  vhxc(mxdnr)                      !<  effective potential in Rydberg (down or total)

! output

  real(REAL64), intent(out)         ::  polydrc(0:4)                     !<  n-th derivatives of Kerker polynomial at r_c
  real(REAL64), intent(out)         ::  ar2p                             !<  first derivatives of second function

! allocatable local variables

  real(REAL64), allocatable         ::  aa(:), rr(:), coe(:), cerror(:)
  real(REAL64), allocatable         ::  rcpn(:)

! local variables

  real(REAL64)      ::  arp
  integer           ::  ka
  integer           ::  lp

  real(REAL64)      ::  arc, brc                                         !  psi(r_cut) and psi'(r_cut)/psi(r_cut)
  real(REAL64)      ::  vrc                                              !  screened all-electron potential v(r_cut) at r_cut
  real(REAL64)      ::  vap, vapp                                        !  v'(r_cut) and v''(r_cut)

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  AI = 137.0359991_REAL64

! counter

  integer      ::  j



  if(jrc < 4) stop

  allocate(aa(7), rr(7), coe(7), cerror(7))

  ka = lo + 1
  if (iso == -1 .and. lo /= 0) ka = -lo

  arc = arjrc

  if (ispp == 'r') then
    arp = ka*arjrc/r(jrc) + (brjrc / (2*AI)) *                           &
                      (ev1 - vionic(jrc)/r(jrc) - vhxc(jrc) + 4*AI*AI)
  else
    arp = brjrc
  endif
  brc = arp / arc

  if (ispp == 'r') then
    ar2p = ka*ar2jrc/r(jrc) + (br2jrc / (2*AI)) *                        &
                      (ev2 - vionic(jrc)/r(jrc) - vhxc(jrc) + 4*AI*AI)
  else
    ar2p = br2jrc
  endif

  aa(1) = vionic(jrc-3)/r(jrc-3) + vhxc(jrc-3)
  aa(2) = vionic(jrc-2)/r(jrc-2) + vhxc(jrc-2)
  aa(3) = vionic(jrc-1)/r(jrc-1) + vhxc(jrc-1)
  aa(4) = vionic(jrc  )/r(jrc  ) + vhxc(jrc  )
  aa(5) = vionic(jrc+1)/r(jrc+1) + vhxc(jrc+1)
  aa(6) = vionic(jrc+2)/r(jrc+2) + vhxc(jrc+2)
  aa(7) = vionic(jrc+3)/r(jrc+3) + vhxc(jrc+3)
  vrc = aa(4)

  rr(1) = r(jrc-3)-r(jrc)
  rr(2) = r(jrc-2)-r(jrc)
  rr(3) = r(jrc-1)-r(jrc)
  rr(4) = ZERO
  rr(5) = r(jrc+1)-r(jrc)
  rr(6) = r(jrc+2)-r(jrc)
  rr(7) = r(jrc+3)-r(jrc)


  call poly_interp(coe, cerror, rr, aa, 6, 2)

  vap   = coe(2)
  vapp  = coe(3)

  deallocate(aa, rr, coe, cerror)

  allocate(rcpn(max(4,lo+1)))

  rcpn(1) = r(jrc)
  do j = 2,max(4,lo+1)
    rcpn(j) = rcpn(j-1)*r(jrc)
  enddo

  lp = lo+1

  polydrc(0) = log(arc/rcpn(lp))
  polydrc(1) = brc - lp/rcpn(1)
  polydrc(2) = vrc - ev1 - 2*lp*polydrc(1)/rcpn(1)                       &
                        - polydrc(1)*polydrc(1)
  polydrc(3) = vap + 2*lp*polydrc(1)/rcpn(2) - 2*lp*polydrc(2)/rcpn(1)   &
                   - 2*polydrc(1)*polydrc(2)
  polydrc(4) = vapp - 4*lp*polydrc(1)/rcpn(3) + 4*lp*polydrc(2)/rcpn(2)  &
                    - 2*lp*polydrc(3)/rcpn(1) - 2*polydrc(2)*polydrc(2)  &
                    - 2*polydrc(1)*polydrc(3)

  deallocate(rcpn)

  return

end subroutine atom_psd_mrpp_atrc
