!>  Identifies the zeroes and extrmes of a wave-function
!>
!>  \author       Sverre Froyen,  Norm Troullier, Jose Luis Martins
!>  \version      6.0.6
!>  \date         1980s, 22 June 2021, 31 October 2021.
!>  \copyright    GNU Public License v2

subroutine atom_atm_zero_extr(ispp, ar, br, nr, r,                       &
      no, lo, iso ,vionic, vhxc, ev,                                     &
      nzero, nextr, rzero, rextr, aextr, bextr,                          &
      mxdnr)

! extracted from atom_atm_orban to consolidate repeated code lines. 31 October 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points

  character(len=1), intent(in)      ::  ispp                             !<  spin polarization  ' ', 's', 'r'

  real(REAL64), intent(in)          ::  ar(mxdnr)                        !<  radial wave-function u = rR(r)  (integral ar**2 = 1) or major component radial wave-function u = rR(r)
  real(REAL64), intent(in)          ::  br(mxdnr)                        !<  d ar / d r    or  minor component (integral ar**2 + br**2 = 1)

  integer, intent(in)               ::  nr                               !<  dimension of the number of radial points

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points

  integer, intent(in)               ::  no                               !<  principal quantum number n
  integer, intent(in)               ::  lo                               !<  angular quantum number l
  integer, intent(in)               ::  iso                              !<  2*spin or 2*(j-l)

  real(REAL64), intent(in)          ::  vionic(mxdnr)                    !<  r*ionic potential for orbital (Ry)

  real(REAL64), intent(in)          ::  vhxc(mxdnr)                      !<  screening potential (Ry)

  real(REAL64), intent(in)          ::  ev                               !<  orbital energy

! output

  integer, intent(out)              ::  nzero                            !<  number of zeroes
  integer, intent(out)              ::  nextr                            !<  number of extrema

  real(REAL64), intent(out)         ::  rzero(no-lo + 2)                 !<  non trivial node of ar
  real(REAL64), intent(out)         ::  rextr(no-lo + 2)                 !<  maximum of wavefunction
  real(REAL64), intent(out)         ::  aextr(no-lo + 2)                 !<  wave-function or major component at rextr
  real(REAL64), intent(out)         ::  bextr(no-lo + 2)                 !<  monor component at rextr

! local

  integer                ::  lp, ka

  real(REAL64)           ::  arp, arpm

! constants

  real(REAL64), parameter   ::  ZERO = 0.0_REAL64
  real(REAL64), parameter   ::  AI = 2*137.0359991_REAL64

! counters

  integer     ::  i


  ka = lo+1
  lp = ka
  if (iso == -1 .and. lo /= 0) ka = -lo

! compute zeroes and extrema

  nzero = 0
  nextr = 0
  rzero(1) = ZERO
  arp = br(2)
  if (ispp == 'r') then
    arp = ka*ar(2)/r(2) + (ev - vionic(2)/r(2)                           &
         - vhxc(2) + AI*AI) * br(2) / AI
  endif

  do i = 3,nr

    if (nextr >= no-lo) exit

    if (ar(i)*ar(i-1) <= ZERO) then
!     zero
      nzero = nzero + 1
      rzero(nzero) = (ar(i)*r(i-1)-ar(i-1)*r(i)) / (ar(i)-ar(i-1))
    endif

    arpm = arp
    arp = br(i)
    if (ispp == 'r') then
        arp = ka*ar(i)/r(i) + (ev - vionic(i)/r(i) - vhxc(i) + AI*AI)    &
                  *br(i) / AI
    endif

    if (arp*arpm <= ZERO) then
!     extremum
      nextr = nextr + 1
      rextr(nextr) = (arp*r(i-1)-arpm*r(i)) / (arp-arpm)
      aextr(nextr) = (ar(i)+ar(i-1))/2                                   &
            - (arp**2+arpm**2) * (r(i)-r(i-1)) / (4*(arp-arpm))
      bextr(nextr) = br(i)
    endif

  enddo

  return

end subroutine atom_atm_zero_extr
