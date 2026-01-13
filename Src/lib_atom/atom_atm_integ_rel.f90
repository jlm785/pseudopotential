!>  Integrates the radial "Dirac" equation using a predictor corrrector method,
!>  needs 5 points to start.
!>
!>  The wavefunction is calculated between nbegin and nend
!>
!>  \author       Sverre Froyen,  Norm Troullier, Jose Luis Martins
!>  \version      6.013
!>  \date         1980s, 22 June 2021
!>  \copyright    GNU Public License v2

subroutine atom_atm_integ_rel(nbegin, nend, r, drdi,                     &
        v, ev, ka, ar, br, fanctp,                                       &
        iowrite, mxdnr)

! written June 2021 based on Sverre Froyen
! code as modified by Norm Troullier

! jlm version 6.013

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  nbegin                           !<  starting point for integration
  integer, intent(in)               ::  nend                             !<  ending point for integration

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  real(REAL64), intent(in)          ::  v(mxdnr)                         !<  effective potential (RYDBERG)
  real(REAL64), intent(in)          ::  ev                               !<  orbital energy

  integer, intent(in)               ::  ka

  real(REAL64), intent(inout)       ::  ar(mxdnr)                        !<  major component radial wave-function u = rR(r)
  real(REAL64), intent(inout)       ::  br(mxdnr)                        !<  minor component (integral ar**2 + br**2 = 1)
  real(REAL64), intent(out)         ::  fanctp

! local variables

  integer                    ::  idir                                    !  indicate if integration is inwards or outwards

  real(REAL64)       ::  evv, evvai2
  real(REAL64)       ::  arp, brp
  real(REAL64)       ::  arc, brc
  real(REAL64)       ::  faj, fbj

  real(REAL64)       ::  faj0,faj1,faj2,faj3,faj4,faj5                   !  6 values of arprime
  real(REAL64)       ::  fbj0,fbj1,fbj2,fbj3,fbj4,fbj5                   !  6 values of brprime
! it would be more elegant to have a small array, but it would not be thread safe.

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64

  real(REAL64), parameter    ::  AI = 137.035999177_REAL64
  real(REAL64), parameter    ::  AI2 = 4*AI*AI

! Adams-Bashforth
  real(REAL64), parameter    :: ABC1 = 190.1_REAL64/72.0_REAL64
  real(REAL64), parameter    :: ABC2 =-138.7_REAL64/36.0_REAL64
  real(REAL64), parameter    :: ABC3 =  10.9_REAL64/3.0_REAL64
  real(REAL64), parameter    :: ABC4 = -63.7_REAL64/36.0_REAL64
  real(REAL64), parameter    :: ABC5 =  25.1_REAL64/72.0_REAL64
! Adams-Moulton
  real(REAL64), parameter    :: AMC0 =  25.1_REAL64/72.0_REAL64
  real(REAL64), parameter    :: AMC1 =  32.3_REAL64/36.0_REAL64
  real(REAL64), parameter    :: AMC2 =  -1.1_REAL64/3.0_REAL64
  real(REAL64), parameter    :: AMC3 =   5.3_REAL64/36.0_REAL64
  real(REAL64), parameter    :: AMC4 =  -1.9_REAL64/72.0_REAL64

! counter

  integer                  ::  j, nn


! sanity checks

  if(nbegin < nend) then
    idir = 1
  else
    idir =-1
  endif

  if(idir == 1) then
    if(nbegin -5 < 1) then
      write(iowrite,'("  atm_integ_nr, nbegin = ",i8," is < 6")') nbegin

      stop

    endif
  else
    if(nbegin+5 > mxdnr) then
      write(iowrite,'("  atm_integ_nr, nbegin = ",i8," is >",i8)') nbegin,mxdnr-5

      stop

    endif
  endif

! startup

  nn = nbegin-idir*5
  if(nn == 1) then
    faj5 = ZERO
    fbj5 = ZERO
  else
    faj5 =  drdi(nn)*ar(nn)*ka / r(nn) + (ev-v(nn)+AI2)*br(nn)*drdi(nn) / (2*AI)
    fbj5 = -drdi(nn)*br(nn)*ka / r(nn) - (ev-v(nn))*ar(nn)*drdi(nn) / (2*AI)
  endif

  nn = nbegin-idir*4
  faj4 =  drdi(nn)*ar(nn)*ka / r(nn) + (ev-v(nn)+AI2)*br(nn)*drdi(nn) / (2*AI)
  fbj4 = -drdi(nn)*br(nn)*ka / r(nn) - (ev-v(nn))*ar(nn)*drdi(nn) / (2*AI)
  nn = nbegin-idir*3
  faj3 =  drdi(nn)*ar(nn)*ka / r(nn) + (ev-v(nn)+AI2)*br(nn)*drdi(nn) / (2*AI)
  fbj3 = -drdi(nn)*br(nn)*ka / r(nn) - (ev-v(nn))*ar(nn)*drdi(nn) / (2*AI)
  nn = nbegin-idir*2
  faj2 =  drdi(nn)*ar(nn)*ka / r(nn) + (ev-v(nn)+AI2)*br(nn)*drdi(nn) / (2*AI)
  fbj2 = -drdi(nn)*br(nn)*ka / r(nn) - (ev-v(nn))*ar(nn)*drdi(nn) / (2*AI)
  nn = nbegin-idir*1
  faj1 =  drdi(nn)*ar(nn)*ka / r(nn) + (ev-v(nn)+AI2)*br(nn)*drdi(nn) / (2*AI)
  fbj1 = -drdi(nn)*br(nn)*ka / r(nn) - (ev-v(nn))*ar(nn)*drdi(nn) / (2*AI)

! Intergration loop.

  do j = nbegin,nend,idir

!   Predictor (Adams-Bashforth).

    evvai2 = ev - v(j) + AI2
    evv = ev - v(j)
    arp = ar(j-idir) + idir*(ABC1*faj1 + ABC2*faj2 + ABC3*faj3 + ABC4*faj4 + ABC5*faj5)
    brp = br(j-idir) + idir*(ABC1*fbj1 + ABC2*fbj2 + ABC3*fbj3 + ABC4*fbj4 + ABC5*fbj5)
    faj0 =  drdi(j)*arp*ka / r(j) + evvai2*brp*drdi(j) / (2*AI)
    fbj0 = -drdi(j)*brp*ka / r(j) - evv*arp*drdi(j) / (2*AI)

!   Corrector (Adams-Moulton).

    arc = ar(j-idir) + idir*(AMC0*faj0 + AMC1*faj1 + AMC2*faj2 + AMC3*faj3 + AMC4*faj4)
    brc = br(j-idir) + idir*(AMC0*fbj0 + AMC1*fbj1 + AMC2*fbj2 + AMC3*fbj3 + AMC4*fbj4)
    faj =  drdi(j)*arc*ka / r(j) + evvai2*brc*drdi(j) / (2*AI)
    fbj = -drdi(j)*brc*ka / r(j) - evv*arc*drdi(j) / (2*AI)

!   Error reduction step.

    ar(j) = arc + idir*AMC0*(faj-faj0)
    br(j) = brc + idir*AMC0*(fbj-fbj0)

!   shift

    faj5 = faj4
    fbj5 = fbj4
    faj4 = faj3
    fbj4 = fbj3
    faj3 = faj2
    fbj3 = fbj2
    faj2 = faj1
    fbj2 = fbj1
    faj1 =  drdi(j)*ar(j)*ka / r(j) + evvai2*br(j)*drdi(j) / (2*AI)
    fbj1 = -drdi(j)*br(j)*ka / r(j) - evv*ar(j)*drdi(j) / (2*AI)

  enddo

  fanctp = faj1

  return

end subroutine atom_atm_integ_rel
