!>  Integrates the radial Schrodinger equation using a predictor corrrector method,
!>  needs 5 points to start.
!>
!>  The wavefunction is calculated between nbegin and nend
!>  For outward intergration (nbegin < nend) rpsi(nbegin-5),...,rpsi(nbegin-1) is needed as input
!>  For inward integration (nbegin >= nend) rpsi(nbegin+1),...,rpsi(nbegin+5) is needed as input
!>  For both cases the corresponding drpsidi is also needed.
!>
!>  To deal with the case v(1) -> infty, vrpsi0 = lim r-> 0 v(r) rpsi(r) is also needed for nbegin = 6 < nend
!>
!>  \author       Sverre Froyen,  Norm Troullier, Jose Luis Martins
!>  \version      6.013
!>  \date         1980s, 22 June 2021
!>  \copyright    GNU Public License v2

subroutine atom_atm_integ_nr(nbegin, nend, drdisq, d2rodr,               &
      v, ev, rpsi, drpsidi, vrpsi0,                                      &
      iowrite, mxdnr)

! written May 2005 based on Sverre Froyen
! code as modified by Norm Troullier
! modified (f90) in 29 May 2012  JLM
! merged inward and and outward codes and
! added checks. 16 July 2019.  JLM

! jlm version 6.01

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input:

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points

  integer, intent(in)               ::  nbegin                           !<  starting point for integration
  integer, intent(in)               ::  nend                             !<  ending point for integration

  real(REAL64), intent(in)          ::  drdisq(mxdnr)                    !<  (d r(i) / d i)**2
  real(REAL64), intent(in)          ::  d2rodr(mxdnr)                    !<  (d^2 r(i) / d i^2) / (d r / d i)

  real(REAL64), intent(in)          ::  v(mxdnr)                         !<  effective potential (RYDBERG)
  real(REAL64), intent(in)          ::  ev                               !<  orbital energy
  real(REAL64), intent(in)          ::  vrpsi0                           !<  special correction for the case r(1)=0, v(1) -> infty. l=0 -2z, l=1 2,l>1 0.

!input and output

  real(REAL64), intent(inout)       ::  rpsi(mxdnr)                      !<  radial wave-function u = rR(r)  (integral rpsi**2 = 1)
  real(REAL64), intent(inout)       ::  drpsidi(mxdnr)                   !<  d rpsi / d i (derivative with respect to index i of r(i).

! local variables

  integer                    ::  idir                                    !  indicate if integration is inwards or outwards

  real(REAL64)               :: faj1,faj2,faj3,faj4,faj5                 !  5 values of rpsi
  real(REAL64)               :: fbj1,fbj2,fbj3,fbj4,fbj5                 !  5 values of drpsidi
  real(REAL64)               :: vev,fb0,fb1
  real(REAL64)               :: arp,brp,arc,brc                          !  predicted and corrected rpsi and drpsidi

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

  integer                  :: j, nn


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
  faj5 = drpsidi(nn)
  if(nn == 1) then
    fbj5 = d2rodr(nn)*drpsidi(nn) + drdisq(nn)*vrpsi0
  else
    fbj5 = d2rodr(nn)*drpsidi(nn) + drdisq(nn)*(v(nn)-ev)*rpsi(nn)
  endif
  nn = nbegin-idir*4
  faj4 = drpsidi(nn)
  fbj4 = d2rodr(nn)*drpsidi(nn) + drdisq(nn)*(v(nn)-ev)*rpsi(nn)
  nn = nbegin-idir*3
  faj3 = drpsidi(nn)
  fbj3 = d2rodr(nn)*drpsidi(nn) + drdisq(nn)*(v(nn)-ev)*rpsi(nn)
  nn = nbegin-idir*2
  faj2 = drpsidi(nn)
  fbj2 = d2rodr(nn)*drpsidi(nn) + drdisq(nn)*(v(nn)-ev)*rpsi(nn)
  nn = nbegin-idir
  faj1 = drpsidi(nn)
  fbj1 = d2rodr(nn)*drpsidi(nn) + drdisq(nn)*(v(nn)-ev)*rpsi(nn)

! intergration loop

  do j = nbegin,nend,idir

!   predictor (Adams-Bashforth)

    vev=v(j)-ev
    arp = rpsi(j-idir) + idir*(ABC1*faj1+ABC2*faj2+ABC3*faj3+           &
            ABC4*faj4+ABC5*faj5)
    brp = drpsidi(j-idir) + idir*(ABC1*fbj1+ABC2*fbj2+ABC3*fbj3+        &
            ABC4*fbj4+ABC5*fbj5)
    fb1 = d2rodr(j)*brp + drdisq(j)*vev*arp

!   corrector (Adams-Moulton)

    arc = rpsi(j-idir) + idir*(AMC0*brp+AMC1*faj1+AMC2*faj2+            &
            AMC3*faj3+AMC4*faj4)
    brc = drpsidi(j-idir) + idir*(AMC0*fb1+AMC1*fbj1+AMC2*fbj2+         &
            AMC3*fbj3+AMC4*fbj4)
    fb0 = d2rodr(j)*brc + drdisq(j)*vev*arc

!   error reduction step

    rpsi(j) = arc + idir*AMC0*(brc-brp)
    drpsidi(j) = brc + idir*AMC0*(fb0-fb1)

!   shift

    faj5 = faj4
    fbj5 = fbj4
    faj4 = faj3
    fbj4 = fbj3
    faj3 = faj2
    fbj3 = fbj2
    faj2 = faj1
    fbj2 = fbj1
    faj1 = drpsidi(j)
    fbj1 = d2rodr(j)*drpsidi(j) + drdisq(j)*vev*rpsi(j)

  enddo

  return

end subroutine atom_atm_integ_nr
