!>  unscreens the pseudopotential
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.0.8
!>  \date         1980s and 1990s, 30 June 2021, 25 May 2022.
!>  \copyright    GNU Public License v2


subroutine atom_psd_unscreen(ifcore, icorr, ispp, nr, r, drdi,           &
    norb, ncore, lo, iso, zo, znuc, zel,                                 &
    cdpsd, cdc, vscr, vpsd,                                              &
    cfac, rcfac, zratio, zion,                                           &
    iowrite, mxdnr, mxdorb)

! Writen 13 April 2018 from the old codes.  The old
! pseudopotential generation subroutines were broken in several subroutines
! and converted to f90.
! Modified (cleaning) 8 July 2021. JLM
! Initializes vpsd. 25 May 2022. JLM

!mmga  modifications from early Sverre code by Manuel Maria Gonzalez Alemany
!njtj  modifications from early Sverre code by Norm Troullier

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter                ::  lc = 4                           !  maximum valence angular momentum

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points
  integer, intent(in)               ::  mxdorb                           !<  dimension of the number of orbitals

  integer, intent(in)               ::  iowrite                          !<  default tape for writing


  integer, intent(in)               ::  ifcore                           !<  0: no partial core correction; 1 : xc partial core correction; 2 : not recommended.
  character(len=2), intent(in)      ::  icorr                            !<  correlation type
  character(len=1), intent(in)      ::  ispp                             !<  spin polarization  ' ', 's', 'r'

  integer, intent(in)               ::  nr                               !<  number of radial grid points
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  integer, intent(in)               ::  norb                             !<  number of orbitals
  integer, intent(in)               ::  ncore                            !<  number of orbitals treated as core
  integer, intent(in)               ::  lo(mxdorb)                       !<  angular quantum number l
  integer, intent(in)               ::  iso(mxdorb)                      !<  spin quantum number
  real(REAL64), intent(in)          ::  zo(mxdorb)                       !<  orbital occupation

  real(REAL64), intent(in)          ::  znuc                             !<  nuclear charge
  real(REAL64), intent(in)          ::  zel                              !<  electron charge

  real(REAL64), intent(in)          ::  cfac                             !<  core correction factor
  real(REAL64), intent(in)          ::  rcfac                            !<  core correction radius

  real(REAL64), intent(in)          ::  vscr(mxdnr,0:lc,-1:1)            !<  screened-pseudo-potential in Rydberg

! output

  real(REAL64), intent(out)         ::  vpsd(mxdnr,0:lc,-1:1)            !<  pseudo-potential in Rydberg

  real(REAL64), intent(out)         ::  zratio
  real(REAL64), intent(out)         ::  zion                             !<  pseudo-potential core charge

! input and output

  real(REAL64), intent(inout)       ::  cdpsd(mxdnr,-1:1)                !<  valence charge density (down or total)
  real(REAL64), intent(inout)       ::  cdc(nr)                          !<  core charge density


! local allocatable arrays

  real(REAL64), allocatable         ::  vhxc(:,:)                        !<  screening potential in Rydberg

  real(REAL64), allocatable         ::  rcut(:)


! local variables

  integer            ::  ncp

  real(REAL64)       ::  etot(10)                                        !  components of total energy

  real(REAL64)       ::  zval, zval2
  real(REAL64)       ::  ac, bc, cc
  integer            ::  icore
  character(len=1)   ::  blank

  real(REAL64)       ::  fcut

  real(REAL64)       ::  vp2z, cdfac

  integer            ::  jcut

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  ECUT = 0.001_REAL64

! counters

  integer     ::  i, j, l


! initializes vpsd

  do i =-1,1
  do l = 0,lc
  do j = 1,mxdnr
    vpsd(j,l,i) = ZERO
  enddo
  enddo
  enddo

! Reset the n quantum numbers to include all valence orbitals.
! Compute the ratio between the valence charge present and the
! valence charge of a neutral atom.
! Transfer pseudo valence charge to charge array

  allocate(rcut(norb))

  ncp = ncore + 1
  zval = ZERO
  zratio = ZERO
  do i = ncp,norb
    zval = zval + zo(i)
  enddo
  zion = zval + znuc - zel
  if (zval  /=  ZERO) zratio = zion/zval


!!!!!!!
!mmga

!  If a core correction is indicated construct pseudo core charge
!  cdc(r) = r^2*exp(ac+bc*r^2+cc*r^4) inside r(icore)
!  if cfac < 0 or the valence charge is zero the full core is used

   if (ifcore  /=  0) then
     ac = ZERO
     bc = ZERO
     cc = ZERO
     icore = 1
     if (cfac > ZERO .and. zratio /= ZERO) then

       if (rcfac <= ZERO) then
         do i = nr,2,-1
           icore = i
           if(ispp == ' ') then
             cdfac = cfac*zratio*cdpsd(i,0)
           else
             cdfac = cfac*zratio*(cdpsd(i, 1)+cdpsd(i,-1))
           endif
           if (cdc(i) > cdfac) exit
         enddo
       else
         do i = nr,2,-1
           icore = i
           if (r(i) <= rcfac ) exit
         enddo
       endif

       call atom_psd_pcc_exp(nr, icore, ac, bc, cc, r, cdc)

     endif

     write(iowrite,*)
     write(iowrite,*)
     write(iowrite,'(" core correction used")')
     write(iowrite,'(" pseudo core inside r =",f7.3)') r(icore)
     write(iowrite,'(" ac =",f7.3," bc =",f7.3," cc =",f7.3)') ac, bc, cc
     write(iowrite,*)
  endif

!mmga
!!!!!!

! End the pseudo core charge.
! Compute the potential due to pseudo valence charge.

!njtj  ***  NOTE  ***
! Spin-polarized potentails should be unscreend with
! spin-polarized valence charge.  This was not
! done in pseudo and pseudok in earlier versions
! of this program.
!njtj    ***  NOTE  ***

  if (ispp == 's') then
    blank='s'
  else
    blank=' '
  endif

  allocate(vhxc(mxdnr,-1:1))

  zval2 = zval
  call atom_atm_velect(0, 1, icorr, blank, ifcore,                       &
      nr, r, drdi, zval, cdpsd, cdc, vhxc, etot,                         &
      iowrite, mxdnr)

  if (ifcore == 2) zion = zion + zval - zval2

! Construct the ionic pseudopotential and find the cutoff,
! ecut should be adjusted to give a reasonable ionic cutoff
! radius, but should not alter the pseudopotential, ie.,
! the ionic cutoff radius should not be inside the pseudopotential
! cutoff radius

  do i = ncp,norb

    do j = 2,nr
      vpsd(j,lo(i),iso(i)) = (vscr(j,lo(i),iso(i)) - vhxc(j,iso(i)))*r(j)
    enddo

    do j = 2,nr
      vp2z = vpsd(j,lo(i),iso(i)) + 2*zion
      if (abs(vp2z) > ECUT) jcut = j
    enddo
    rcut(i-ncore) = r(jcut)
    do j = jcut,nr
      fcut = exp(-5*(r(j)-r(jcut)))
      vpsd(j,lo(i),iso(i)) = - 2*zion + fcut * (vpsd(j,lo(i),iso(i))+2*zion)
    enddo

  enddo

  deallocate(rcut)
  deallocate(vhxc)

  return
end subroutine atom_psd_unscreen
