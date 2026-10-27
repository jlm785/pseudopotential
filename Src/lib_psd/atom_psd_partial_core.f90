!>  Applies the partial core correction.
!>  radius is when rho_core = cfac*rho_pseudo
!>
!>  \author       Norm Troullier, Jose Luis Martins, Manuel Maria Alemany
!>  \version      6.1.0
!>  \date         24 March 2026, 27 October 2026.
!>  \copyright    GNU Public License v2


subroutine atom_psd_partial_core(ifcore, ispp, cfac, rcfac, cfactau,     &
    nr, r, drdi,                                                         &
    norb, ncore, zo, znuc, zel, zratio, zion, zval,                      &
    cdpsd, cdc_total, cdc, tauc_total, tauc,                             &
    iowrite, mxdnr, mxdorb)

! extracted from atom_psd_unscreen, 5 February 2026. JLM
! added flexibility to tau correction radius. 24 March 2026. JLM
! full core kinetic energy density if cfactau < 0 (as cfac). 27 October 2026. claude

!mmga  modifications from early Sverre code by Manuel Maria Gonzalez Alemany

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points
  integer, intent(in)               ::  mxdorb                           !<  dimension of the number of orbitals

  integer, intent(in)               ::  iowrite                          !<  default tape for writing


  integer, intent(in)               ::  ifcore                           !<  0: no partial core correction; 1 : xc partial core correction; 2 : not recommended.
  character(len=1), intent(in)      ::  ispp                             !<  spin polarization  ' ', 's', 'r'

  integer, intent(in)               ::  nr                               !<  number of radial grid points
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  real(REAL64), intent(in)          ::  cdpsd(mxdnr,-1:1)                !<  valence charge density (down or total)
  real(REAL64), intent(in)          ::  cdc_total(mxdnr)                 !<  total core charge density
  real(REAL64), intent(in)          ::  tauc_total(mxdnr)                !<  core kinetic energy density (total)

  integer, intent(in)               ::  norb                             !<  number of orbitals
  integer, intent(in)               ::  ncore                            !<  number of orbitals treated as core
  real(REAL64), intent(in)          ::  zo(mxdorb)                       !<  orbital occupation

  real(REAL64), intent(in)          ::  znuc                             !<  nuclear charge
  real(REAL64), intent(in)          ::  zel                              !<  electron charge

  real(REAL64), intent(in)          ::  cfac                             !<  core correction factor for charge density
  real(REAL64), intent(in)          ::  rcfac                            !<  default core correction radius
  real(REAL64), intent(in)          ::  cfactau                          !<  core correction factor for kinetic energy density

! output

  real(REAL64), intent(out)         ::  zratio
  real(REAL64), intent(out)         ::  zion                             !<  pseudo-potential core charge
  real(REAL64), intent(out)         ::  zval                             !<  valence charge

  real(REAL64), intent(out)         ::  tauc(mxdnr)                      !<  core kinetic energy density (partial or total)

  real(REAL64), intent(out)         ::  cdc(mxdnr)                       !<  4 pi r^2 core charge density (total or partial)

! local variables

  integer            ::  ncp

  real(REAL64)       ::  ac, bc, cc
  integer            ::  icore

  real(REAL64)       ::  cdfac
  real(REAL64)       ::  xsum

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  PI = 4*atan(ONE)
  real(REAL64), parameter    ::  EPS = 0.001_REAL64

! counters

  integer     ::  i

! Compute the ratio between the valence charge present and the
! valence charge of a neutral atom.

  ncp = ncore + 1
  zval = ZERO
  zratio = ZERO
  do i = ncp,norb
    zval = zval + zo(i)
  enddo
  zion = zval + znuc - zel
  if (zval  /=  ZERO) zratio = zion/zval

  do i = 1,nr
    cdc(i) = cdc_total(i)
  enddo

!!!!!!!
!mmga

!  If a core correction is indicated construct pseudo core charge
!  cdc_partial(r) = r^2*exp(ac+bc*r^2+cc*r^4) inside r(icore)
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
          if (cdc_total(i) > cdfac) exit
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

! core kinetic energy

  xsum = ZERO
  do i = 1,nr
    tauc(i) = 4*PI*tauc_total(i)*r(i)*r(i)
    xsum = xsum + tauc(i)*drdi(i)
  enddo

! if cfactau < 0 or the valence charge is zero the full core kinetic
! energy density is used (tests with very large cutoffs)

  if(xsum > EPS) then

    write(iowrite,'(" core kinetic energy =",f10.2)') xsum

    if(cfactau > ZERO .and. zratio /= ZERO) then

      do i = nr,2,-1
        icore = i
        if(ispp == ' ') then
          cdfac = cfactau*zratio*cdpsd(i,0)
        else
          cdfac = cfactau*zratio*(cdpsd(i, 1)+cdpsd(i,-1))
        endif
        if (cdc_total(i) > cdfac) exit
      enddo

      call atom_psd_pcc_exp(nr, icore, ac, bc, cc, r, tauc)

      write(iowrite,*)
      write(iowrite,'(" kinetic energy density truncated")')
      write(iowrite,'(" inside r =",f7.3)') r(icore)
      write(iowrite,'(" ac =",f7.3," bc =",f7.3," cc =",f7.3)') ac, bc, cc
      write(iowrite,*)

      xsum = ZERO
      do i = 1,nr
        xsum = xsum + tauc(i)*drdi(i)
      enddo

      write(iowrite,'(" truncated core kinetic energy =",f10.2)') xsum
      write(iowrite,*)

    else

      write(iowrite,*)
      write(iowrite,'(" full core kinetic energy density used")')
      write(iowrite,*)

    endif

  endif

  return

end subroutine atom_psd_partial_core
