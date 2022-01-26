!>  Finds local exchange energy density and potential
!>
!>  \author       J.M.Soler
!>  \version      6.013
!>  \date         January 97, 22 June 2021
!>  \copyright    GNU Public License v2

subroutine atom_xc_x( irel, nspin, dens, epsx, vx )

! *****************************************************************
!  Adapted by J.M.Soler from routine velect of Froyen's
!    pseudopotential generation program. Madrid, Jan'97. Version 0.5.
! **** Input ******************************************************
! INTEGER IREL    : relativistic-exchange switch (0=no, 1=yes)
! INTEGER nspin     : spin-polarizations (1=>unpolarized, 2=>polarized)
! REAL*8  dens(nspin) : total (nsp=1) or spin (nsp=2) electron density
! **** Output *****************************************************
! REAL*8  epsx      : exchange energy density
! REAL*8  vx(nspin) : (spin-dependent) exchange potential
! **** Units ******************************************************
! Densities in electrons/Bohr**3
! Energies in Hartrees
! *****************************************************************

! Adapted for to f90, April 2018
! cleanup, July 2019. JLM
! jlm  version 6.00


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  irel                             !<  relativistic exchange? (0=>no, 1=>yes)
  integer, intent(in)               ::  nspin                            !<  nspin=1 => unpolarized; nspin=2 => polarized

  real(REAL64), intent(in)          ::  dens(nspin)                      !<  charge density

! output

  real(REAL64), intent(out)         ::  epsx                             !<  exchange energy densities
  real(REAL64), intent(out)         ::  vx(nspin)                        !<  exchange potential

! local

  real(REAL64)    ::  d                                                  !  total density
  real(REAL64)    ::  z, fz, fzp                                         !  polarization etc...
  real(REAL64)    ::  rs                                                 !  Wigner parameter
  logical         ::  lzero

  real(REAL64)    ::  alp                                                !  X-alpha parameter                                          !  total density

  real(REAL64)    ::  epsxp, vxp
  real(REAL64)    ::  epsxf, vxf

  real(REAL64)    ::  beta, sb, alb, cst

! constants

  real(REAL64), parameter   ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  real(REAL64), parameter   ::  PI   = 3.14159265358979312_REAL64
  real(REAL64), parameter   ::  TRD  = ONE/3,  FTRD = 4*TRD

  real(REAL64), parameter   ::  TFTM = (2*ONE)**FTRD - 2*ONE
  real(REAL64), parameter   ::  CRS  = (3 / (4*PI))**TRD
  real(REAL64), parameter   ::  A0 = (4/(9*PI))**TRD

  real(REAL64), parameter   ::  C014 = 0.014_REAL64


! X-alpha parameter:

  alp = 2 * TRD

  lzero = .FALSE.

! Find density and polarization

  if (nspin == 2) then
    d = dens(1) + dens(2)

    if (d <= ZERO) then
      epsx = ZERO
      vx(1) = ZERO
      vx(2) = ZERO
      lzero = .TRUE.
    else
      z = (dens(1) - dens(2)) / d
      fz = ((1+z)**FTRD + (1-z)**FTRD-2) / TFTM
      fzp = FTRD*((1+z)**TRD - (1-z)**TRD) / TFTM
    endif
  else
    d = dens(1)
    if (d <= ZERO) then
      epsx = ZERO
      vx(1) = ZERO
      lzero = .TRUE.
    else
      z = ZERO
      fz = ZERO
      fzp = ZERO
    endif
  endif

  if(.not. lzero) then

    rs = CRS / d**TRD

    vxp = -3*alp / (2*PI*A0*rs)
    epsxp = 3*vxp / (4*ONE)

    if (irel == 1) then
      beta = C014 / rs
      sb = sqrt(ONE+beta*beta)
      alb = log(beta+sb)
      vxp = vxp * (-ONE/2 + 3*alb / (2*beta*sb))
      cst = (beta*sb-alb) / (beta*beta)
      epsxp = epsxp * (ONE - 3*cst*cst / 2)
    endif

    vxf = (2*ONE)**TRD * vxp
    epsxf = (2*ONE)**TRD * epsxp
    if (nspin == 2) then
      vx(1) = vxp + fz*(vxf-vxp) + (1-z)*fzp*(epsxf-epsxp)
      vx(2) = vxp + fz*(vxf-vxp) - (1+z)*fzp*(epsxf-epsxp)
      epsx = epsxp + fz*(epsxf-epsxp)
    else
      vx(1) = vxp
      epsx = epsxp
    endif

  endif

  return

end subroutine atom_xc_x
