!>  Implements the Perdew-Zunger parameterization of Ceperley-Alder exchange and
!>  correlation. Ref: Perdew & Zunger, Phys. Rev. B 23 5075 (1981).
!>
!>  \author       J.M.Soler
!>  \version      6.013
!>  \date         January 97, 22 June 2021
!>  \copyright    GNU Public License v2

subroutine atom_xc_pzc( nspin, dens, epsc, vc )

! *****************************************************************
!  Adapted by J.M.Soler from routine velect of Froyen's
!    pseudopotential generation program. Madrid, Jan'97. Version 0.5.
! **** Units *******************************************************
! Densities in electrons/Bohr**3
! Energies in Hartrees
! *****************************************************************

! converted to f90, April 2018
! removed exchange
! cleanup and new interface, July 2019. JLM
! jlm  version 6.00


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  nspin                            !<  nspin=1 => unpolarized; nspin=2 => polarized

  real(REAL64), intent(in)          ::  dens(nspin)                      !<  charge density

! output

  real(REAL64), intent(out)         ::  epsc                             !<  correlation energy densities
  real(REAL64), intent(out)         ::  vc(nspin)                        !<  correlation potential

! local

  real(REAL64)    ::  d                                                  !  total density
  real(REAL64)    ::  z, fz, fzp                                         !  polarization etc...
  real(REAL64)    ::  rs                                                 !  Wigner parameter
  logical         ::  lzero

  real(REAL64)    ::  sqrs, rslog
  real(REAL64)    ::  te, be

  real(REAL64)    ::  epscp, vcp
  real(REAL64)    ::  epscf, vcf

! constants

  real(REAL64), parameter   ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  real(REAL64), parameter   ::  PI   = 3.14159265358979312_REAL64
  real(REAL64), parameter   ::  TRD  = ONE/3,  FTRD = 4*TRD

  real(REAL64), parameter   ::  TFTM = (2*ONE)**FTRD - 2*ONE
  real(REAL64), parameter   ::  CRS  = (3 / (4*PI))**TRD

  real(REAL64), parameter   ::  CON2 = 0.008_REAL64 / 3
  real(REAL64), parameter   ::  CON3 = 0.3502_REAL64 / 3
  real(REAL64), parameter   ::  CON4 = 0.0504_REAL64 / 3
  real(REAL64), parameter   ::  CON5 = 0.0028_REAL64 / 3
  real(REAL64), parameter   ::  CON6 = 0.1925_REAL64 / 3
  real(REAL64), parameter   ::  CON7 = 0.0206_REAL64 / 3
  real(REAL64), parameter   ::  CON8 = 9.7867_REAL64 / 6
  real(REAL64), parameter   ::  CON9 = 1.0444_REAL64 / 3
  real(REAL64), parameter   ::  CON10 = 7.3703_REAL64 / 6
  real(REAL64), parameter   ::  CON11 = 1.3336_REAL64 / 3

  real(REAL64), parameter   ::  C1P053 = 1.0529_REAL64
  real(REAL64), parameter   ::  C3334 = 0.3334_REAL64
  real(REAL64), parameter   ::  C2846 = 0.2846_REAL64
  real(REAL64), parameter   ::  C1P398 = 1.3981_REAL64
  real(REAL64), parameter   ::  C2611 = 0.2611_REAL64
  real(REAL64), parameter   ::  C1686 = 0.1686_REAL64
  real(REAL64), parameter   ::  C0622 = 0.0622_REAL64
  real(REAL64), parameter   ::  C004 = 0.004_REAL64
  real(REAL64), parameter   ::  C096 = 0.096_REAL64
  real(REAL64), parameter   ::  C0232 = 0.0232_REAL64
  real(REAL64), parameter   ::  C0311 = 0.0311_REAL64
  real(REAL64), parameter   ::  C0014 = 0.0014_REAL64
  real(REAL64), parameter   ::  C0538 = 0.0538_REAL64
  real(REAL64), parameter   ::  C0096 = 0.0096_REAL64

! counters

  integer         ::  isp


  lzero = .FALSE.

! Find density and polarization

  if (nspin == 2) then
    d = dens(1) + dens(2)

    if (d <= ZERO) then
      epsc = ZERO
      vc(1) = ZERO
      vc(2) = ZERO
      lzero = .TRUE.
    else
      z = (dens(1) - dens(2)) / d
      fz = ((1+z)**FTRD + (1-z)**FTRD-2) / TFTM
      fzp = FTRD*((1+z)**TRD - (1-z)**TRD) / TFTM
    endif
  else
    d = dens(1)
    if (d <= ZERO) then
      epsc = ZERO
      vc(1) = ZERO
      lzero = .TRUE.
    else
      z = ZERO
      fz = ZERO
      fzp = ZERO
    endif
  endif

  if(.not. lzero) then

    rs = CRS / d**TRD

!   Perdew-Zunger correlation

    if (rs > ONE) then

      sqrs = sqrt(rs)
      te = ONE + CON10*sqrs + CON11*rs
      be = ONE + C1P053*sqrs + C3334*rs
      epscp = -C2846 / be
      vcp = epscp*te / be
      te = ONE + CON8*sqrs + CON9*rs
      be = ONE + C1P398*sqrs + C2611*rs
      epscf = -C1686 / be
      vcf = epscf*te / be

    else

      rslog = log(rs)
      epscp= (C0622+C004*rs)*rslog - C096 - C0232*rs
      vcp = (C0622+CON2*rs)*rslog - CON3 - CON4*rs
      epscf = (C0311+C0014*rs)*rslog - C0538 - C0096*rs
      vcf = (C0311+CON5*rs)*rslog - CON6 - CON7*rs

    endif

!   Find up and down potentials

    if (nspin == 2) then
      epsc = epscp + fz*(epscf-epscp)
      vc(1) = vcp + fz*(vcf-vcp) + (1-z)*fzp*(epscf-epscp)
      vc(2) = vcp + fz*(vcf-vcp) - (1+z)*fzp*(epscf-epscp)
    else
      epsc = epscp
      vc(1) = vcp
    endif

!   Change from Rydbergs to Hartrees

    epsc = epsc / 2
    do isp = 1,nspin
      vc(isp) = vc(isp) / 2
    enddo

  endif

  return

end subroutine atom_xc_pzc
