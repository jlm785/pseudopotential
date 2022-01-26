!>  implements Perdew-Burke-Ernzerhof Generalized-Gradient-Approximation.
!>  Ref: J.P.Perdew, K.Burke & M.Ernzerhof, PRL 77, 3865 (1996)
!>
!>  \author       L.C.Balbas and J.M.Soler
!>  \version      6.013
!>  \date         December 96, 22 June 2021
!>  \copyright    GNU Public License v2

  subroutine atom_xc_pbe( irel, nspin, dens, gdens,                      &
                   ex, ec, dexdd, decdd, dexdgd, decdgd )

! *********************************************************************
! Written by L.C.Balbas and J.M.Soler. December 1996. Version 0.5.
! ******** INPUT ******************************************************
! INTEGER IREL           : Relativistic-exchange switch (0=No, 1=Yes)
! INTEGER NSPIN          : Number of spin polarizations (1 or 2)
! REAL*8  DENS(NSPIN)    : Total electron density (if NSPIN=1) or
!                           spin electron density (if NSPIN=2)
! REAL*8  GDENS(3,NSPIN) : Total or spin density gradient
! ******** OUTPUT *****************************************************
! REAL*8  EX             : Exchange energy density
! REAL*8  EC             : Correlation energy density
! REAL*8  DEXDD(NSPIN)   : Partial derivative
!                           d(DensTot*Ex)/dDens(ispin),
!                           where DensTot = Sum_ispin( DENS(ispin) )
!                          For a constant density, this is the
!                          exchange potential
! REAL*8  DECDD(NSPIN)   : Partial derivative
!                           d(DensTot*Ec)/dDens(ispin),
!                           where DensTot = Sum_ispin( DENS(ispin) )
!                          For a constant density, this is the
!                          correlation potential
! REAL*8  DEXDGD(3,NSPIN): Partial derivative
!                           d(DensTot*Ex)/d(GradDens(i,ispin))
! REAL*8  DECDGD(3,NSPIN): Partial derivative
!                           d(DensTot*Ec)/d(GradDens(i,ispin))
! ********* UNITS ****************************************************
! Lengths in Bohr
! Densities in electrons per Bohr**3
! Energies in Hartrees
! Gradient vectors in cartesian coordinates
! ********* ROUTINES CALLED ******************************************
! EXCHNG, PW92C
! ********************************************************************


! converted to f90, July 2018
! cleanup and new interface, July 2019. JLM
! jlm  version 6.00

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter          ::  iowrite = 6                       !  default tape for writing

! input

  integer, intent(in)               ::  irel                        !<  relativistic exchange? (0=>no, 1=>yes)

  integer, intent(in)               ::  nspin                       !<  nspin=1 => unpolarized; nspin=2 => polarized

  real(REAL64), intent(in)          ::  dens(nspin)                 !<  charge density
  real(REAL64), intent(in)          ::  gdens(3,nspin)              !<  charge density

! output

  real(REAL64), intent(out)         ::  ex                          !<  exchange energy density
  real(REAL64), intent(out)         ::  ec                          !<  correlation energy density
  real(REAL64), intent(out)         ::  dexdd(nspin)                !<  d rho*ex / d dens(spin)
  real(REAL64), intent(out)         ::  decdd(nspin)                !<  d rho*ec / d dens(spin)
  real(REAL64), intent(out)         ::  dexdgd(3,nspin)             !<  d rho*ex / d dgens(i,spin)
  real(REAL64), intent(out)         ::  decdgd(3,nspin)             !<  d rho*ec / d dgens(i,spin)

! local

!  real(REAL64)    ::  rho(nspin)

! constants

  real(REAL64), parameter   ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  real(REAL64), parameter   ::  PI   = 3.14159265358979312_REAL64
  real(REAL64), parameter   ::  THD  = ONE/3,  TWOTHD  = 2*THD

  real(REAL64), parameter   ::  DENMIN = 0.1_REAL64**12             !      lower bounds of density to avoid divisions by zero
  real(REAL64), parameter   ::  GDMIN = 0.1_REAL64**12              !      lower bounds of density gradient to avoid divisions by zero

  real(REAL64), parameter   ::  BETA = 0.066725_REAL64
  real(REAL64), parameter   ::  XGAMMA = (ONE - log(2*ONE)) / (PI*PI)
  real(REAL64), parameter   ::  XMU = BETA * PI*PI / 3
  real(REAL64), parameter   ::  XKAPPA = 0.804_REAL64

! counters

  integer         ::  is, ix

! internal variables

  real(REAL64)  ::  a,  d(2), dadd, decudd
  real(REAL64)  ::  df1dd, df2dd, df3dd, df4dd, df1dgd, df3dgd, df4dgd
  real(REAL64)  ::  dfcdd(2), dfcdgd(3,2), dfdd, dfdgd, dfxdd(2), dfxdgd(3,2)
  real(REAL64)  ::  dhdd, dhdgd, dkfdd, dksdd, dpdd, dpdz, drsdd
  real(REAL64)  ::  ds, dsdd, dsdgd, dt, dtdd, dtdgd, dzdd(2)
  real(REAL64)  ::  ecunif, exunif
  real(REAL64)  ::  f, f1, f2, f3, f4, fc, fx
  real(REAL64)  ::  gd(3,2), gdm(2), gdms, gdmt, gds, gdt(3)
  real(REAL64)  ::  h,  kf, kfs, ks,  phi,  rs, s
  real(REAL64)  ::  t,  vcunif(2), vxunif, zeta

! jlm    variables to avoid ifort interface complaints

    real(REAL64)  ::  dsjlm(2),vxjlm(2)

!  jlm

! translate density and its gradient to new variables

  if (nspin == 1) then
    d(1) = dens(1) / 2
    d(2) = d(1)
    dt = max( DENMIN, dens(1) )
    do ix = 1,3
      gd(ix,1) = gdens(ix,1) / 2
      gd(ix,2) = gd(ix,1)
      gdt(ix) = gdens(ix,1)
    enddo
  else
    d(1) = dens(1)
    d(2) = dens(2)
    dt = max( DENMIN, dens(1)+dens(2) )
    do ix = 1,3
      gd(ix,1) = gdens(ix,1)
      gd(ix,2) = gdens(ix,2)
      gdt(ix) = gdens(ix,1) + gdens(ix,2)
    enddo
  endif
  gdm(1) = gd(1,1)*gd(1,1) + gd(2,1)*gd(2,1) + gd(3,1)*gd(3,1)
  gdm(1) = sqrt(gdm(1))
  gdm(2) = gd(1,2)*gd(1,2) + gd(2,2)*gd(2,2) + gd(3,2)*gd(3,2)
  gdm(2) = sqrt(gdm(2))
  gdmt   = gdt(1)*gdt(1)  + gdt(2)*gdt(2)  + gdt(3)*gdt(3)
  gdmt   = sqrt(gdmt)
  gdmt = max( GDMIN, gdmt )

! find local correlation energy and potential

  call atom_xc_pw92c( 2, d, ecunif, vcunif )
!  call pw92c( 2, d, ecunif, vcunif )

! find total correlation energy

  rs = ( 3 / (4*PI*dt) )**THD
  kf = (3 * PI*PI * dt)**THD
  ks = sqrt( 4 * kf / PI )
  zeta = ( d(1) - d(2) ) / dt
  zeta = max( -ONE+DENMIN, zeta )
  zeta = min(  ONE-DENMIN, zeta )
  phi = ( (1+zeta)**TWOTHD + (1-zeta)**TWOTHD ) / 2
  t = gdmt / (2 * phi * ks * dt)
  f1 = ecunif / XGAMMA / phi**3
  f2 = exp(-f1)
  a = BETA / XGAMMA / (f2-1)
  f3 = t*t + a * t*t*t*t
  f4 = BETA/XGAMMA * f3 / (1 + a*f3)
  h = XGAMMA * phi*phi*phi * log( 1 + f4 )
  fc = ecunif + h

! find correlation energy derivatives

  drsdd = - THD * rs / dt
  dkfdd =   THD * kf / dt
  dksdd = ks * dkfdd / (2 * kf)
  dzdd(1) =   1 / dt - zeta / dt
  dzdd(2) = - 1 / dt - zeta / dt
  dpdz = THD * ( 1/(1+zeta)**THD - 1/(1-zeta)**THD )

  do is = 1,2
    decudd = ( vcunif(is) - ecunif ) / dt
    dpdd = dpdz * dzdd(is)
    dtdd = - t * ( dpdd/phi + dksdd/ks + 1/dt )
    df1dd = f1 * ( decudd/ecunif - 3*dpdd/phi )
    df2dd = - f2 * df1dd
    dadd = - a * df2dd / (f2-1)
    df3dd = (2*t + 4*a*t*t*t) * dtdd + dadd * t*t*t*t
    df4dd = f4 * ( df3dd/f3 - (dadd*f3+a*df3dd)/(1+a*f3) )
    dhdd = 3 * h * dpdd / phi
    dhdd = dhdd + XGAMMA * phi*phi*phi * df4dd / (1+f4)
    dfcdd(is) = vcunif(is) + h + dt * dhdd

    do ix = 1,3
      dtdgd = (t / gdmt) * gdt(ix) / gdmt
      df3dgd = dtdgd * ( 2 * t + 4 * a * t**3 )
      df4dgd = f4 * df3dgd * ( 1/f3 - a/(1+a*f3) )
      dhdgd = XGAMMA * phi*phi*phi * df4dgd / (1+f4)
      dfcdgd(ix,is) = dt * dhdgd
    enddo
  enddo

! find exchange energy and potential

  fx = 0
  do is = 1,2
    ds   = max( DENMIN, 2 * d(is) )
    gdms = max( GDMIN, 2 * gdm(is) )
    kfs = (3 * PI*PI * ds)**THD
    s = gdms / (2 * kfs * ds)
    f1 = 1 + XMU * s*s / XKAPPA
    f = 1 + XKAPPA - XKAPPA / f1
!   jlm
    dsjlm(1) = ds
    call atom_xc_x( irel, 1, dsjlm, exunif, vxjlm )
!    call xc_x( irel, 1, ds, exunif, vxunif )
!    call exchng( irel, 1, dsjlm, exunif, vxjlm )
!    call exchng( irel, 1, ds, exunif, vxunif )
    vxunif = vxjlm(1)
!   jlm
    fx = fx + ds * exunif * f

    dkfdd = THD * kfs / ds
    dsdd = s * ( -dkfdd/kfs - 1/ds )
    df1dd = 2 * (f1-1) * dsdd / s
    dfdd = XKAPPA * df1dd / (f1*f1)
    dfxdd(is) = vxunif * f + ds * exunif * dfdd

    do ix = 1,3
      gds = 2 * gd(ix,is)
      dsdgd = (s / gdms) * gds / gdms
      df1dgd = 2 * XMU * s * dsdgd / XKAPPA
      dfdgd = XKAPPA * df1dgd / f1**2
      dfxdgd(ix,is) = ds * exunif * dfdgd
    enddo

  enddo
  fx = fx / (2 * dt)

! set output arguments

  ex = fx
  ec = fc
  do is = 1,nspin
    dexdd(is) = dfxdd(is)
    decdd(is) = dfcdd(is)
    do ix = 1,3
      dexdgd(ix,is) = dfxdgd(ix,is)
      decdgd(ix,is) = dfcdgd(ix,is)
    enddo
  enddo

  return

end subroutine atom_xc_pbe
