!>  Implements the Perdew-Wang'92 local correlation (beyond rpa).
!>  ref: J.P.Perdew & Y.Wang, PRB, 45, 13244 (1992)
!>
!>  \author       L.C.Balbas and J.M.Soler
!>  \version      6.013
!>  \date         December 96, 22 June 2021
!>  \copyright    GNU Public License v2

subroutine atom_xc_pw92c( nspin, dens, ec, vc )

! ********************************************************************
! Written by L.C.Balbas and J.M.Soler. dec'96.  version 0.5.
! ********* Units ****************************************************
! densities in electrons per bohr**3
! energies in hartrees
! ********************************************************************

! converted to f90, June 2018
! cleanup and new interface, July 2019. JLM
! jlm  version 6.00


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  nspin                            !<  nspin=1 => unpolarized; nspin=2 => polarized

  real(REAL64), intent(in)          ::  dens(nspin)                      !<  charge density

! output

  real(REAL64), intent(out)         ::  ec                               !<  correlation energy densities
  real(REAL64), intent(out)         ::  vc(nspin)                        !<  correlation potential

! local

  real(REAL64)    ::  b, c
  real(REAL64)    ::  rs, zeta                                           !  Wigner-Seitz radius, polarization
  real(REAL64)    ::  f, g(0:2)
  real(REAL64)    ::  dtot

  real(REAL64)    ::  dbdrs, dcdrs
  real(REAL64)    ::  drsdd, dzdd(2)
  real(REAL64)    ::  decdrs, decdz, decdd(2)
  real(REAL64)    ::  dfdz, dgdrs(0:2)

! constants


  real(REAL64), parameter  ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  real(REAL64), parameter  ::  PI   = 3.14159265358979312_REAL64

  real(REAL64), parameter  ::  THD  = ONE/3,  FOUTHD = 4*THD
  real(REAL64), parameter  ::  CRS  = (3 / (4*PI))**THD
  real(REAL64), parameter  ::  CZ = ONE / ((2*ONE)**FOUTHD - 2*ONE)
  real(REAL64), parameter  ::  FPP0 = 8 * CZ / (9 * ONE)

! fix lower bound of density to avoid division by zero

  real(REAL64), parameter   ::  DENMIN = 0.1_REAL64 ** 12

! parameters from Table I of Perdew & Wang, PRB, 45, 13244 (92)

  real(REAL64), parameter   ::  P(0:2) = (/ ONE, ONE, ONE /)
  real(REAL64), parameter   ::  A(0:2) =                                 &
         (/ 0.031091_REAL64, 0.015545_REAL64, 0.016887_REAL64 /)
  real(REAL64), parameter   ::  ALPHA1(0:2) =                            &
         (/ 0.21370_REAL64, 0.20548_REAL64, 0.11125_REAL64  /)
  real(REAL64), parameter   ::  BETA1(0:2) =                             &
         (/ 7.5957_REAL64,  14.1189_REAL64,  10.357_REAL64 /)
  real(REAL64), parameter   ::  BETA2(0:2) =                             &
         (/ 3.5876_REAL64,   6.1977_REAL64,   3.6231_REAL64 /)
  real(REAL64), parameter   ::  BETA3(0:2) =                             &
         (/ 1.6382_REAL64,   3.3662_REAL64,   0.88026_REAL64 /)
  real(REAL64), parameter   ::  BETA4(0:2) =                             &
         (/ 0.49294_REAL64,  0.62517_REAL64,  0.49671_REAL64 /)

! counters

  integer         ::  ig


! find rs and zeta and derivatives drs/ddens and dzeta/ddens

  if (nspin  ==  1) then
    dtot = max( denmin, dens(1) )
    zeta = ZERO
    rs = CRS / dtot**THD
    drsdd = - THD * rs / dtot
    dzdd(1) = ZERO
  else
    dtot = max( denmin, dens(1)+dens(2) )
    zeta = ( dens(1) - dens(2) ) / dtot
    rs = CRS / dtot**THD
    drsdd = - THD * rs / dtot
    dzdd(1) =   ONE / dtot - zeta / dtot
    dzdd(2) = - ONE / dtot - zeta / dtot
  endif

! find eps_c(rs,0)=g(0), eps_c(rs,1)=g(1) and -alpha_c(rs)=g(2)
! using eq.(10) of cited reference (perdew & wang, prb, 45, 13244 (92))

  do ig = 0,2
    b = BETA1(ig) * sqrt(rs)    +                                        &
        BETA2(ig) * rs          +                                        &
        BETA3(ig) * rs*sqrt(rs) +                                        &
        BETA4(ig) * rs**(P(ig)+1)
    dbdrs = BETA1(ig) / (2*sqrt(rs))             +                       &
            BETA2(ig)                            +                       &
            BETA3(ig) * 3 * sqrt(rs) / (2 * ONE) +                       &
            BETA4(ig) * (P(ig)+1) * rs**P(ig)
    c = 1 + 1 / (2 * A(ig) * b)
    dcdrs = - (c-1) * dbdrs / b
    g(ig) = - 2 * A(ig) * ( 1 + ALPHA1(ig)*rs ) * log(c)
    dgdrs(ig) = - 2*A(ig) * ( ALPHA1(ig) * log(c) +                      &
                             (1+ALPHA1(ig)*rs) * dcdrs / c )
  enddo

! find f''(0) and f(zeta) from eq.(9)

  f = ( (1+zeta)**FOUTHD + (1-zeta)**FOUTHD - 2 ) * CZ
  dfdz = FOUTHD * ( (1+zeta)**THD - (1-zeta)**THD ) * CZ

! find eps_c(rs,zeta) from eq.(8)

  ec = g(0) - g(2) * f / FPP0 * (1-zeta**4) +                            &
     (g(1)-g(0)) * f * zeta**4
  decdrs = dgdrs(0) - dgdrs(2) * f / FPP0 * (1-zeta**4) +                &
         (dgdrs(1)-dgdrs(0)) * f * zeta**4
  decdz = - g(2) / FPP0 * ( dfdz*(1-zeta**4) - f*4*zeta**3 ) +           &
         (g(1)-g(0)) * ( dfdz*zeta**4 + f*4*zeta**3 )

! find correlation potential

  if (nspin  ==  1) then
    decdd(1) = decdrs * drsdd
    vc(1) = ec + dtot * decdd(1)
  else
    decdd(1) = decdrs * drsdd + decdz * dzdd(1)
    decdd(2) = decdrs * drsdd + decdz * dzdd(2)
    vc(1) = ec + dtot * decdd(1)
    vc(2) = ec + dtot * decdd(2)
  endif

end subroutine atom_xc_pw92c
