!>  Calculates the kinetic energy density for a given orbital
!>
!>  \author       Emmi Gareis
!>  \version      6.1.0
!>  \date         2025, 3 January 2026.
!>  \copyright    GNU Public License v2

subroutine atom_atm_tau(ispp, ar, br, nr, r,                             &
      lo, iso, ev, v, tau,                                               &
      mxdnr)

! Written in atom_atm_orban by Emmi Gareis in collaboration with Ingo Schelter and Stephan Kümmel.
! Separated in an autonomous subroutine, 3 January 2025. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points

  character(len=1), intent(in)      ::  ispp                             !<  spin polarization  ' ', 's', 'r'

  real(REAL64), intent(in)          ::  ar(mxdnr)                        !<  radial wave-function u = rR(r)  (integral ar**2 = 1) or major component radial wave-function u = rR(r)
  real(REAL64), intent(in)          ::  br(mxdnr)                        !<  d ar / d r    or  minor component (integral ar**2 + br**2 = 1)

  integer, intent(in)               ::  nr                               !<  dimension of the number of radial points

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points

  integer, intent(in)               ::  lo                               !<  angular quantum number l
  integer, intent(in)               ::  iso                              !<  2*spin or 2*(j-l)

  real(REAL64), intent(in)          ::  ev                               !<  orbital energy

  real(REAL64), intent(in)          ::  v(mxdnr)                         !<  KS-potential

! output

  real(REAL64), intent(out)         ::  tau(mxdnr)                       !<  orbital core kinetic energy density (gradient version)

! local alocatable arrays

  real(REAL64), allocatable         ::  arder(:)                         !  d ar / d r derivative of major spinor component
  real(REAL64), allocatable         ::  brder(:)                         !  d br / d r derivative of minor spinor component

! local

  real(REAL64)           ::  jo                                          !  spin-orbit-coupling related quantum number j
  integer                ::  lambda                                      !  eigenvalue of the spin-orbit-coupling related operator K acting on the spherical spinors
  real(REAL64)           ::  aror2                                       !  ar / r^2

  real(REAL64)           ::  tauSecondTermValue                          !  contribution to tau from l-quantum number from theta and phi gradient contribution

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  AI = 137.035999177_REAL64               ! speed of light in hartree atomic units (1/alpha)
  real(REAL64), parameter    ::  PI = 4*atan(ONE)

! counters

  integer          ::  i

  tauSecondTermValue = (((lo+1)*(lo+1) + lo+1)*(2*lo+3)) / (4*PI)

  allocate(arder(mxdnr))
  allocate(brder(mxdnr))

  if (ispp == 'r') then

    jo = abs(lo + (iso * ONE)/2)

    if (iso == 1) then              !relativistic case for spin +1/2 orbital
      lambda = -(lo + 1)
    else                            !relativistic case for spin -1/2 orbital
      lambda = lo
    endif
    do i = 2, nr

!     contribution of i-th orbital to tau in the relativistic case

      arder(i) = -(lambda/r(i))*ar(i) + (ev - v(i) + 4*AI*AI) * (br(i)/(2*AI))
      brder(i) = +(lambda/r(i))*br(i) - (ev - v(i)) * (ar(i)/(2*AI))
      tau(i) = (AI* (2*jo +1)/(4*PI)) *                                        &
                 (  ((ar(i)/(r(i)*r(i)))*(-brder(i) + lambda*br(i)/r(i)))      &
                  + ((br(i)/(r(i)*r(i)))*( arder(i) + lambda*ar(i)/r(i)))  )
    enddo

  else

    do i = 2,nr

!     contribution of i-th orbital to tau in the nonrelativistic case

      aror2 = ar(i)/(r(i)*r(i))
      tau(i) = ((br(i)/r(i) - aror2)*(br(i)/r(i) - aror2))*(2*lo +1)/(4*PI) +  &
              aror2*aror2*tauSecondTermValue
    enddo

  endif

  tau(1) = ZERO

  deallocate(arder)
  deallocate(brder)

  return

end subroutine atom_atm_tau
