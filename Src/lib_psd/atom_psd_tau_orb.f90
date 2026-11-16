!>  Calculates 4 pi r^2 times the (positive definite) kinetic energy
!>  density of an occupied orbital with radial function u = rR(r),
!>  non-relativistic form, the derivative of u is calculated numerically
!>
!>  4 pi r^2 tau = zo [ (du/dr - u/r)^2 + l(l+1) (u/r)^2 ]
!>
!>  Its integral is zo times the kinetic energy (Ry) of the orbital,
!>  the same normalization as the core kinetic energy density of atom_atm_tau.
!>
!>  \author       Jose Luis Martins, claude
!>  \version      6.2.0
!>  \date         16 November 2026.
!>  \copyright    GNU Public License v2

subroutine atom_psd_tau_orb(nr, r, drdi, u, lo, zo, tau,                 &
      mxdnr)

! Written for the plots of the valence kinetic energy densities. 16 November 2026. claude+JLM

! The derivative du/dr = (du/di) / (dr/di) uses central differences in
! the grid index i (7 points, fewer near the ends of the grid).

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points

  integer, intent(in)               ::  nr                               !<  number of radial points
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  real(REAL64), intent(in)          ::  u(mxdnr)                         !<  radial wave-function u = rR(r) (integral u**2 = 1)
  integer, intent(in)               ::  lo                               !<  angular quantum number l
  real(REAL64), intent(in)          ::  zo                               !<  orbital occupation

! output

  real(REAL64), intent(out)         ::  tau(mxdnr)                       !<  4 pi r^2 kinetic energy density of the orbital

! local variables

  real(REAL64)           ::  dudi                                        !  d u / d i
  real(REAL64)           ::  dg                                          !  r d(u/r)/dr = du/dr - u/r
  real(REAL64)           ::  uor                                         !  u / r

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64

! counters

  integer          ::  i


  tau(1) = ZERO

  do i = 2,nr

    if(i > 3 .and. i < nr-2) then
      dudi = (      (u(i+3) - u(i-3)) - 9*(u(i+2) - u(i-2))              &
               + 45*(u(i+1) - u(i-1)) ) / 60
    elseif(i > 2 .and. i < nr-1) then
      dudi = ( -(u(i+2) - u(i-2)) + 8*(u(i+1) - u(i-1)) ) / 12
    elseif(i < nr) then
      dudi = (u(i+1) - u(i-1)) / 2
    else
      dudi = u(i) - u(i-1)
    endif

    uor = u(i) / r(i)
    dg = dudi / drdi(i) - uor

    tau(i) = zo * (dg*dg + (lo*(lo+1))*uor*uor)

  enddo

  return

end subroutine atom_psd_tau_orb
