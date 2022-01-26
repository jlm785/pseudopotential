!>  Finds the eigenvalue ev, the wavefunction rpsi (r u(r))
!>  and the derivative drpsidr = d(rpsi)/dr.
!>  It uses an input guess for the eigenvalue (from atm_dsolv1 for example).
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.2
!>  \date         9 Septenber 2021
!>  \copyright    GNU Public License v2

subroutine atom_atm_difnrl(nr, r, drdi, d2rodr, v, rpsi, drpsidr,        &
      n, l, ev, iflag, tol,                                              &
      iowrite, mxdnr)

! interface to two subroutines depending if the requested state is
! bound.


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points

  integer, intent(in)               ::  nr                               !<  number of radial points

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i
  real(REAL64), intent(in)          ::  d2rodr(mxdnr)                    !<  (d^2 r(i) / d i^2) / (d r(i) / d i)

  real(REAL64), intent(in)          ::  v(mxdnr)                         !<  effective potential (RYDBERG)

  integer, intent(in)               ::  n                                !<  principal quantum number n
  integer, intent(in)               ::  l                                !<  angular quantum number l

  real(REAL64), intent(in)          ::  tol                              !<  precision of the solution tol ~1.0e-10

! output

  real(REAL64), intent(out)         ::  rpsi(mxdnr)                      !<  radial wave function
  real(REAL64), intent(out)         ::  drpsidr(mxdnr)                   !<  d rpsi / d i

  integer, intent(out)              ::  iflag                            !<  iflag = 0:success; iflag = 1: failed to converge; iflag > 3: major error

! input and output

  real(REAL64), intent(inout)       ::  ev                               !<  orbital energy, guess on input, accurate on output

! local variables

  real(REAL64)              ::  einf
  real(REAL64)              ::  revi
  integer                   ::  nrevi

  integer                   ::  nodes

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64

! counters

  integer  ::  j

! checks if the orbital is bound.


  einf = v(nr)
  revi = r(nr)

  call atom_atm_difnrl_one(nr, r, drdi, d2rodr, v, rpsi, drpsidr,        &
      l, einf, revi, nrevi, iflag,                                       &
      iowrite, mxdnr)

  if(nrevi /= nr .or. iflag /= 0) then
    write(6,*)
    write(6,*) '  atom_atm_difnrl: could not find if orbital is bound'
    write(6,*)

    return

  endif

  nodes = 0
  do j = 3,nr
    if(rpsi(j-1)*rpsi(j) < ZERO) nodes = nodes + 1
  enddo

  if(nodes > (n - l - 1)) then

    call atom_atm_difnrl_bound(nr, r, drdi, d2rodr, v, rpsi, drpsidr,    &
        n, l, ev, iflag, tol,                                            &
        iowrite, mxdnr)

  else

    ev = einf
    call atom_atm_difnrl_zr(nr, r, drdi, d2rodr, v, rpsi, drpsidr,       &
        n, l,  ev, iflag, revi, nrevi, tol,                              &
        iowrite, mxdnr)

  endif

  return

  end subroutine atom_atm_difnrl
