!>  Integrates the relativistic 2 component "Dirac" equation.
!>  Finds the eigenvalue ev, the major and minor component
!>  of the wavefunction, ar and br. (r u(r))
!>  It uses an intial guess for the eigenvalues from dsolv1 or older value.
!>
!>  \author       Sverre Froyen, Norm Troullier, Jose Luis Martins
!>  \version      6.0.3
!>  \date         1980s, 12 September 2021.
!>  \copyright    GNU Public License v2

subroutine atom_atm_difrel(nr, r, drdi, v, ar, br,                       &
          n, l, iso, znuc, vzero, ev, iflag, tol,                        &
          iowrite, mxdnr)

! interface to two subroutines depending if the requested state is
! bound.
! so->iso, 12 September 2021. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points

  integer, intent(in)               ::  nr                               !<  number of radial points
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  real(REAL64), intent(in)          ::  v(mxdnr)                         !<  effective potential (RYDBERG)

  integer, intent(in)               ::  n                                !<  principal quantum number n
  integer, intent(in)               ::  l                                !<  angular quantum number l
  integer, intent(in)               ::  iso                              !<  2*spin or 2*(j-l)

  real(REAL64), intent(in)          ::  znuc                             !<  nuclear charge
  real(REAL64), intent(in)          ::  vzero                            !<  effective potential in Rydberg

  real(REAL64), intent(in)          ::  tol                              !<  precision of the solution tol ~1.0e-10

! output

  real(REAL64), intent(out)         ::  ar(mxdnr)                        !<  major component radial wave-function u = rR(r)
  real(REAL64), intent(out)         ::  br(mxdnr)                        !<  minor component (integral ar**2 + br**2 = 1)

  integer, intent(out)              ::  iflag                            !<  iflag = 0: success; iflag = 1: failed to converge; iflag > 3: major error

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

  call atom_atm_difrel_one(nr, r, drdi, v, ar, br,                       &
      l, iso, znuc, vzero, einf, revi, nrevi, iflag,                     &
      iowrite, mxdnr)

  if(nrevi /= nr .or. iflag /= 0) then
    write(6,*)
    write(6,*) '  atom_atm_difrel: could not find if orbital is bound'
    write(6,*)

    return

  endif

! nodes of major component

  nodes = 0
  do j = 3,nrevi
    if(ar(j-1)*ar(j) < ZERO) nodes = nodes + 1
  enddo

  if(nodes > (n - l - 1)) then

    call atom_atm_difrel_bound(nr, r, drdi, v, ar, br,                   &
          n, l, iso, znuc, vzero, ev, iflag, tol,                        &
          iowrite, mxdnr)

  else

    call atom_atm_difrel_zr(nr, r, drdi, v, ar, br,                      &
          n, l, iso, znuc, vzero, ev, iflag, revi, nrevi, tol,           &
          iowrite, mxdnr)

  endif

  return

end subroutine atom_atm_difrel
