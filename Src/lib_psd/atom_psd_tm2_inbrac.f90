!>  Looks for a root of v''(0) inside the initial interval (x1,x2)
!>  for gamma when v''(0) has the same sign at both ends,
!>  which happens when the interval contains two roots.
!>  On success (x1,x2) is narrowed to the sign change closest to x1,
!>  the initial guess of gamma.
!>
!>  \author       claude+J.L.Martins
!>  \version      6.1.0
!>  \date         26 October 2026.
!>  \copyright    GNU Public License v2

subroutine atom_psd_tm2_inbrac(x1, x2, lfound,                           &
        jrc, r, drdi, bjin, rcpn, cdrc, lo, rpsi, delta, bj,             &
        mxdnr)

! written for atom_psd_tm2_eqsolve, the outward expansion of
! atom_psd_tm2_zrbac cannot find roots inside the interval. 26 October 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension radial grid

  integer, intent(in)               ::  jrc                              !<  r(jrc) = r_cut

  integer, intent(in)               ::  lo                               !<  angular momentum l
  real(REAL64), intent(in)          ::  cdrc                             !<  int_0^r_cut rho(r) dr

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  real(REAL64), intent(in)          ::  bjin(5)                          !<  right member of the equation
  real(REAL64), intent(in)          ::  rcpn(12)                         !<  r_cut^n

! output

  logical, intent(out)              ::  lfound                           !<  a sign change was found inside (x1,x2)
  real(REAL64), intent(out)         ::  delta                            !<  value of delta
  real(REAL64), intent(out)         ::  bj(5)                            !<  polyr coefficients

! input and output

  real(REAL64), intent(inout)       ::  x1, x2                           !<  values of gamma, x1 is the initial guess

  real(REAL64), intent(inout)       ::  rpsi(mxdnr)                      !<  r*psi

! local variables

  real(REAL64)      ::  f1, f2                                           !  v''(0) at the ends
  real(REAL64)      ::  xa, fa, xb, fb                                   !  consecutive sample points
  real(REAL64)      ::  xguess, dist, distmin
  real(REAL64)      ::  xlo, xhi
  logical           ::  la                                               !  gamfn succeeded at xa

  real(REAL64)      ::  x1_in, x2_in                                     !  original interval

  integer    ::  ifl

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  integer, parameter         ::  NSAMP = 20                              !  number of subintervals

! counters

  integer  ::  j


  lfound = .FALSE.

  x1_in = x1
  x2_in = x2

  if(x1 == x2) return

  call atom_psd_tm2_gamfn(x1, f1, ifl,                                   &
      jrc, r, drdi, bjin, rcpn, cdrc, lo, rpsi, delta, bj,               &
      mxdnr)
  if(ifl /= 0) return

  call atom_psd_tm2_gamfn(x2, f2, ifl,                                   &
      jrc, r, drdi, bjin, rcpn, cdrc, lo, rpsi, delta, bj,               &
      mxdnr)
  if(ifl /= 0) return

! the interval already brackets a root, nothing to do

  if(f1*f2 < ZERO) return

! samples the interval, keeps the sign change closest to the guess

  xguess = x1
  distmin = huge(distmin)

  xa = x1
  fa = f1
  la = .TRUE.
  do j = 1,NSAMP
    if(j == NSAMP) then
      xb = x2
      fb = f2
      ifl = 0
    else
      xb = x1 + j*(x2-x1)/NSAMP
      call atom_psd_tm2_gamfn(xb, fb, ifl,                               &
          jrc, r, drdi, bjin, rcpn, cdrc, lo, rpsi, delta, bj,           &
          mxdnr)
    endif

    if(ifl == 0) then
      if(la .and. fa*fb < ZERO) then
        dist = abs((xa+xb)/2 - xguess)
        if(dist < distmin) then
          distmin = dist
          xlo = xa
          xhi = xb
          lfound = .TRUE.
        endif
      endif
      xa = xb
      fa = fb
      la = .TRUE.
    else
      la = .FALSE.
    endif
  enddo

  if(lfound) then
    x1 = xlo
    x2 = xhi
  write(6,*)
  write(6,'("   The interval ",2g12.3," contains an even number of roots")') x1_in, x2_in
  write(6,'("   Check the wave-functions for l = ",i5," and check if core radius should be larger")') lo
  write(6,*)
  endif

  return

end subroutine atom_psd_tm2_inbrac
