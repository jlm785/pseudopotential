!>  Estimates the error of the partial core kinetic energy density when its
!>  plane-wave expansion is truncated at |q| < k,
!>
!>  err(k) = sqrt( omega (1/(2 pi^2)) int_k^infinity q^2 |tau_c(q)|^2 dq )
!>
!>  By Parseval the integral is int |tau_c(r) - tau_c^(k)(r)|^2 d^3r, so err(k)
!>  is an upper bound (Cauchy-Schwarz) of int |tau_c(r) - tau_c^(k)(r)| d^3r
!>  in a volume omega (Ry).  It decreases monotonically with k.
!>  The volume should be comparable to the volume per atom in a solid.
!>
!>  \author       Jose Luis Martins, claude
!>  \version      6.2.0
!>  \date         16 November 2026.
!>  \copyright    GNU Public License v2

subroutine atom_kb_tauc_error(nql, delql, taucft, omega, err, nqmax)

! Written for the convergence plot of the partial core kinetic energy density. 16 November 2026. claude+JLM

! A sharp truncation in q leaves a ripple ~ cos(kr)/r^2 in tau_c(r), so the integral
! of |tau_c - tau_c^(k)| over all space diverges and a finite volume omega is needed.
! The tail beyond the last point of the grid (q_max) is not zero.  It is estimated
! with the asymptotic q^-5 decay of tau_c(q) (third derivative discontinuous at the
! radius of the partial core), |tau_c|^2 = A / q^10 with A averaged over the last
! QFIT of the grid.  If tau_c(q) decays faster the estimate is too large (conservative).
! With the r**6 term of the partial core (lr6 in atom_psd_partial_core) the decay is q^-6:
! |tau_c|^2 = A / q^12, q**12 in afit and afit / (9*qmax**9) for the tail.

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  nqmax                            !<  dimension of fourier grid points

  integer, intent(in)               ::  nql                              !<  number of points for the Fourier grid
  real(REAL64), intent(in)          ::  delql                            !<  spacing of points in Fourier grid
  real(REAL64), intent(in)          ::  taucft(0:nqmax)                  !<  Fourier transform of 4*pi*r**2 kinetic energy density of core (partial)
  real(REAL64), intent(in)          ::  omega                            !<  reference volume (bohr^3)

! output

  real(REAL64), intent(out)         ::  err(0:nqmax)                     !<  error of tau_c truncated at q(j) = j*delql (Ry)

! local variables

  real(REAL64)       ::  q, qmax
  real(REAL64)       ::  afit                                            !  amplitude of the q^-10 tail of |tau_c|^2
  real(REAL64)       ::  xsum
  integer            ::  jfit, nfit

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  real(REAL64), parameter    ::  PI = 4*atan(ONE)
  real(REAL64), parameter    ::  QFIT = 4.0_REAL64                       !  range of q used for the tail

! counters

  integer     ::  j


  qmax = nql*delql

! amplitude of the asymptotic tail

  jfit = max(1, nql - nint(QFIT/delql))
  nfit = 0
  afit = ZERO
  do j = jfit,nql
    q = j*delql
    afit = afit + (q**10)*taucft(j)*taucft(j)
    nfit = nfit + 1
  enddo
  afit = afit / nfit

! cumulative integral from the end of the grid (trapezoidal rule) plus the tail

  xsum = afit / (7*qmax**7)
  err(nql) = xsum
  do j = nql-1,0,-1
    q = j*delql
    xsum = xsum + delql*( q*q*taucft(j)*taucft(j)                        &
                       + (q+delql)*(q+delql)*taucft(j+1)*taucft(j+1) ) / 2
    err(j) = xsum
  enddo

  do j = 0,nql
    err(j) = sqrt(omega*err(j) / (2*PI*PI))
  enddo

  return

end subroutine atom_kb_tauc_error
