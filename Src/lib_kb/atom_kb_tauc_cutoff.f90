!>  Fourier transform of the partial core kinetic energy density, the error
!>  of its plane-wave expansion truncated at |q| < k (atom_kb_tauc_error)
!>  in the volume of the sphere with 90% of the valence charge, and the
!>  smallest k for which that error is below a tolerance.
!>
!>  \author       Jose Luis Martins, claude
!>  \version      6.2.0
!>  \date         17 November 2026.
!>  \copyright    GNU Public License v2

subroutine atom_kb_tauc_cutoff(nr, r, drdi, cdv, tauc,                   &
         nql, delql, toltau, taucft, omega, err, qtau,                   &
         mxdnr, nqmax)

! Fourier transform moved from atom_kb_pot_four and volume from atom_kb_plot,
! cutoff of the partial core kinetic energy density. 17 November 2026. claude+JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points
  integer, intent(in)               ::  nqmax                            !<  dimension of fourier grid points

  integer, intent(in)               ::  nr                               !<  number of points in the radial grid
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points   r(i) = a*(exp(b*i)-1), i=1,...,nr
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i
  real(REAL64), intent(in)          ::  cdv(mxdnr)                       !<  4*pi*r**2 charge density of valence
  real(REAL64), intent(in)          ::  tauc(mxdnr)                      !<  4*pi*r**2 kinetic energy density of core (partial)

  integer, intent(in)               ::  nql                              !<  number of points for the Fourier grid
  real(REAL64), intent(in)          ::  delql                            !<  spacing of points in Fourier grid
  real(REAL64), intent(in)          ::  toltau                           !<  tolerance for the truncation error (Ry)

! output

  real(REAL64), intent(out)         ::  taucft(0:nqmax)                  !<  Fourier transform of 4*pi*r**2 kinetic energy density of core (partial)
  real(REAL64), intent(out)         ::  omega                            !<  volume of the sphere with a fraction FVAL of the valence charge (bohr^3)
  real(REAL64), intent(out)         ::  err(0:nqmax)                     !<  error of tau_c truncated at q(j) = j*delql (Ry)
  real(REAL64), intent(out)         ::  qtau                             !<  smallest q with err <= toltau (1/bohr)

! allocatable arrays

  real(REAL64), allocatable         ::  yp(:)

! local variables

  real(REAL64)       ::  rval, xsum, xtot                                !  radius with FVAL of the valence charge
  real(REAL64)       ::  qmax

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  real(REAL64), parameter    ::  PI = 4*atan(ONE)
  real(REAL64), parameter    ::  FVAL = 0.9_REAL64                       !  fraction of valence charge that defines the volume for the tau_c error

! counters

  integer     ::  j


! Fourier transform of the kinetic energy density

  allocate(yp(0:nqmax))

  do j = 0,nql
    yp(j) = j*delql
  enddo

  taucft = ZERO
!  call ft_radial(nr-1, nql, r, drdi, yp, tauc, taucft, 0)
  call ft_charge(nr-1, nql, r, drdi, yp, tauc, taucft)

! extrapolates to zero

  taucft(0) = (4*taucft(1) - taucft(2)) / (3*ONE)

  deallocate(yp)

! volume of the sphere with a fraction FVAL of the valence charge
! (r_WS ~ 0.84 r_val for elemental solids)

  xtot = ZERO
  do j = 1,nr
    xtot = xtot + cdv(j)*drdi(j)
  enddo
  xsum = ZERO
  rval = r(nr)
  do j = 1,nr
    xsum = xsum + cdv(j)*drdi(j)
    if(xsum > FVAL*xtot) then
      rval = r(j)
      exit
    endif
  enddo
  omega = 4*PI*rval**3 / 3

! error of tau_c(r) if its expansion is truncated at q

  call atom_kb_tauc_error(nql, delql, taucft, omega, err, nqmax)

! smallest q of the grid with err <= toltau (err decreases with q, zero without core).
! Beyond the grid the q^-5 tail of atom_kb_tauc_error gives err(q) = err(qmax) (qmax/q)**(7/2)

  qmax = nql*delql

  if(err(nql) > toltau) then
    qtau = qmax * (err(nql)/toltau)**(2*ONE/7)
  else
    qtau = qmax
    do j = nql,0,-1
      if(err(j) > toltau) exit
      qtau = j*delql
    enddo
  endif

  return

end subroutine atom_kb_tauc_cutoff
