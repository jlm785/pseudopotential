!>  Perturbative polarization orbital (SIESTA recipe):  solves the
!>  inhomogeneous radial equation
!>
!>     - u'' + (veff - ev) u = - r**POLPOW u_src ,   u(0) = 0,  u(r(nrc)) = 0
!>
!>  where u_src = r R_{l-1} is a confined orbital of angular momentum l-1
!>  with eigenvalue ev, and veff (Rydberg) is the effective potential for
!>  angular momentum l, including l(l+1)/r**2 and the same confinement
!>  as u_src.  The result is normalized and positive near the origin.
!>  On exit ev_pol = <u|H|u>.
!>
!>  POLPOW = 1 is the first order change of u_src in an electric field
!>  (l-1 -> l component), the equation of the SIESTA papers.  The SIESTA
!>  code (5.4.2, atom.F, subroutine polarization) integrates the source
!>  -r u_src in the grid index without the factor (dr/di)**(3/2), which
!>  on its logarithmic grid (dr/di proportional to r + b) is the source
!>  - r**(-1/2) u_src, POLPOW = -1/2, a more compact orbital.  For Ge
!>  (cpw DZP, 40 Ha) POLPOW = 1, 1/2, 0, -1/2, -1 gave E_AO - E_PW =
!>  103, 88, 79, 74, 73 meV/atom.  POLPOW = -1 makes u ~ r**(l+1) ln r.
!>
!>  The equation is linear and solved directly (no shooting):  Numerov in
!>  the grid index for y = u / sqrt(dr/di), outward from the origin and
!>  inward from r(nrc), matched at two consecutive points near RMFAC r(nrc).
!>  Assumes the logarithmic grid r(i) = a (exp(b(i-1)) - 1), d2rodr = b.
!>
!>  E. Artacho, D. Sanchez-Portal, P. Ordejon, A. Garcia, and J. M. Soler,
!>  Phys. Stat. Sol. (b) 215, 809 (1999).
!>
!>  \author       claude+J.L.Martins
!>  \version      6.1.0
!>  \date         5 November 2026.
!>  \copyright    GNU Public License v2

subroutine atom_kb_polarization(nr, r, drdi, d2rodr, veff, ev, rpsi_src, &
         nrc, rpsi_pol, ev_pol, iflag,                                   &
         mxdnr)

! written 5 November 2026. claude+JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points

  integer, intent(in)               ::  nr                               !<  number of points in the radial grid
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points   r(i) = a*(exp(b*(i-1))-1), i=1,...,nr
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i
  real(REAL64), intent(in)          ::  d2rodr(mxdnr)                    !<  (d^2 r / d i^2) / (d r / d i)

  real(REAL64), intent(in)          ::  veff(mxdnr)                      !<  effective potential for l (Ry), veff(1) not used
  real(REAL64), intent(in)          ::  ev                               !<  eigenvalue of the source orbital (Ry)
  real(REAL64), intent(in)          ::  rpsi_src(mxdnr)                  !<  source orbital u_src = r R_{l-1}
  integer, intent(in)               ::  nrc                              !<  u(r(nrc)) = 0, u = 0 beyond

! output

  real(REAL64), intent(out)         ::  rpsi_pol(mxdnr)                  !<  polarization orbital u = r R_l, normalized
  real(REAL64), intent(out)         ::  ev_pol                           !<  <u|H|u> (Ry)
  integer, intent(out)              ::  iflag                            !<  0: success; 1: bad nrc; 2: singular matching

! local allocatable arrays

  real(REAL64), allocatable         ::  t(:)                             !  (dr/di)**2 (veff - ev) + b**2/4
  real(REAL64), allocatable         ::  s(:)                             !  (dr/di)**(3/2) r u_src
  real(REAL64), allocatable         ::  yh(:), yp(:)                     !  outward homogeneous and particular solutions
  real(REAL64), allocatable         ::  zh(:), zp(:)                     !  inward homogeneous and particular solutions

! local variables

  integer                           ::  nm                               !  matching at nm and nm+1
  real(REAL64)                      ::  det, alpha, beta, d1, d2
  real(REAL64)                      ::  unorm, umax, ursrc
  integer                           ::  ll

! constants

  real(REAL64), parameter           ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter           ::  RMFAC = 0.9_REAL64               !  matching radius / r(nrc) (start of the confinement)
  real(REAL64), parameter           ::  POLPOW = -0.5_REAL64             !  power of r in the source (1: electric field, -1/2: SIESTA code)
  real(REAL64), parameter           ::  EPS = 1.0E-12_REAL64

! counters

  integer                           ::  j


  iflag = 0
  rpsi_pol(:) = ZERO
  ev_pol = ev

  nm = 0
  do j = 2,nrc
    if(r(j) <= RMFAC*r(nrc)) nm = j
  enddo
  if(nrc > nr .or. nm < 10 .or. nm > nrc-4) then
    iflag = 1

    return

  endif

  allocate(t(nrc), s(nrc))
  allocate(yh(nrc), yp(nrc), zh(nrc), zp(nrc))

  t(1) = ZERO
  s(1) = ZERO
  do j = 2,nrc
    t(j) = drdi(j)*drdi(j)*(veff(j) - ev) + d2rodr(j)*d2rodr(j)/4
    s(j) = drdi(j)*sqrt(drdi(j)) * r(j)**POLPOW * rpsi_src(j)
  enddo
  s(nrc) = ZERO

! outward integration (u(0) = 0, regular solution ~ r**(l+1))

  yh(1) = ZERO
  yh(2) = UM
  yp(1) = ZERO
  yp(2) = ZERO
  do j = 2,nm
    yh(j+1) = (2*yh(j)*(UM + 5*t(j)/12) - yh(j-1)*(UM - t(j-1)/12))       &
              / (UM - t(j+1)/12)
    yp(j+1) = (2*yp(j)*(UM + 5*t(j)/12) - yp(j-1)*(UM - t(j-1)/12)        &
              + (s(j+1) + 10*s(j) + s(j-1))/12) / (UM - t(j+1)/12)
  enddo

! inward integration (u(r(nrc)) = 0)

  zh(nrc) = ZERO
  zh(nrc-1) = UM
  zp(nrc) = ZERO
  zp(nrc-1) = ZERO
  do j = nrc-1,nm+1,-1
    zh(j-1) = (2*zh(j)*(UM + 5*t(j)/12) - zh(j+1)*(UM - t(j+1)/12))       &
              / (UM - t(j-1)/12)
    zp(j-1) = (2*zp(j)*(UM + 5*t(j)/12) - zp(j+1)*(UM - t(j+1)/12)        &
              + (s(j+1) + 10*s(j) + s(j-1))/12) / (UM - t(j-1)/12)
  enddo

! yp + alpha yh = zp + beta zh at nm and nm+1

  d1 = zp(nm) - yp(nm)
  d2 = zp(nm+1) - yp(nm+1)
  det = zh(nm)*yh(nm+1) - yh(nm)*zh(nm+1)
  if(abs(det) < EPS*abs(yh(nm+1)*zh(nm))) then
    iflag = 2
    deallocate(t, s, yh, yp, zh, zp)

    return

  endif
  alpha = (zh(nm)*d2 - zh(nm+1)*d1) / det
  beta = (yh(nm)*d2 - yh(nm+1)*d1) / det

  do j = 1,nm
    rpsi_pol(j) = (yp(j) + alpha*yh(j)) * sqrt(drdi(j))
  enddo
  do j = nm+1,nrc
    rpsi_pol(j) = (zp(j) + beta*zh(j)) * sqrt(drdi(j))
  enddo
  rpsi_pol(1) = ZERO
  rpsi_pol(nrc) = ZERO

! norm and <u|H - ev|u> = - <u|r**POLPOW|u_src> (Simpson, as in atom_kb_orbitals)

  unorm = ZERO
  ursrc = ZERO
  ll = 4
  do j = 2,nrc
    unorm = unorm + ll*rpsi_pol(j)*rpsi_pol(j)*drdi(j)
    ursrc = ursrc + ll*rpsi_pol(j)*r(j)**POLPOW*rpsi_src(j)*drdi(j)
    ll = 6 - ll
  enddo
  unorm = unorm / 3
  ursrc = ursrc / 3

  ev_pol = ev - ursrc / unorm

! normalizes, positive near the origin

  umax = ZERO
  do j = 2,nrc
    umax = max(umax, abs(rpsi_pol(j)))
  enddo
  do j = 2,nrc
    if(abs(rpsi_pol(j)) > umax/100) exit
  enddo
  unorm = sqrt(unorm)
  if(rpsi_pol(j) < ZERO) unorm = -unorm
  do j = 1,nrc
    rpsi_pol(j) = rpsi_pol(j) / unorm
  enddo

  deallocate(t, s, yh, yp, zh, zp)

  return

end subroutine atom_kb_polarization
