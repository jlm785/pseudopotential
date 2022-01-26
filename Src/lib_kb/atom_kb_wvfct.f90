!>  Calculates the pseudo-wavefunctions
!>
!>  \author       N. Troullier, J.L.Martins
!>  \version      6.0.4
!>  \date         November 1990, May 2012, July 2015, 13 October 2021.
!>  \copyright    GNU Public License v2

subroutine atom_kb_wvfct(npot, lo, irel, nr, r, drdi, d2rodr,            &
        vionic, vscreen, ev, rpsi,                                       &
        iowrite, mxdl, mxdnr)

! mxdl, iowrite, 20 September 2021. JLM
! atom_atm_trieig.  13 October 2021. JLM


  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum
  integer, intent(in)               ::  mxdnr                            !<  dimension of radial grid points.

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  npot(-1:1)                       !<  number of orbitals to be calculated
  integer, intent(in)               ::  lo(mxdl+1,-1:1)                  !<  angular momentum of orbital
  character(len=3), intent(in)      ::  irel                             !<  flag for relativistic (r) or spin-polarized (s) original calculation

  integer, intent(in)               ::  nr                               !<  number of points in the radial grid
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points   r(i) = a*(exp(b*i)-1), i=1,...,nr
  real(REAL64), intent(in)          ::  drdi (mxdnr)                     !<  (d r(i) / d i)
  real(REAL64), intent(in)          ::  d2rodr(mxdnr)                    !<  (d^2 r(i) / d i^2) / (d r / d i)
  real(REAL64), intent(in)          ::  vionic(mxdnr,0:mxdl,-1:1)        !<  r*pseudopotential, -1: j=l-1/2;  0: average;  1:  j=l+1/2
  real(REAL64), intent(in)          ::  vscreen(mxdnr)                   !<  screening potential

! output

  real(REAL64), intent(out)         ::  ev(0:mxdl,-1:1)                  !<  eigenvalues (l,2j-2l)
  real(REAL64), intent(out)         ::  rpsi(mxdnr,0:mxdl,-1:1)          !<  r*wavefunctions (r(i),l,2j-2l).

! allocatable local variables

  real(REAL64), allocatable         ::  y(:), yp(:)                      !  for atom_atm_trieig
  real(REAL64), allocatable         ::  z(:,:)                           !  for atom_atm_trieig
  real(REAL64), allocatable         ::  v(:)                             !  effective potential
  real(REAL64), allocatable         ::  drpsidr(:)                       !  d r*psi / dr


! local variables

  real(REAL64)             ::  c1, c2
  integer                  ::  l, lp, llp, jmax, iflag

  real(REAL64)             ::  e(1)                                      !  only one eigenvalue needed

  integer                  ::  nrm

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  TOL = 1.0E-12_REAL64

! counters

  integer               ::  i, j, k


  allocate(y(mxdnr),yp(mxdnr))
  allocate(z(mxdnr,1))
  allocate(v(mxdnr))
  allocate(drpsidr(mxdnr))


!  up potentials

  c2 = -ONE/d2rodr(1)**2
  c1 = -2*c2 + ONE/4

! check this out     WARNING
  C1 = C1*D2RODR(1)**2
  C2 = C2*D2RODR(1)**2

  jmax = 0
  if(irel == 'rel') jmax = 1

  do j = -jmax,jmax
    do i = 1,npot(j)
      l = lo(i,j)
      lp = l+1
      llp = lp*(lp-1)
!     set up hamiltonian matrix for kinetic energy,
!     only the diagonal depends on the potential.

      y(1)  = c1 / drdi(2)**2
      yp(1)  = ZERO
      do k = 3,nr
        y(k-1)  = c1 / drdi(k)**2
        yp(k-1)  = c2 / (drdi(k)*drdi(k-1))
      enddo

      v(1) = ZERO
      do k = 2,nr
        v(k) = (vionic(k,l,j)  + llp/r(k)) / r(k)  + vscreen(k)
      enddo

!     add the potential

      do k = 2,nr
        y(k-1) = y(k-1) + (vionic(k,l,j)+llp/r(k)) / r(k) + vscreen(k)
      enddo

!     diagonalize and find wave function and store in rpsi().

      nrm = nr-1

      call atom_atm_trieig(nrm, y, yp, 1, e, z, iowrite, mxdnr)

      ev(l,j) =  e(1)
      call atom_atm_difnrl(nr, r, drdi, d2rodr, v, rpsi(:,l,j), drpsidr,      &
          l+1, l, ev(l,j), iflag, TOL,                                        &
          iowrite, mxdnr)
    enddo
  enddo

  deallocate(y, yp)
  deallocate(z)
  deallocate(v)
  deallocate(drpsidr)

  return

end subroutine atom_kb_wvfct
