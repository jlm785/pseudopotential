!>  Calculates the error and jacobian of the equations of the
!>  improved scheme of N. Troullier and J. L. Martins
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.0.6
!>  \date         3 November 2021.
!>  \copyright    GNU Public License v2

subroutine atom_psd_tm2_func(r, drdi, jrc, lo,                           &
      polydrc, cdrc, bkrk,                                               &
      errsq, fvec, xj,                                                   &
      mxdnr)


 implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter          ::  ntm2 = 6                              !  order of tm2 polynomial

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension  of the number of radial points

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  integer, intent(in)               ::  jrc                              !<  r(jrc) = r_cut

  real(REAL64), intent(in)          ::  cdrc                             !<  int_0^r_cut rho(r) dr
  real(REAL64), intent(in)          ::  polydrc(0:ntm2-2)                !<  n-th derivatives of Kerker polynomial ar rc

  integer, intent(in)               ::  lo                               !<  angular quantum number l

  real(REAL64), intent(in)          ::  bkrk(0:ntm2)                     !<  coefficient of polynomial

! output

  real(REAL64), intent(out)         ::  errsq                            !<  sum of square errors
  real(REAL64), intent(out)         ::  fvec(0:ntm2)                     !<  error in the solution of equation
  real(REAL64), intent(out)         ::  xj(0:ntm2,0:ntm2)                !<  jacobian


! local variables

  real(REAL64)         ::  cdps                                          !  int_0^r_cut rho_psd(r) dr

! local allocatable arrays

  real(REAL64), allocatable         ::  rcpn(:)                          !  rc^n
  real(REAL64), allocatable         ::  rpsi(:)                          !  r*psi
  real(REAL64), allocatable         ::  aux(:)

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64

! counter

  integer      ::  j, k


  do j = 0,ntm2
  do k = 0,ntm2
    xj(k,j) = ZERO
  enddo
  enddo

  allocate(rcpn(0:2*ntm2))
  allocate(rpsi(jrc))
  allocate(aux(jrc))

  rcpn(0) = ONE
  rcpn(1) = r(jrc)
  do j = 2,2*ntm2
    rcpn(j) = rcpn(j-1)*r(jrc)
  enddo

  do j = 0,ntm2
    xj(0,j) = ONE
  enddo
  do k = 1,ntm2-2
    do j = (k+1)/2,ntm2
      xj(k,j) = (2*j-k+1)*xj(k-1,j)
    enddo
  enddo
  do k = 0,ntm2-2
    do j = (k+1)/2,ntm2
      xj(k,j) = xj(k,j)*rcpn(2*j-k)
    enddo
  enddo

! error in the polynomial derivatives

  do k = 0,ntm2-2
    fvec(k) = ZERO
    do j = (k+1)/2,ntm2
      fvec(k) = fvec(k) + bkrk(j)*xj(k,j)
    enddo
    fvec(k) = fvec(k) - polydrc(k)
  enddo

! V''(0) = 0

  fvec(ntm2-1) = bkrk(1)*bkrk(1) + (2*lo+5)*bkrk(2)
  do j = 0,ntm2
    xj(ntm2-1,j) = ZERO
  enddo
  xj(ntm2-1,1) = 2*bkrk(1)
  xj(ntm2-1,2) = (2*lo+5)*ONE

! error in norm conservation

  call atom_psd_krk_psi(ntm2, bkrk, .TRUE., jrc, r, lo, rpsi, jrc)

  do k = 1,jrc
    aux(k) = rpsi(k)*rpsi(k)
  enddo

  call atom_atm_integ(cdps, aux, drdi, jrc)

  fvec(ntm2) = cdps - cdrc

  xj(ntm2,0) = 2*cdps

  do j = 1,ntm2

    do k = 1,jrc
      aux(k) = aux(k)*r(k)*r(k)
    enddo

    call atom_atm_integ(cdps, aux, drdi, jrc)

    xj(ntm2,j) = 2*cdps

  enddo

! square modulo of error vector

  errsq = ZERO
  do k = 0,ntm2
    errsq = errsq + fvec(k)*fvec(k)
  enddo

  deallocate(rcpn)
  deallocate(rpsi)
  deallocate(aux)

  return

end subroutine atom_psd_tm2_func
