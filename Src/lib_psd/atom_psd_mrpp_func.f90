!>  Calculates the error and jacobian of the equations of the
!>  improved scheme of N. Troullier and J. L. Martins
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.0.6
!>  \date         3 November 2021.
!>  \copyright    GNU Public License v2

subroutine atom_psd_mrpp_func(r, drdi, d2rodr, jrc, lo,                  &
      polydrc, cdrc, cdrc2, ar2, ar2p, ev1, ev2, bkrk,                   &
      errsq, fvec, xj,                                                   &
      iowrite, mxdnr)


 implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter          ::  MRPP = 8                              !  order of tm2 polynomial

! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdnr                            !<  dimension  of the number of radial points

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i
  real(REAL64), intent(in)          ::  d2rodr(mxdnr)                    !<  (d^2 r(i) / d i^2) / (d r / d i)

  integer, intent(in)               ::  jrc                              !<  r(jrc) = r_cut

  real(REAL64), intent(in)          ::  polydrc(0:MRPP-4)                !<  n-th derivatives of Kerker polynomial ar rc
  real(REAL64), intent(in)          ::  cdrc                             !<  int_0^r_cut rho(r) dr
  real(REAL64), intent(in)          ::  cdrc2                            !<  int_0^r_cut rho(r) dr second orbital

  real(REAL64), intent(in)          ::  ar2                              !<  r*psi(jrc) for second orbital
  real(REAL64), intent(in)          ::  ar2p                             !<  d r*psi / d r (jrc) for second orbital

  real(REAL64), intent(in)          ::  ev1                              !<  orbital energy of main orbital
  real(REAL64), intent(in)          ::  ev2                              !<  orbital energy of second orbital

  integer, intent(in)               ::  lo                               !<  angular quantum number l

  real(REAL64), intent(in)          ::  bkrk(0:MRPP)                     !<  coefficient of polynomial

! output

  real(REAL64), intent(out)         ::  errsq                            !<  sum of square errors
  real(REAL64), intent(out)         ::  fvec(0:MRPP)                     !<  error in the solution of equation
  real(REAL64), intent(out)         ::  xj(0:MRPP,0:MRPP)                !<  jacobian


! local variables

  real(REAL64)         ::  cdps                                          !  int_0^r_cut rho_psd(r) dr
  real(REAL64)         ::  cdps2                                         !  int_0^r_cut rho_psd(r) dr  second orbital

  integer              ::  nrevi, iflag

  real(REAL64)         ::  scal
  real(REAL64)         ::  hstep

  real(REAL64)         ::  bkrk_p(0:MRPP), bkrk_m(0:MRPP)

  real(REAL64)         ::  br_p, br_m
  real(REAL64)         ::  cdps2_p, cdps2_m

! local allocatable arrays

  real(REAL64), allocatable         ::  rcpn(:)                          !  rc^n
  real(REAL64), allocatable         ::  rpsi(:)                          !  r*psi
  real(REAL64), allocatable         ::  aux(:)

  real(REAL64), allocatable         ::  vj(:)
  real(REAL64), allocatable         ::  ar(:)
  real(REAL64), allocatable         ::  br(:)


! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64

! counter

  integer      ::  j, k


  do j = 0,MRPP
  do k = 0,MRPP
    xj(k,j) = ZERO
  enddo
  enddo

  allocate(rcpn(0:2*MRPP))
  allocate(rpsi(jrc))
  allocate(aux(jrc))

  allocate(vj(jrc))
  allocate(ar(jrc))
  allocate(br(jrc))

  rcpn(0) = ONE
  rcpn(1) = r(jrc)
  do j = 2,2*MRPP
    rcpn(j) = rcpn(j-1)*r(jrc)
  enddo

  do j = 0,MRPP
    xj(0,j) = ONE
  enddo
  do k = 1,MRPP-4
    do j = (k+1)/2,MRPP
      xj(k,j) = (2*j-k+1)*xj(k-1,j)
    enddo
  enddo
  do k = 0,MRPP-4
    do j = (k+1)/2,MRPP
      xj(k,j) = xj(k,j)*rcpn(2*j-k)
    enddo
  enddo

! error in the polynomial derivatives

  do k = 0,MRPP-4
    fvec(k) = ZERO
    do j = (k+1)/2,MRPP
      fvec(k) = fvec(k) + bkrk(j)*xj(k,j)
    enddo
    fvec(k) = fvec(k) - polydrc(k)
  enddo

! V''(0) = 0

  fvec(MRPP-3) = bkrk(1)*bkrk(1) + (2*lo+5)*bkrk(2)
  do j = 0,MRPP
    xj(MRPP-3,j) = ZERO
  enddo
  xj(MRPP-3,1) = 2*bkrk(1)
  xj(MRPP-3,2) = (2*lo+5)*ONE

! error in norm conservation

  call atom_psd_krk_psi(MRPP, bkrk, .TRUE., jrc, r, lo, rpsi, jrc)

  do k = 1,jrc
    aux(k) = rpsi(k)*rpsi(k)
  enddo

  call atom_atm_integ(cdps, aux, drdi, jrc)

  fvec(MRPP-2) = cdps - cdrc

  xj(MRPP-2,0) = 2*cdps

  do j = 1,MRPP

    do k = 1,jrc
      aux(k) = aux(k)*r(k)*r(k)
    enddo

    call atom_atm_integ(cdps, aux, drdi, jrc)

    xj(MRPP-2,j) = 2*cdps

  enddo

! second orbital

! effective potential

  call atom_psd_krk_v(MRPP, bkrk, .TRUE., jrc, r, lo, ev1, vj, mxdnr)

  if(lo > 0) then
    do k = 2,jrc
      vj(k) = vj(k) + (lo*(lo+1)) / (r(k)*r(k))
    enddo
  endif

! wave-function

  call atom_atm_difnrl_one(jrc, r, drdi, d2rodr, vj, ar, br,             &
      lo, ev2, r(jrc), nrevi, iflag,                                     &
      iowrite, mxdnr)

  scal = ar2 / ar(jrc)
  do k = 1,jrc
    ar(k) = ar(k)*scal
    br(k) = br(k)*scal
  enddo

  do k = 1,jrc
    aux(k) = ar(k)*ar(k)
  enddo

  call atom_atm_integ(cdps2, aux, drdi, jrc)

  fvec(MRPP-1) = br(jrc) - ar2p
  fvec(MRPP) = cdps2 - cdrc2

! jacobian by finite differences

  xj(MRPP-1,0) = ZERO
  xj(MRPP,0) = ZERO

!   WRITE(6,*) '  VJ(JRC) REFERENCE = ',VJ(JRC)

  do j = 1,MRPP

    hstep = 0.0001 / (rcpn(j)*j*j)

    do k = 0,MRPP
      bkrk_p(k) = bkrk(k)
      bkrk_m(k) = bkrk(k)
    enddo
    bkrk_p(j) = bkrk(j) + hstep
    bkrk_m(j) = bkrk(j) - hstep

    call atom_psd_krk_v(MRPP, bkrk_p, .TRUE., jrc, r, lo, ev1, vj, mxdnr)

    if(lo > 0) then
      do k = 2,jrc
        vj(k) = vj(k) + (lo*(lo+1)) / (r(k)*r(k))
      enddo
    endif

!     WRITE(6,*) '  J = ',J
!     WRITE(6,*) '  VJ(JRC) PLUS = ',VJ(JRC)

    call atom_atm_difnrl_one(jrc, r, drdi, d2rodr, vj, ar, br,           &
      lo, ev2, r(jrc), nrevi, iflag,                                     &
      iowrite, mxdnr)

    scal = ar2 / ar(jrc)
    do k = 1,jrc
      ar(k) = ar(k)*scal
      br(k) = br(k)*scal
    enddo

    do k = 1,jrc
      aux(k) = ar(k)*ar(k)
    enddo

    call atom_atm_integ(cdps2_p, aux, drdi, jrc)

    br_p = br(jrc)

    call atom_psd_krk_v(MRPP, bkrk_m, .TRUE., jrc, r, lo, ev1, vj, mxdnr)

    if(lo > 0) then
      do k = 2,jrc
        vj(k) = vj(k) + (lo*(lo+1)) / (r(k)*r(k))
      enddo
    endif

!     WRITE(6,*) '  VJ(JRC) MINUS = ',VJ(JRC)

    call atom_atm_difnrl_one(jrc, r, drdi, d2rodr, vj, ar, br,           &
      lo, ev2, r(jrc), nrevi, iflag,                                     &
      iowrite, mxdnr)

    scal = ar2 / ar(jrc)
    do k = 1,jrc
      ar(k) = ar(k)*scal
      br(k) = br(k)*scal
    enddo

    do k = 1,jrc
      aux(k) = ar(k)*ar(k)
    enddo

    call atom_atm_integ(cdps2_m, aux, drdi, jrc)

    br_m = br(jrc)

    xj(MRPP-1,j) = (br_p - br_m) / (2*hstep)

    xj(MRPP,j) = (cdps2_p - cdps2_m) / (2*hstep)

  enddo

! square modulo of error vector

  errsq = ZERO
  do k = 0,MRPP
    errsq = errsq + fvec(k)*fvec(k)
  enddo

!   WRITE(6,'(20X,"  ERRO = ",E18.5)') ERRSQ

  deallocate(rcpn)
  deallocate(rpsi)
  deallocate(aux)

  deallocate(vj)
  deallocate(ar)
  deallocate(br)

  return

end subroutine atom_psd_mrpp_func
