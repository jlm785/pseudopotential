!>  fills the arrays with the log derivatives of the true (all-elecron)
!>  or semi-local (pseudo-) potential
!>
!>  \author       Norm Troullier, Manuel Maria Alemany, Jose Luis Martins
!>  \version      6.0.3
!>  \date         1980s and 1990s, 21 August 2021, 26 September 2021.
!>  \copyright    GNU Public License v2

subroutine atom_plot_ln_nokb(ispp, npoint, numtotal, ehist, dlogx,       &
       nr, r, drdi, d2rodr,                                              &
       nsc, lo, iso, znuc,                                               &
       lmax_v, vionic, vhxc,                                             &
       iowrite, mxdnr, mxdl, mxdsc, mxdpts)


! extracted from lnplot.f and converted to f90. JLM
! mxdsc, etc. 25 September 2021. JLM



  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points
  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum
  integer, intent(in)               ::  mxdsc                            !<  dimension of maximum scattering channels
  integer, intent(in)               ::  mxdpts                           !<  dimension of number of histogram points

  character(len=1), intent(in)      ::  ispp                             !<  spin polarization  ' ', 's', 'r'
  integer, intent(in)               ::  npoint                           !<  log derivatives calculated at r(npoint)

  integer, intent(in)               ::  numtotal                         !<  number of energy in histogram
  real(REAL64), intent(in)          ::  ehist(mxdpts)                    !<  energy in histogram

  integer, intent(in)               ::  nr                               !<  dimension of the number of radial points
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i
  real(REAL64), intent(in)          ::  d2rodr(mxdnr)                    !<  (d^2 r(i) / d i^2) / (d r / d i)

  integer, intent(in)               ::  nsc                              !<  number of scattering channels
  integer, intent(in)               ::  lo(mxdsc)                        !<  angular quantum number l
  integer, intent(in)               ::  iso(mxdsc)                       !<  2*spin or 2*(j-l)

  real(REAL64), intent(in)          ::  znuc                             !<  nuclear charge

  integer, intent(in)               ::  lmax_v                           !<  maximum angular momentum in potential
  real(REAL64), intent(in)          ::  vionic(mxdnr,0:mxdl,-1:1)        !<  r*ionic potential in Rydberg

  real(REAL64), intent(in)          ::  vhxc(mxdnr,-1:1)                 !<  screening potential in Rydberg

! output

  real(REAL64), intent(out)         ::  dlogx(mxdpts,mxdsc)              !<  logarithmic derivatives

! local variables

  integer        ::  iflag

  real(REAL64)   ::  vzero
  real(REAL64)   ::  dtmp

  integer        ::  nord, nd, ns

  integer        ::  nrevi
  integer        ::  l, lvio

! allocatable work arrays

  real(REAL64), allocatable         ::  v(:)                             !  work array (potential)
  real(REAL64), allocatable         ::  ar(:), br(:)                     !  work array (wave-functions or derivatives)

  real(REAL64), allocatable         ::  xin(:), yin(:)
  real(REAL64), allocatable         ::  xd(:), dxd(:)

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64

! counters

  integer     ::  i, j, numev


! order and number of derivatives for interpolation

  nord = 6
  nd = 1

  ns = (nord+3)/2

! allocations and initialization

  allocate(v(mxdnr))
  allocate(ar(mxdnr), br(mxdnr))

  allocate(xin(nord+1),yin(nord+1))
  allocate(xd(nd+1),dxd(nd+1))

  do i = 1,nsc
    do numev = 1,mxdpts
      dlogx(numev,i) = ZERO
   enddo
  enddo

  do i = 1,nsc

    l = lo(i)
    lvio = min(l, lmax_v)
    if(ispp == 'r' ) then
      do j = 2,nr
        v(j) = vionic(j,lvio,iso(i))/r(j) + vhxc(j,iso(i))
      enddo
      vzero = vhxc(1,iso(i))
    else
      do j = 2,nr
        v(j) = vionic(j,lvio,iso(i))/r(j) + vhxc(j,iso(i)) + l*(l+1)/(r(j)*r(j))
      enddo
    endif


    do numev = 1,numtotal

      if (ispp == 'r') then

        call atom_atm_difrel_one(nr, r, drdi, v, ar, br,                 &
            lo(i), iso(i), znuc, vzero,                                  &
            ehist(numev), r(npoint+nord-ns+1), nrevi, iflag,             &
            iowrite, mxdnr)

        if(nrevi == npoint+nord-ns+1 .and. iflag == 0) then
          do j = 1,nord+1
            xin(j) = r(npoint - ns + j) - r(npoint)
            yin(j) = (ar(npoint - ns + j)*ar(npoint - ns + j) +          &
                       br(npoint - ns + j)*br(npoint - ns + j))
            yin(j) = sqrt(yin(j))/r(npoint - ns + j)
          enddo
          call poly_interp(xd, dxd, xin, yin, nord, nd)
          dtmp = xd(2)/xd(1)
        else
          dtmp = ZERO
        endif

        dlogx(numev,i) = dtmp

      else

        call atom_atm_difnrl_one(nr, r, drdi, d2rodr, v, ar, br,         &
            lo(i), ehist(numev), r(npoint), nrevi, iflag,                &
            iowrite, mxdnr)

        if(nrevi == npoint .and. iflag == 0) then
          dtmp = br(npoint)/ar(npoint) - ONE/r(npoint)
        else
          dtmp = ZERO
        endif

        dlogx(numev,i) = dtmp

      endif

    enddo

  enddo


  deallocate(v)
  deallocate(ar, br)

  deallocate(xin,yin)
  deallocate(xd,dxd)

  return

end subroutine atom_plot_ln_nokb
