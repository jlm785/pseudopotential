!>  Fills the arrays with the log derivatives for later plotting
!>  for the Kleinman-Bylander type of pseudopotential
!>
!>  \author       Norm Troullier, Manuel Maria Alemany, Jose Luis Martins
!>  \version      6.013
!>  \date         1980s and 1990s, 21 August 2021. 25 September 2021.
!>  \copyright    GNU Public License v2

subroutine atom_plot_ln_kb(npoint, numtotal, ehist, dlogx,               &
       nr, r, drdi, d2rodr,                                              &
       nsc, lo, iso,                                                     &
       lmax_v, vkbproj, vionic, vhxc, vlocal, inorm,                     &
       iowrite, mxdnr, mxdl, mxdsc, mxdpts)


! extracted from lnplot.f and converted fo f90
! mxdsc. 25 September 2021. JLM



  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points
  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum
  integer, intent(in)               ::  mxdsc                            !<  dimension of maximum scattering channels
  integer, intent(in)               ::  mxdpts                           !<  dimension of number of histogram points

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

  integer, intent(in)               ::  lmax_v                           !<  maximum angular momentum in potential
  real(REAL64), intent(in)          ::  vionic(mxdnr,0:mxdl,-1:1)        !<  r*ionic potential in Rydberg

  real(REAL64), intent(in)          ::  vhxc(mxdnr,-1:1)                 !<  screening potential in Rydberg

  real(REAL64), intent(in)          ::  vkbproj(mxdnr,0:mxdl,-1:1)       !<  KB projector
  real(REAL64), intent(in)          ::  vlocal(mxdnr)                    !<  local potential
  integer, intent(in)               ::  inorm(0:mxdl,-1:1)               !<  normalizatio

! output

  real(REAL64), intent(out)         ::  dlogx(mxdpts,mxdsc)              !<  logarithmic derivatives

! local variables

  integer        ::  iflag

  real(REAL64)   ::  prowav

  integer        ::  nrevi
  integer        ::  l, lvio

  logical        ::  ltry                                                !  successful integration
  logical        ::  lsafe                                               !  potential not influenced by viod

  integer        ::  ntrymax, n_error
  real(REAL64)   ::  alpha
  real(REAL64)   ::  err

! allocatable work arrays

  real(REAL64), allocatable         ::  v(:)                             !  work array (potential with KB)
  real(REAL64), allocatable         ::  vs(:)                            !  work array (local potential)
  real(REAL64), allocatable         ::  vref(:)                          !  work array (semi-local potential)
  real(REAL64), allocatable         ::  ag(:), agp(:)                    !  work array (projections)

  real(REAL64), allocatable         ::  rpsi(:), drpsidr(:)              !  wave-functions and derivatives
  real(REAL64), allocatable         ::  rpsinew(:)                       !  new wave-functions

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  EPS = 1.0E-4_REAL64, SMALL = 1.0E-20_REAL64

! counters

  integer     ::  i, j, numev, ntry


! order and number of derivatives for interpolation


  ntrymax = 1000

! allocations and initialization

  allocate(v(mxdnr))
  allocate(vs(mxdnr))
  allocate(vref(mxdnr))

  allocate(rpsi(mxdnr), drpsidr(mxdnr))
  allocate(rpsinew(mxdnr))

  allocate(ag(mxdnr), agp(mxdnr))

  do i = 1,nsc
    do numev = 1,mxdpts
      dlogx(numev,i) = ZERO
   enddo
  enddo

  do i = 1,nsc

    l = lo(i)
    lvio = min(l,lmax_v)
    do j = 2,npoint
      vs(j) = vlocal(j) + vhxc(j,iso(i)) + l*(l+1)/(r(j)*r(j))
      v(j) = vs(j)
      vref(j) = vionic(j,lvio,iso(i))/r(j) + vhxc(j,iso(i))
    enddo

    if(inorm(lvio,iso(i)) /= 0) then

      ag(1) = ZERO
      do j = 2,npoint
        ag(j) = vkbproj(j,lvio,iso(i))*drdi(j)*r(j)
      enddo

      ltry = .FALSE.

      do numev = 1,numtotal

        if(.not. ltry) then

          call atom_atm_difnrl_one(nr, r, drdi, d2rodr, v, rpsi, drpsidr,   &
              lo(i), ehist(numev), r(npoint), nrevi, iflag,                 &
              iowrite, mxdnr)

          if(nrevi == npoint .and. iflag == 0) then
          else
            do j = nrevi,npoint
              rpsi(j) = SMALL
            enddo
          endif
        endif

        n_error = 0
        alpha = ONE/2

        do ntry = 1,ntrymax

          prowav = ZERO
!         SHOULD FIX THIS
          do j = 1,npoint-4,4
            prowav = prowav +  7*(ag(j  )*rpsi(j  )+ag(j+4)*rpsi(j+4))       &
                            + 32*(ag(j+1)*rpsi(j+1)+ag(j+3)*rpsi(j+3))       &
                            + 12* ag(j+2)*rpsi(j+2)
          enddo
          prowav = 2*prowav/(45*ONE)

          do j = 2,npoint
            agp(j) = inorm(lvio,iso(i))*prowav*vkbproj(j,lvio,iso(i))*r(j)/rpsi(j)
          enddo
          do j = 2,npoint
            v(j) = vs(j) + agp(j)
!           tries to avoid instabilities
            lsafe = .TRUE.
            if(v(j) > vref(j) + 10*ONE) then
              v(j) = vref(j) + 10*ONE
              lsafe = .FALSE.
            elseif(v(j) < vref(j) - 10*ONE) then
              v(j) = vref(j) - 10*ONE
              lsafe = .FALSE.
            endif
          enddo

          call atom_atm_difnrl_one(nr, r, drdi, d2rodr, v, rpsinew, drpsidr,    &
             lo(i), ehist(numev), r(npoint), nrevi, iflag,                      &
             iowrite, mxdnr)

          if(nrevi == npoint .and. iflag == 0 .and. lsafe) then
            dlogx(numev,i) = drpsidr(npoint)/rpsinew(npoint) - ONE/r(npoint)
            ltry = .TRUE.
          else
            dlogx(numev,i) = ZERO
            do j = nrevi+1,npoint
              rpsinew(j) = ZERO
            enddo
            ltry = .FALSE.
          endif

          err = ZERO
          do j = 2,npoint
            err = max(abs(rpsinew(j)-rpsi(j)),err)
          enddo

          if (err < EPS) then

            exit

          else

            do j = 2,npoint
              rpsi(j) = (alpha*rpsinew(j) + (ONE-alpha)*rpsi(j)) + SMALL
            enddo

!           old code was inconsistent

            n_error = n_error + 1
            if (n_error == 10) then
              n_error = 0
              alpha = 0.9*alpha

              if (alpha <= ZERO) then
                write(6,*)
                write(6,*) '  Stopped in atom_plot_ln_kb:'
                write(6,*) '  can not reach self-consistency at energy:',  &
                         ehist(i),' (Ry)'

                 STOP

              endif
            endif
          endif

        enddo

      enddo

    else

      do numev = 1,numtotal

        call atom_atm_difnrl_one(nr, r, drdi, d2rodr, vs, rpsi, drpsidr,   &
            lo(i), ehist(numev), r(npoint), nrevi, iflag,                  &
            iowrite, mxdnr)

        if(nrevi == npoint .and. iflag == 0) then
          dlogx(numev,i) = drpsidr(npoint)/rpsi(npoint) - ONE/r(npoint)
        else
          dlogx(numev,i) = ZERO
        endif

      enddo

    endif

  enddo


  deallocate(v)
  deallocate(vs)
  deallocate(vref)

  deallocate(rpsi, drpsidr)
  deallocate(rpsinew)

  deallocate(ag, agp)

  return

end subroutine atom_plot_ln_kb
