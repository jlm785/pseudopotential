!>  writes the file with data fot later plotting
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.2.0
!>  \date         1980s and 1990s, 30 June 2021, 16 November 2026.
!>  \copyright    GNU Public License v2

subroutine atom_psd_plot(itplot, filename, nr, r, drdi,                  &
              nindx, lcmax, rpsi_ps, rpsi_ae, vpsd, zion,                &
              ifcore, cdv_ae, cdv_ps, cdc_total, cdc,                    &
              tauv_ae, tauv_ps, tauc_total, tauc,                        &
              mxdnr, mxdlc, mxdnw)

!   njtj  ***  plotting routines ***
!   potrw is called to save a usefull number of points
!   of the pseudowave function to make a plot.  The
!   info is written to the current plot.dat file.
!   psd_plot_wt is called to fourier transform the the pseudo
!   wave function and save it to the current plot.dat file.
!   The calls to 1)psd_plot_ran take the fourier transform of
!   the potential and saves it in the current plot.dat file,
!   2)psd_plot_rv saves the potential in the current plot.dat file
!   3)zion is saved to the current plot.dat file wtih a
!   marker 'zio' for latter plotting

! Adapted from the old code, 30 June 2021. JLM
! so->iso, vpsd, rpsi, 15 September 2021. JLM
! indx, ist, 3 November 2021. JLM
! charge and kinetic energy densities (valence, core), without 4 pi r^2. 16 November 2026. claude+JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension  of the number of radial points
  integer, intent(in)               ::  mxdlc                            !<  maximum angular momentum scattering channel
  integer, intent(in)               ::  mxdnw                            !<  dimension of number of wave-functions same l

  integer, intent(in)               ::  itplot                           !<  default plot file
  character(len=*)                  ::  filename                         !<  name of plot file

  integer, intent(in)               ::  nr                               !<  number of radial points
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r / d i

  integer, intent(in)               ::  nindx(0:mxdlc,-1:1)              !<  number of orbitals for scattering channel l
  integer, intent(in)               ::  lcmax                            !<  maximum ang mom in scattering channels

  real(REAL64), intent(in)      ::  rpsi_ps(mxdnr,mxdnw,0:mxdlc,-1:1)    !<  r*pseudo-wave-function
  real(REAL64), intent(in)      ::  rpsi_ae(mxdnr,mxdnw,0:mxdlc,-1:1)    !<  r*pseudo-wave-function

  real(REAL64), intent(in)          ::  vpsd(mxdnr,0:mxdlc,-1:1)         !<  pseudo-potential in Rydberg

  real(REAL64), intent(in)          ::  zion                             !<  pseudo-potential ionic charge

  integer, intent(in)               ::  ifcore                           !<  0: no partial core correction; 1 : xc partial core correction; 2 : not recommended.
  real(REAL64), intent(in)          ::  cdv_ae(mxdnr)                    !<  4 pi r^2 all-electron valence charge density
  real(REAL64), intent(in)          ::  cdv_ps(mxdnr)                    !<  4 pi r^2 pseudo valence charge density
  real(REAL64), intent(in)          ::  cdc_total(mxdnr)                 !<  4 pi r^2 core charge density
  real(REAL64), intent(in)          ::  cdc(mxdnr)                       !<  4 pi r^2 partial core charge density (or total)
  real(REAL64), intent(in)          ::  tauv_ae(mxdnr)                   !<  4 pi r^2 all-electron valence kinetic energy density
  real(REAL64), intent(in)          ::  tauv_ps(mxdnr)                   !<  4 pi r^2 pseudo valence kinetic energy density
  real(REAL64), intent(in)          ::  tauc_total(mxdnr)                !<  core kinetic energy density (without 4 pi r^2)
  real(REAL64), intent(in)          ::  tauc(mxdnr)                      !<  4 pi r^2 partial core kinetic energy density (or total)

! local allocatable arrays

  real(REAL64), allocatable         ::  aux(:)
  real(REAL64), allocatable         ::  den(:)                           !  charge or kinetic energy density (without 4 pi r^2)

! local variables

  integer        ::  ist                               !  sign/scaling
  integer        ::  nrplot
  real(REAL64)   ::  ch                                !  charge up to RMAX
  integer        ::  nrden                             !  last point of the density plots
  real(REAL64)   ::  xint                              !  core charge or kinetic energy

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  RMAX = 15.0_REAL64
  real(REAL64), parameter    ::  QMAX = 40.0_REAL64
  real(REAL64), parameter    ::  QSTEP = 0.1_REAL64
  real(REAL64), parameter    ::  RSTEP = 0.05_REAL64
  real(REAL64), parameter    ::  RMAXD = 10.0_REAL64
  real(REAL64), parameter    ::  RSTEPD = 0.005_REAL64
  real(REAL64), parameter    ::  SMALL = 1.0E-6_REAL64
  real(REAL64), parameter    ::  PI = 4*atan(1.0_REAL64)

! counter

  integer     ::  is, i, j, l


!njtj    ***  plotting routines ***

  open(unit=itplot, file=trim(filename), form='FORMATTED', status='UNKNOWN')

  allocate(aux(nr))

  do j = 1,nr
    nrplot = j
    if(r(j) > RMAX) exit
  enddo

  do is = -1,1
  do l = 0,lcmax
    if(nindx(l,is) > 0) then
      do j = 1,nindx(l,is)

!       Plots only if charge mostly inside RMAX

        do i = 1,nrplot
          aux(i) = rpsi_ps(i,j,l,is)*rpsi_ps(i,j,l,is)
        enddo
        call atom_atm_integ(ch, aux, drdi, nrplot)

        if(ch > 0.2) then

          ist = 1
          if(rpsi_ae(nrplot,j,l,is)*rpsi_ps(nrplot,j,l,is) < 0) ist = -1

          call atom_plot_one(1, nrplot, r, rpsi_ae(:,j,l,is), RSTEP, l,  &
              ist, 'w   ','t   ', .TRUE., itplot, mxdnr)

          call atom_plot_one(1, nrplot, r, rpsi_ps(:,j,l,is), RSTEP, l,  &
              1, 'w   ','p   ', .TRUE., itplot, mxdnr)

          call atom_psd_plot_wt(nr, r, drdi, rpsi_ps(:,j,l,is), QMAX,    &
            QSTEP, l, 1, itplot, mxdnr)

        endif

      enddo

    endif
  enddo
  enddo

  do is = -1,1
  do l = 0,lcmax
    if(nindx(l,is) > 0) then

      do j = 2,nr
        aux(j) = vpsd(j,l,is)/r(j)
      enddo

!njtj    ***  plotting routines ***

      call atom_psd_plot_ran(r, nr, aux, QMAX, QSTEP, l+1, zion,         &
          itplot, mxdnr)

      call atom_plot_one(2, nrplot, r, aux, RSTEP, l, 1,                 &
          'vn  ','    ', .TRUE., itplot, mxdnr)

!njtj    ***  user should adjust for their needs  ***

    endif
  enddo
  enddo

! charge densities and kinetic energy densities, without the 4 pi r^2 factor
! that would hide a cusp at the nucleus.  From the first point after r = 0 with a
! fine step, exponential format (the core kinetic energy density of heavy atoms
! is very large close to the nucleus).
! markers d (charge density) or t (kinetic energy density), v (valence) or c (core),
! a (all-electron) or p (pseudo, or partial core)

  do j = 1,nr
    nrden = j
    if(r(j) > RMAXD) exit
  enddo

  allocate(den(nr))
  den(1) = ZERO

  do j = 2,nr
    den(j) = cdv_ae(j) / (4*PI*r(j)*r(j))
  enddo
  call atom_plot_one_es(2, nrden, r, den, RSTEPD, 'dva', itplot, mxdnr)
  do j = 2,nr
    den(j) = cdv_ps(j) / (4*PI*r(j)*r(j))
  enddo
  call atom_plot_one_es(2, nrden, r, den, RSTEPD, 'dvp', itplot, mxdnr)

  call atom_atm_integ(xint, cdc_total, drdi, nr)
  if(xint > SMALL) then
    do j = 2,nr
      den(j) = cdc_total(j) / (4*PI*r(j)*r(j))
    enddo
    call atom_plot_one_es(2, nrden, r, den, RSTEPD, 'dca', itplot, mxdnr)
    if(ifcore /= 0) then
      do j = 2,nr
        den(j) = cdc(j) / (4*PI*r(j)*r(j))
      enddo
      call atom_plot_one_es(2, nrden, r, den, RSTEPD, 'dcp', itplot, mxdnr)
    endif
  endif

  do j = 2,nr
    den(j) = tauv_ae(j) / (4*PI*r(j)*r(j))
  enddo
  call atom_plot_one_es(2, nrden, r, den, RSTEPD, 'tva', itplot, mxdnr)
  do j = 2,nr
    den(j) = tauv_ps(j) / (4*PI*r(j)*r(j))
  enddo
  call atom_plot_one_es(2, nrden, r, den, RSTEPD, 'tvp', itplot, mxdnr)

  do j = 1,nr
    aux(j) = 4*PI*r(j)*r(j)*tauc_total(j)
  enddo
  call atom_atm_integ(xint, aux, drdi, nr)
  if(xint > SMALL) then
    call atom_plot_one_es(2, nrden, r, tauc_total, RSTEPD, 'tca', itplot, mxdnr)
    do j = 2,nr
      den(j) = tauc(j) / (4*PI*r(j)*r(j))
    enddo
    call atom_plot_one_es(2, nrden, r, den, RSTEPD, 'tcp', itplot, mxdnr)
  endif

  deallocate(den)

! marker before data in contrast to other markers!

  write(itplot,'(1x,"marker zio")')
  write(itplot,'(2x,f5.2)') zion

!njtj    ***  user should adjust for their needs  ***

  close(unit=itplot)

  deallocate(aux)

  return

end subroutine atom_psd_plot
