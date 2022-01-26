!>  writes the file with data fot later plotting
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.0.6
!>  \date         1980s and 1990s, 30 June 2021, 3 November 2021.
!>  \copyright    GNU Public License v2

subroutine atom_psd_plot(itplot, filename, nr, r, drdi,                  &
              nindx, lcmax, rpsi_ps, rpsi_ae, vpsd, zion,                &
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

! local allocatable arrays

  real(REAL64), allocatable         ::  aux(:)

! local variables

  integer        ::  ist                               !  sign/scaling
  integer        ::  nrplot
  real(REAL64)   ::  ch                                !  charge up to RMAX

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  RMAX = 15.0_REAL64
  real(REAL64), parameter    ::  QMAX = 40.0_REAL64
  real(REAL64), parameter    ::  QSTEP = 0.1_REAL64
  real(REAL64), parameter    ::  RSTEP = 0.05_REAL64

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

! marker before data in contrast to other markers!

  write(itplot,'(1x,"marker zio")')
  write(itplot,'(2x,f5.2)') zion

!njtj    ***  user should adjust for their needs  ***

  close(unit=itplot)

  deallocate(aux)

  return

end subroutine atom_psd_plot
