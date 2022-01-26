!>  writes the file with data fot later plotting with plot_kb_show
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.0.4
!>  \date         4 September 2021. 6 October 2021.
!>  \copyright    GNU Public License v2

subroutine atom_kb_plot(itplot, filename, irel, nr, r, vlocal,           &
                    rpsi_ps, rpsi_b, nqbas, delqbas, ektot,              &
                    npot, lo, norbas, lo_b,                              &
                    inorm, vkbproj, nqnl, delql, vkbprft,                &
                    mxdnr, nqmax, mxdl, mxdbas)

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

! kinetic (ektot). 6 October 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimenasion of radial grid points.
  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum
  integer, intent(in)               ::  nqmax                            !<  dimension of fourier grid points for potential
  integer, intent(in)               ::  mxdbas                           !<  maximum number of basis functions

  integer, intent(in)               ::  itplot                           !<  default plot file
  character(len=*)                  ::  filename                         !<  name of plot file

  character(len=3), intent(in)      ::  irel                             !<  flag for relativistic (r) or spin-polarized (s) original calculation

  integer, intent(in)               ::  nr                               !<  number of radial points
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid

  real(REAL64), intent(in)          ::  vlocal(mxdnr)                    !<  r*local pseudopotential

  real(REAL64), intent(in)          ::  rpsi_ps(mxdnr,0:mxdl)            !<  r*pseudo-wave-function
  real(REAL64), intent(in)          ::  rpsi_b(mxdnr,mxdbas)             !<  r*basis-function

  integer, intent(in)               ::  nqbas                            !<  number of points in the fourier grid for basis
  real(REAL64), intent(in)          ::  delqbas                          !<  spacing of the fourier grid for basis
  real(REAL64), intent(in)          ::  ektot(0:nqmax)                   !<  Cumulative kineti energy of wave-functions in Fourier space

  integer, intent(in)               ::  npot(-1:1)                       !<  number of orbitals to be calculated
  integer, intent(in)               ::  lo(mxdl+1,-1:1)                  !<  angular momentum of orbital
  integer, intent(in)               ::  norbas                           !<  number of basis functions
  integer, intent(in)               ::  lo_b(mxdbas)                     !<  angular momentum of basis function

  integer, intent(in)               ::  inorm(0:mxdl,-1:1)               !<  sign of denominator of KB operator
  real(REAL64), intent(in)          ::  vkbproj(mxdnr,0:mxdl,-1:1)       !<  kb-projector

  integer, intent(in)               ::  nqnl                             !<  number of points for the non-local Fourier grid
  real(REAL64), intent(in)          ::  delql                            !<  spacing of points in Fourier grid
  real(REAL64), intent(in)          ::  vkbprft(0:nqmax,0:mxdl,-1:1)     !<  Fourier transform of kb-projector

! local variables

  integer        ::  nrplot, nrp_v, nqp
  integer        ::  l                                 !  angular momentum

  real(REAL64), allocatable         ::  q(:)                             !  wave-vector
  real(REAL64), allocatable         ::  ekdif(:)                         !  kinetic energy difference
  real(REAL64), allocatable         ::  vk(:)                            !  local potential

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  SMALL = 1.0E-8_REAL64

  real(REAL64), parameter    ::  RMAX = 10.0_REAL64
  real(REAL64), parameter    ::  QMAX = 12.0_REAL64
  real(REAL64), parameter    ::  QSTEP = 0.25_REAL64
  real(REAL64), parameter    ::  RSTEP = 0.05_REAL64


! counter

  integer     ::  i, j


  open(unit=itplot, file=trim(filename), form='FORMATTED', status='UNKNOWN')

! plot local potential

  do j = 1,nr
    nrplot = j
    if(r(j) > RMAX) exit
  enddo

  allocate(vk(nrplot))

  do j = 2,nrplot
    vk(j) = vlocal(j) / r(j)
  enddo

  call atom_plot_one(2, nrplot, r, vk, RSTEP/2, -1, 1,                   &
           'loca','l   ', .TRUE., itplot, mxdnr)

  deallocate(vk)

! plot basis versus wave-function

  do i = 1,npot(0)

    l = lo(i,0)
    call atom_plot_one(1, nrplot, r, rpsi_ps(:,l), RSTEP, l, 1,          &
           'wvf ','    ', .TRUE., itplot, mxdnr)

  enddo

  do i = 1,norbas

    l = lo_b(i)
    call atom_plot_one(1, nrplot, r, rpsi_b(:,i), RSTEP, l, 1,           &
           'bas ','    ', .TRUE., itplot, mxdnr)

  enddo

! plot kinetic energy convergence

  allocate(q(0:nqmax))
  allocate(ekdif(0:nqmax))

  do i = 0,nqbas
    q(i) = i*delqbas
  enddo

  do i = 0,nqbas
    ekdif(i) = ektot(nqbas) - ektot(i)
  enddo

  do i = nqbas,0,-1
    nqp = i
    if(ekdif(i) > SMALL) exit
  enddo

  call atom_plot_one(1, nqp, q, ekdif, QSTEP, 2, 1,                      &
           'ekin','    ', .TRUE., itplot, nqmax)

  deallocate(ekdif)

! plot KB operators

  if(irel == 'rel') then

    do i = 1,npot(1)

      l = lo(i, 1)
      do j = nrplot,1,-1
        nrp_v = j
        if(abs(vkbproj(j,l, 1)) > SMALL) exit
      enddo

      call atom_plot_one(1, nrp_v, r, vkbproj(:,l, 1), RSTEP,            &
            l, inorm(l, 1), 'p   ','r   ', .TRUE.,                       &
            itplot, mxdnr)

    enddo

    do i = 1,npot(-1)

      l = lo(i,-1)
      do j = nrplot,1,-1
        nrp_v = j
        if(abs(vkbproj(j,l,-1)) > SMALL) exit
      enddo

      call atom_plot_one(1, nrp_v, r, vkbproj(:,l,-1), RSTEP,            &
            l, inorm(l,-1), 'p   ','r   ', .TRUE.,                       &
            itplot, mxdnr)

    enddo

  else

    do i = 1,npot(0)

      l = lo(i, 0)
      do j = nrplot,1,-1
        nrp_v = j
        if(abs(vkbproj(j,l, 0)) > SMALL) exit
      enddo

      call atom_plot_one(1, nrp_v, r, vkbproj(:,l, 0), RSTEP,            &
            l, inorm(l, 0), 'p   ','r   ', .TRUE.,                       &
            itplot, mxdnr)

    enddo

  endif


! plot Transforms of KB operators

  do j = 0,nqnl
    nqp = j
    q(j) = j*delql
    if(q(j) > QMAX) exit
  enddo

  if(irel == 'rel') then

    do i = 1,npot(1)

      l = lo(i, 1)

      call atom_plot_one(1, nqp, q, vkbprft(:,l, 1), QSTEP,              &
            l, inorm(l, 1), 'p   ','q   ', .TRUE.,                       &
            itplot, nqmax)

    enddo

    do i = 1,npot(-1)

      l = lo(i,-1)

      call atom_plot_one(1, nqp, q, vkbprft(:,l,-1), QSTEP,              &
            l, inorm(l,-1), 'p   ','q   ', .TRUE.,                       &
            itplot, nqmax)

    enddo

  else

    do i = 1,npot(0)

      l = lo(i, 0)

      call atom_plot_one(1, nqp, q, vkbprft(:,l, 0), QSTEP,              &
            l, inorm(l, 0), 'p   ','q   ', .TRUE.,                       &
            itplot, nqmax)

    enddo

  endif


  close(unit=itplot)

  deallocate(q)

  return

end subroutine atom_kb_plot
