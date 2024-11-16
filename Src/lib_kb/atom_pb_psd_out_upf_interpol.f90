!>  Interpolates from the atom code grid to the UPF grid,
!>  the local and non-local pseudopotentials, the wavefunctions
!>  and charge densities (core and valence)
!>
!>  \author       J.L.Martins
!>  \version      6.0.9
!>  \date         15 November 2024.
!>  \copyright    GNU Public License v2

subroutine atom_kb_psd_out_upf_interp(nr, r, nrupf, rupf, lso,           &
      vlocal, lmax_pot, vkbproj, lmax_psi, rpsi, cdv, cdc,               &
      vlocupf, vnlupf, chi, cdvupf, cdcupf,                              &
      mxdl, mxdnr)

! Extracted from atom_kb_psd_out_upf. 15 November 2024. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! dimensions:

  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum components
  integer, intent(in)               ::  mxdnr                            !<  dimension of radial grid points.

! input:

  integer, intent(in)               ::  nr                               !<  number of points in the radial grid
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points   r(i) = a*(exp(b*(i-1))-1), i=1,...,nr

  integer, intent(in)               ::  nrupf                            !<  number of points in the UPF grid
  real(REAL64), intent(in)          ::  rupf(nrupf)                      !<  radial grid points   rupf(i) = a*(exp(b*(i-1))), i=1,...,nrupf

  logical, intent(in)               ::  lso                              !<  indicates that spin-orbit is used

  real(REAL64), intent(in)          ::  vlocal(mxdnr)                    !<  r*local pseudopotential
  real(REAL64), intent(in)          ::  cdc(mxdnr)                       !<  4*pi*r**2 charge density of valence.  kb_conv is not spin polarized...
  real(REAL64), intent(in)          ::  cdv(mxdnr)                       !<  4*pi*r**2 charge density of core

  integer                           ::  lmax_pot                         !<  maximum angular momentum in potential
  real(REAL64), intent(in)          ::  vkbproj(mxdnr,0:mxdl,-1:1)       !<  kb-projector(r(i),l,2j-2l)

  integer, intent(in)               ::  lmax_psi                         !<  maximum angular momentum in wave-functions
  real(REAL64), intent(in)          ::  rpsi(mxdnr,0:mxdl,-1:1)          !<  wavefunctions (r(i),l,2j-2l).

! output

  real(REAL64), intent(out)         ::  vlocupf(nrupf)                   !<  local pseudopotential (not r*vlocal) for UPF
  real(REAL64), intent(out)         ::  vnlupf(nrupf,0:mxdl,-1:1)        !<  kb-projectors (multiplied by r) for UPF
  real(REAL64), intent(out)         ::  chi(nrupf,0:mxdl,-1:1)           !<  wavefunctions (r(i),l,2j-2l) for UPF (r*psi)
  real(REAL64), intent(out)         ::  cdcupf(nrupf)                    !<  charge density of core.
  real(REAL64), intent(out)         ::  cdvupf(nrupf)                    !<  4*pi*r**2 charge density of valence.

! spline interpolation

  real(REAL64), allocatable        ::  yp(:),ypp(:),wlu(:,:)
  real(REAL64), allocatable        ::  ypi(:),yppi(:)

! other variables

  integer                          ::  ierr
  integer                          ::  jmin, jmax

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  PI4 = 16*atan(1.0_REAL64)

! counters

  integer       ::   i, j, l


  allocate(yp(nr),ypp(nr),wlu(nr,3))
  allocate(ypi(nrupf),yppi(nrupf))

! interpolates local potential

  call splift (r, vlocal, yp, ypp, nr, wlu, ierr, 0, ZERO,ZERO,ZERO,ZERO)

  if(ierr /= 1) then
    write(6,*) '   STOPPED in atom_kb_psd_out_upf_interp'
    write(6,*) '   error in splift vlocal'
    stop
  endif

  call splint (r, vlocal, ypp, nr, rupf, vlocupf, ypi, yppi, nrupf, ierr)

  if(ierr /= 1) then
    write(6,*) '   STOPPED in atom_kb_psd_out_upf_interp'
    write(6,*) '   error in splint vlocal'
    stop
  endif

  do i = 1,nrupf
    vlocupf(i) = vlocupf(i)/rupf(i)
  enddo


! interpolates KB projectors

  do l = 0,lmax_pot

    if(lso) then
      jmin = -1
      jmax = 1
      if(l == 0) jmin = 0
    else
      jmin = 0
      jmax = 0
    endif

    do j = jmin,jmax

      call splift (r, vkbproj(:,l,j), yp, ypp, nr, wlu, ierr, 1, ZERO,ZERO,ZERO,ZERO)

      if(ierr /= 1) then
        write(6,*) '   STOPPED in atom_kb_psd_out_upf_interp'
        write(6,*) '   error in splift vkbproj  ', l
        stop
      endif

      call splint (r, vkbproj(:,l,j), ypp, nr, rupf, vnlupf(1,l,j), ypi, yppi, nrupf, ierr)

      if(ierr /= 1) then
        write(6,*) '   STOPPED in atom_kb_psd_out_upf_interp'
        write(6,*) '   error in splint vkbproj  ', l
        stop
      endif

      do i = 1,nrupf
        vnlupf(i,l,j) = vnlupf(i,l,j)*rupf(i)
      enddo

    enddo

  enddo

! interpolates wave-functions

  do l = 0,lmax_psi

    if(lso) then
      jmin = -1
      jmax = 1
      if(l == 0) jmin = 0
    else
      jmin = 0
      jmax = 0
    endif

    do j = jmin,jmax

      call splift (r, rpsi(:,l,j), yp, ypp, nr, wlu, ierr, 1, ZERO,ZERO,ZERO,ZERO)

      if(ierr /= 1) then
        write(6,*) '   STOPPED in atom_kb_psd_out_upf_interp'
        write(6,*) '   error in splift rpi', l
        stop
      endif

      call splint (r, rpsi(:,l,j), ypp, nr, rupf, chi(:,l,j), ypi, yppi, nrupf, ierr)

      if(ierr /= 1) then
        write(6,*) '   STOPPED in atom_kb_psd_out_upf_interp'
        write(6,*) '   error in splint rpsi  ', l
        stop
      endif

    enddo

  enddo

! interpolates charge densities

  call splift (r, cdv, yp, ypp, nr, wlu, ierr, 1, ZERO,ZERO,ZERO,ZERO)

  if(ierr /= 1) then
    write(6,*) '   STOPPED in atom_kb_psd_out_upf_interp'
    write(6,*) '   error in splift valence charge'
    stop
  endif

  call splint (r, cdv, ypp, nr, rupf, cdvupf, ypi, yppi, nrupf, ierr)

  if(ierr /= 1) then
    write(6,*) '   STOPPED in atom_kb_psd_out_upf_interp'
    write(6,*) '   error in splint valence charge'
    stop
  endif


  call splift (r, cdc, yp, ypp, nr, wlu, ierr, 1, ZERO,ZERO,ZERO,ZERO)

  if(ierr /= 1) then
    write(6,*) '   STOPPED in atom_kb_psd_out_upf_interp'
    write(6,*) '   error in splift core charge'
    stop
  endif

  call splint (r, cdc, ypp, nr, rupf, cdcupf, ypi, yppi, nrupf, ierr)

  if(ierr /= 1) then
    write(6,*) '   STOPPED in atom_kb_psd_out_upf_interp'
    write(6,*) '   error in splint core charge'
    stop
  endif

  do i = 1,nrupf
    cdcupf(i) = cdcupf(i) / (PI4*rupf(i)*rupf(i))
  enddo


  return

end subroutine atom_kb_psd_out_upf_interp

