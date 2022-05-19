!>  Calculates the screening potential and total charge
!>
!>  \author       N. Troullier, J.L.Martins
!>  \version      6.0.8
!>  \date         90s, May 2012, 19 May 2022.
!>  \copyright    GNU Public License v2

subroutine atom_kb_screen(nr, r, drdi, cdv, cdc, icorr, ifcore, totvel, vscreen,   &
       iowrite, mxdnr)

! adapted from the old program jlm 22/5/2012
! kb_conv program does not have spin-polarization
! mxdnr, iowrite, 20 September 2021. JLM
! deallocate, 19 May 2022. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension of radial grid points.

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  nr                               !<  number of points in the radial grid
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points   r(i) = a*(exp(b*i)-1), i=1,...,nr
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  real(REAL64), intent(in)          ::  cdc(mxdnr)                       !<  4*pi*r**2 charge density of core
  real(REAL64), intent(in)          ::  cdv(mxdnr)                       !<  4*pi*r**2 charge density of valence.  kb_conv is not spin polarized...
  character(len=2), intent(in)      ::  icorr                            !<  correlation used in the calculation
  integer, intent(in)               ::  ifcore                           !<  0 no partial core correction, 1 partial xc, 2 partial hartree

! output

  real(REAL64), intent(out)         ::  totvel                           !<  total number of valence electrons
  real(REAL64), intent(out)         ::  vscreen(mxdnr)                   !<  screening potential

! allocatable arrays

  real(REAL64), allocatable         ::  y(:), yp(:), ypp(:)
  real(REAL64), allocatable         ::  cdtmp(:,:), vscrtmp(:,:)
  real(REAL64), allocatable         ::  w(:,:)

! local variables

  integer                           ::  ierr
  real(REAL64), dimension(10)       ::  etot(10)

! parameters

  real(REAL64), parameter               ::  ZERO = 0.0_REAL64

! counters

  integer                        ::  i


  allocate(y(mxdnr),yp(mxdnr),ypp(mxdnr))
  allocate(cdtmp(mxdnr,-1:1),vscrtmp(mxdnr,-1:1))
  allocate(w(nr,3))

  y(1) = ZERO
  do i = 2,nr
    y(i) = cdv(i)
  enddo

  call splift(r, y, yp, ypp, nr, w, ierr, 0, ZERO, ZERO, ZERO, ZERO)

  if (ierr /= 1) then
    write(iowrite,'("   error in kb_screen, splift ierr = ",i2)') ierr

    stop

  endif

  call spliq(r, y, yp, ypp, nr, ZERO, r(nr), 1, totvel, ierr)

  if (ierr /= 1) then
    write(iowrite,'(1x,"   error in kb_screen, spliq ierr = ",i2)') ierr

    stop

  endif

!  find el-el potential

  do i = 1,nr
    cdtmp(i,-1) = cdv(i) / 2
    cdtmp(i, 0) = cdv(i)
    cdtmp(i, 1) = cdv(i) / 2
  enddo

  call atom_atm_velect(0, 0, icorr, ' ', ifcore,                         &
      nr, r, drdi, totvel, cdtmp, cdc, vscrtmp, etot,                    &
      iowrite, mxdnr)

  do i = 1,nr
    vscreen(i) = vscrtmp(i,0)
  enddo

  deallocate(y,yp,ypp)
  allocate(cdtmp,vscrtmp)
  allocate(w)

  return

end subroutine atom_kb_screen
