!>  Find the 1st and 2nd eigenvalues for all angular
!>  momentum using just the local potential.  Print
!>  out results, this is used to find if any ghost
!>  states exist for the potential.
!>  See Gonze, Kackell, and Scheffler, Phy. Rev. B. 41, 12264 (1990)
!>
!>  \author       J.L.Martins
!>  \version      6.0.8
!>  \date         early 90s, May 2012, 22 May 2022.
!>  \copyright    GNU Public License v2

subroutine atom_kb_ghost_test(npot, lo, llocal, irel, nr, r, drdi, d2rodr,    &
         vscreen, vlocal, ev, inorm, irayps,                                  &
         iowrite, mxdl, mxdnr)

! adapted from the old program jlm 22/5/2012
! converted to f90 10/6/2012
! mxdnr, mxdl, 17 August 2021. JLM
! Tolerance for scattering states. 22 May 2022. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum components
  integer, intent(in)               ::  mxdnr                            !<  dimension of radial grid points.

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  npot(-1:1)                       !<  number of orbitals to be calculated
  integer, intent(in)               ::  lo(mxdl+1,-1:1)                  !<  angular momentum of orbital
  integer, intent(in)               ::  llocal                           !< angular momentum for local potential (negative: maximum of l-dependent)

  character(len=3), intent(in)      ::  irel                             !<  flag for relativistic (r) or spin-polarized (s) original calculation

  integer, intent(in)               ::  nr                               !<  number of points in the radial grid
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points   r(i) = a*(exp(b*i)-1), i=1,...,nr
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  (d r(i) / d i)
  real(REAL64), intent(in)          ::  d2rodr(mxdnr)                    !<  (d^2 r(i) / d i^2) / (d r / d i)
  real(REAL64), intent(in)          ::  vscreen(mxdnr)                   !<  screening potential
  real(REAL64), intent(in)          ::  vlocal(mxdnr)                    !<  r*local pseudopotential
  real(REAL64), intent(in)          ::  ev(0:mxdl,-1:1)                  !<  eigenvalues (l,2j-2l)
  integer, intent(in)               ::  inorm(0:mxdl,-1:1)               !<  sign of denominator of KB operator

! input and output

  character(len=10), intent(inout)     :: irayps(4)                      !  type of pseudopotential

! local variables

  integer                        ::  l, jmax, iflag
  real(REAL64)                   ::  evt
  real(REAL64)                   ::  delta

! allocatable arrays

  real(REAL64), allocatable      ::  ev0(:), ev1(:)
  real(REAL64), allocatable      ::  rpsi(:), drpsidr(:)
  real(REAL64), allocatable      ::  v(:)

! parameters

  real(REAL64), parameter          :: ZERO = 0.0_REAL64
  real(REAL64), parameter          :: TOL = 1.0E-06_REAL64

! counters

  integer                        :: i, j, k


  allocate(ev0(0:mxdl),ev1(0:mxdl))
  allocate(rpsi(mxdnr),drpsidr(mxdnr))
  allocate(v(mxdnr))

  j = 0
  do i = 1,npot(j)

    l = lo(i,j)
    if (llocal /= l) then
      evt = ev(l,j)

      v(1) = zero
      do k = 2,nr
        v(k) = (vlocal(k)  + (l*(l+1)) /r(k)) / r(k)  + vscreen(k)
      enddo

      ev0(l) = evt
      call atom_atm_difnrl(nr, r, drdi, d2rodr, v, rpsi, drpsidr,        &
       l+1, l, ev0(l), iflag, TOL,                                       &
       iowrite, mxdnr)

      ev1(l) = evt
      call atom_atm_difnrl(nr, r, drdi, d2rodr, v, rpsi, drpsidr,        &
       l+2, l, ev1(l), iflag, TOL,                                       &
       iowrite, mxdnr)

    endif

  enddo

  jmax = 0
  if(irel == 'rel') jmax = 1

  write(iowrite,*)
  write(iowrite,*)
  write(iowrite,'("  Ghost State Test (X. Gonze et al)")')
  write(iowrite,*)
  write(iowrite,'("  l     j",4x,"0-node-eigen    1-node-eigen ",        &
         &   "      True-eigen   inorm")')
  write(iowrite,*)

  do j = -jmax,jmax

    do i = 1,npot(j)
      l = lo(i,j)

      if (llocal /= l) then
        if(j == 0) then
          write(iowrite,'(1x,i2,9x,3(f12.6,4x),i2)')                     &
             l, ev0(l), ev1(l), ev(l,j), inorm(l,j)
        else
          write(iowrite,'(1x,i2,2x,i3,"/2  ",3(f12.6,4x),i2)')           &
             l, 2*l+j, ev0(l), ev1(l), ev(l,j), inorm(l,j)
        endif

        delta = ZERO
        if(ev(l,j) > ZERO) delta = TOL

        if( (inorm(l,j) < 0 .and. ev(l,j) > ev0(l) + delta) .or.         &
             (inorm(l,j) > 0 .and. ev(l,j) > ev1(l) + delta) ) then
          write(iowrite,'(//,"  WARNING:     GHOST STATE   GHOST!!!",//)')
          if(iowrite /= 6) then
            write(6,'(//,"  WARNING:     GHOST STATE   GHOST!!!",//)')
          endif
          irayps(1) = ' G H O S T'
          irayps(2) = '  S T A T '
          irayps(3) = 'E  P R E S'
          irayps(4) = ' E N T !!!'
        endif

      endif

    enddo

  enddo

  deallocate(ev0,ev1)
  deallocate(rpsi,drpsidr)
  deallocate(v)

  return

end subroutine atom_kb_ghost_test
