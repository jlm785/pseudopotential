!>  Printout data about the orbital.  KB version
!>
!>  \author       Sverre Froyen,  Norm Troullier, Jose Luis Martins
!>  \version      6.0.6
!>  \date         1980s, 22 June 2021. 31 October 2021.
!>  \copyright    GNU Public License v2

  subroutine atom_kb_test_orban(iorb, ar, br, nr, r, drdi,              &
      norb, no, lo, zo, iso ,v, ev, ek, ep,                             &
      iowrite, mxdnr, mxdorb, mxdl)

! converted to f90, June 2018
! cleanup and new interface, July 2019. JLM
! mxdnr, mxdnr. 19 September 2021. JLM
! atom_atm_zero_extr. 31 october 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdorb                           !<  dimension of the number of orbitals
  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points
  integer, intent(in)               ::  mxdl                             !<  dimension angular momentum

  integer, intent(in)               ::  iorb                             !<  orbital index

  real(REAL64), intent(in)          ::  ar(mxdnr)                        !<  radial wave-function u = rR(r)  (integral ar**2 = 1) or major component radial wave-function u = rR(r)
  real(REAL64), intent(in)          ::  br(mxdnr)                        !<  d ar / d r    or  minor component (integral ar**2 + br**2 = 1)

  integer, intent(in)               ::  nr                               !<  dimension of the number of radial points

  real(REAL64), intent(in)          ::  r(mxdnr)                       !<  radial grid points
  real(REAL64), intent(in)          ::  drdi(mxdnr)                    !<  d r(i) / d i

  integer, intent(in)               ::  norb                             !<  number of orbitals

  integer, intent(in)               ::  no(mxdorb)                       !<  principal quantum number n
  integer, intent(in)               ::  lo(mxdorb)                       !<  angular quantum number l
  integer, intent(in)               ::  iso(mxdorb)                      !<  2*spin or 2*(j-l)
  real(REAL64), intent(in)          ::  zo(mxdorb)                       !<  orbital occupation

  real(REAL64), intent(in)          ::  v(mxdnr)                         !<  effective potential in Rydberg

  real(REAL64), intent(in)          ::  ev(mxdorb)                       !<  orbital energy

! output

  real(REAL64), intent(out)         ::  ek(mxdorb)                       !<  orbital kinetic energy
  real(REAL64), intent(out)         ::  ep(mxdorb)                       !<  orbital potential energy

! local alocatable arrays

  real(REAL64), allocatable         ::  rzero(:)
  real(REAL64), allocatable         ::  rextr(:)
  real(REAL64), allocatable         ::  aextr(:)
  real(REAL64), allocatable         ::  bextr(:)

  real(REAL64), allocatable         ::  vionic(:)
  real(REAL64), allocatable         ::  vhxc(:)

! local

  character(len=10)      ::  name

  integer                ::  lp
  integer                ::  nn
  integer                ::  nzero                                      !  number of zeroes
  integer                ::  nextr                                      !  number of extrema

  real(REAL64)           ::  ar2, br2
  real(REAL64)           ::  sa2

  integer                ::  i90, i99                                   !  radial index for 90% or 99% of charge
  integer                ::  ll, llp

! counters

  integer     ::  i

! constants

  real(REAL64), parameter   ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64

  if(iorb > norb) then
    write(iowrite,*) '   error in atom_kb_test_orban, orbital index ',   &
        iorb, ' is greater than the number of orbitals', norb

    STOP

  endif

  if(lo(iorb) > mxdl) then
    write(iowrite,*) '   error in atom_kb_test_orban, angular momentum ',  &
        lo(iorb),' is greater than mxdl', mxdl

    STOP

  endif

! allocate arrays

  nn = no(iorb)-lo(iorb) + 2

  allocate(rzero(nn))
  allocate(rextr(nn))
  allocate(aextr(nn))
  allocate(bextr(nn))

! vionic and vhxc are not used if non-relativistic (neither is bextr)

  allocate(vionic(mxdnr))
  allocate(vhxc(mxdnr))

! compute zeroes and extrema

  call atom_atm_zero_extr(' ', ar, br, nr, r,                            &
      no(iorb), lo(iorb), iso(iorb), vionic, vhxc, ev(iorb),             &
      nzero, nextr, rzero, rextr, aextr, bextr,                          &
      mxdnr)

  deallocate(vionic)
  deallocate(vhxc)

! find orbital kinetic and potential energy
! the potential part includes only the interaction with
! the nuclear part

  ek(iorb) = br(1)*br(1)*drdi(1)
  ep(iorb) = ZERO
  sa2 = ZERO
  lp = lo(iorb)+1
  llp = lo(iorb)*lp

  ll = 2
  if (2*(nr/2) == nr) ll = 4
  i90 = nr
  i99 = nr

  do i = nr,2,-1
    ar2 = ar(i)*ar(i)
    br2 = br(i)*br(i)
    ek(iorb) = ek(iorb) + ll * (br2 + ar2*llp/(r(i)*r(i)))*drdi(i)
    ep(iorb) = ep(iorb) + ll * ar2*v(i)*drdi(i)
    ll = 6 - ll

    if (sa2 <= 0.1) then
      sa2 = sa2 + ar2*drdi(i)
      if (sa2 <= 0.01) i99 = i
      i90 = i
    endif

  enddo

  ek(iorb) = ek(iorb) / 3
  ep(iorb) = ep(iorb) / 3

! printout

  write(iowrite,*)
  write(iowrite,'(" n =",i2,"  l =",i2,"  s =",f4.1)')                   &
            no(iorb), lo(iorb), 0.5*iso(iorb)

  name = 'a extr    '
  write(iowrite,'(8x,a10,2x,8f8.3)') name,(aextr(i),i=1,nextr)
  name = 'r extr    '
  write(iowrite,'(8x,a10,2x,8f8.3)') name,(rextr(i),i=1,nextr)
  name = 'r zero    '
  write(iowrite,'(8x,a10,2x,8f8.3)') name,(rzero(i),i=1,nzero)
  name = 'r 90/99 % '
  write(iowrite,'(8x,a10,2x,8f8.3)') name,r(i90),r(i99)
  if (ev(iorb) == ZERO) then
    if (zo(iorb) /= ZERO) then
      write(iowrite,'(8x,"WARNING: This orbital is not bound",           &
        & " and contains ",f6.4," electrons!!")') zo(iorb)
    else
      write(iowrite,'(8x,"WARNING:  This orbital is not bound!")')
    endif
  endif


  deallocate(rzero)
  deallocate(rextr)
  deallocate(aextr)
  deallocate(bextr)

  return
  end subroutine atom_kb_test_orban
