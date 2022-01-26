!>  Printout data about the orbital.
!>
!>  \author       Sverre Froyen,  Norm Troullier, Jose Luis Martins
!>  \version      6.0.6
!>  \date         1980s, 22 June 2021, 31 October 2021.
!>  \copyright    GNU Public License v2

  subroutine atom_atm_orban(ispp, ar, br, nr, r, drdi,                   &
      no, lo, zo, iso ,vionic, vhxc, ev, ek, ep,                         &
      iowrite, mxdnr, mxdl)

! converted to f90, June 2018
! cleanup and new interface, July 2019. JLM
! jlm  version 6.00
! vionic, so->iso. single orbital. 12 September 2021. JLM
! call to atom-atom_atm_zero_extr. 31 october 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points
  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum l

  character(len=1), intent(in)      ::  ispp                             !<  spin polarization  ' ', 's', 'r'

  real(REAL64), intent(in)          ::  ar(mxdnr)                        !<  radial wave-function u = rR(r)  (integral ar**2 = 1) or major component radial wave-function u = rR(r)
  real(REAL64), intent(in)          ::  br(mxdnr)                        !<  d ar / d r    or  minor component (integral ar**2 + br**2 = 1)

  integer, intent(in)               ::  nr                               !<  dimension of the number of radial points

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  integer, intent(in)               ::  no                               !<  principal quantum number n
  integer, intent(in)               ::  lo                               !<  angular quantum number l
  integer, intent(in)               ::  iso                              !<  2*spin or 2*(j-l)
  real(REAL64), intent(in)          ::  zo                               !<  orbital occupation

  real(REAL64), intent(in)          ::  vionic(mxdnr)                    !<  r*ionic potential for orbital (Ry)

  real(REAL64), intent(in)          ::  vhxc(mxdnr)                      !<  screening potential (Ry)

  real(REAL64), intent(in)          ::  ev                               !<  orbital energy

! output

  real(REAL64), intent(out)         ::  ek                               !<  orbital kinetic energy
  real(REAL64), intent(out)         ::  ep                               !<  orbital potential energy

! local alocatable arrays

  real(REAL64), allocatable         ::  rzero(:)
  real(REAL64), allocatable         ::  rextr(:)
  real(REAL64), allocatable         ::  aextr(:)
  real(REAL64), allocatable         ::  bextr(:)

! local

  character(len=10)      ::  name

  integer                ::  lp
  integer                ::  nn
  integer                ::  nzero                                      !  number of zeroes
  integer                ::  nextr                                      !  number of extrema

  real(REAL64)           ::  ar2, br2
  real(REAL64)           ::  deni, sa2

  integer                ::  i90, i99                                   !  radial index for 90% or 99% of charge
  integer                ::  ll, llp

! constants

  real(REAL64), parameter   ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64

! counters

  integer     ::  i


  if(lo > mxdl) then
    write(iowrite,*) '   error in atom_atm_orban, angular momentum ',    &
        lo,' is greater than mxdl', mxdl

    stop

  endif

! allocate arrays

  nn = no-lo + 2

  allocate(rzero(nn))
  allocate(rextr(nn))
  allocate(aextr(nn))
  allocate(bextr(nn))

! finds zeroes and extrema

  call atom_atm_zero_extr(ispp, ar, br, nr, r,                           &
      no, lo, iso ,vionic, vhxc, ev,                                     &
      nzero, nextr, rzero, rextr, aextr, bextr,                          &
      mxdnr)


! find orbital kinetic and potential energy
! the potential part includes only the interaction with
! the nuclear part

  ek = br(1)*br(1)*drdi(1)
  ep = ZERO
  sa2 = ZERO
  lp = lo+1
  llp = lo*lp

  ll = 2
  if (2*(nr/2) == nr) ll=4
  i90 = nr
  i99 = nr

  do i = nr,2,-1
    ar2 = ar(i)*ar(i)
    br2 = br(i)*br(i)
    deni = ar2
    if (ispp == 'r') deni = deni + br2
    ek = ek + ll * (br2 + ar2*llp/r(i)**2)*drdi(i)
    ep = ep + ll * deni*vionic(i)*drdi(i) / r(i)
    ll = 6 - ll

    if (sa2 <= 0.1) then
      sa2 = sa2 + deni*drdi(i)
      if (sa2 <= 0.01) i99 = i
      i90 = i
    endif

  enddo

  ek = ek / 3
  ep = ep / 3
  if (ispp == 'r') ek = ZERO

! printout

  write(iowrite,*)
  write(iowrite,'(" n =",i2,"  l =",i2,"  s =",f4.1)')                   &
            no, lo, iso*0.5

  name = 'a extr    '
  write(iowrite,'(8x,a10,2x,8f8.3)') name,(aextr(i),i=1,nextr)
  if (ispp == 'r') then
    name = 'b extr    '
    write(iowrite,'(8x,a10,2x,8f8.3)') name,(bextr(i),i=1,nextr)
  endif
  name = 'r extr    '
  write(iowrite,'(8x,a10,2x,8f8.3)') name,(rextr(i),i=1,nextr)
  name = 'r zero    '
  write(iowrite,'(8x,a10,2x,8f8.3)') name,(rzero(i),i=1,nzero)
  name = 'r 90/99 % '
  write(iowrite,'(8x,a10,2x,8f8.3)') name,r(i90),r(i99)
  if (ev == ZERO) then
    if (zo /= ZERO) then
      write(iowrite,'(8x,"WARNING: This orbital is not bound",           &
        & " and contains ",f6.4," electrons!!")') zo
    else
      write(iowrite,'(8x,"WARNING:  This orbital is not bound!")')
    endif
  endif

! njtj  ***  plotting routines  ***
! jlm   uncomment dimension declaration for wk1,etc
!
!   Save plotting information to current plot.dat file
! (unit = 3),  User must specify what orbital
!  is to be saved(or all).
!
!     ist=1
!     if (ar(nr-80) < 0.0) ist=-1
!     call potrw(ar,r,nr-85,lo,1,ist)
!     call wtrans(ar,r,nr,drdi,lo,ist,wk1)
!
!      do i=2,nr
!        v(i)=vionic(i)/r(i)+vhxc(i)
!      enddo
!      zion=4
!      dimension wk1(1000),wk2(1000),wk3(1000),v(1000)
!      call potran(lo+1,v,r,nr,zion,wk1,wk2,wk3)
!      call potrv(v,r,nr,lo)
!
! njtj  ***  user should adjust for their needs  ***

  deallocate(rzero)
  deallocate(rextr)
  deallocate(aextr)
  deallocate(bextr)

  return
  end subroutine atom_atm_orban
