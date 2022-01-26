!>     writes result of self-consistent calculation to file
!>      'datafile.dat' for latter use
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.0.3
!>  \date         May 1, 1991, 22 June 2021, 12 September 2021.
!>  \copyright    GNU Public License v2

subroutine atom_atm_datout(ioae, fileae, itype, icorr, ispp,             &
         nr, a, b, r, drdi,                                              &
         nameat, norb, ncore, no, lo, iso, zo, znuc, zel,                &
         cdc, vionic, vhxc, ev,                                          &
         mxdnr, mxdorb, mxdl)

!  ***********************************************************
!  *                                                         *
!  *  Users may want to remove or modify this routine        *
!  *  depending on their needs.                              *
!  *                                                         *
!  *  Version dated May 1, 1991                              *
!  *  njtj                                                   *
!  *                                                         *
!  ***********************************************************

! converted to f90, April 2018
! cleanup and new interface, July 2019. JLM
! added date and version, 27 August 2021. JLM
! vionic, so->iso. 12 September 2021. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  ioae                             !<  default tape for writing
  character(len=*), intent(in)      ::  fileae                           !<  name of default tape for all-electron results

  integer, intent(in)               ::  mxdorb                           !<  dimension of the number of orbitals
  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points
  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum l

  integer, intent(in)               ::  nr                               !<  dimension of the number of radial points

  integer, intent(in)               ::  norb                             !<  dimension of the number of orbitals

  integer, intent(in)               ::  itype                            !<  type of calculation (-1 signals end of calculation)
  character(len=2), intent(in)      ::  icorr                            !<  correlation type
  character(len=1), intent(in)      ::  ispp                             !<  spin polarization  ' ', 's', 'r'

  real(REAL64), intent(in)          ::  a                                !<  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64), intent(in)          ::  b                                !<  r(i) = a*(exp(b*(i-1))-1)

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbot of the atom

  integer, intent(in)               ::  ncore                            !<  number of orbitals treated as core

  integer, intent(in)               ::  no(mxdorb)                       !<  principal quantum number n
  integer, intent(in)               ::  lo(mxdorb)                       !<  angular quantum number l
  integer, intent(in)               ::  iso(mxdorb)                      !<  2*spin or 2*(j-l)
  real(REAL64), intent(in)          ::  zo(mxdorb)                       !<  orbital occupation

  real(REAL64), intent(in)          ::  znuc                             !<  nuclear charge
  real(REAL64), intent(in)          ::  zel                              !<  electron charge

  real(REAL64), intent(in)          ::  cdc(mxdnr)                       !<  core charge density (total)

  real(REAL64), intent(in)          ::  vionic(mxdnr,0:mxdl,-1:1)        !<  r*ionic potential in Rydberg; -1:  j=l-1/2 or s= -1/2.  0:  average or spinless.  1:  j=l+1/2 or s=1/2

  real(REAL64), intent(in)          ::  vhxc(mxdnr,-1:1)                 !<  screening potential (Ry)

  real(REAL64), intent(in)          ::  ev(mxdorb)                       !<  orbital energy

! local variables, maybe should be input...

  character(len=9)                  ::  bdate                            !  date the subroutine was called
  character(len=5)                  ::  vers                             !  version of the code

! other variables

  integer                           ::  lmax                             !  maximum angular momentum l
  integer                           ::  id

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, PFIVE = 0.5_REAL64

! counter

  integer       ::  i, l


  if(ioae == 0) RETURN

  if(ispp == ' ') then
    id = 0
  else
    id = 1
  endif

! gets the date and version.   New in 6.014

  call zedate(bdate)
  call atom_version(vers)

  lmax = 0
  do i = 1,norb
    if(lmax < lo(i)) lmax = lo(i)
  enddo

! Open and write out data to file datafile.dat.

  open (unit=ioae, file=trim(fileae), status='unknown', form='unformatted')

  write(ioae) itype, icorr, ispp, nr, a, b, bdate, vers
  write(ioae) (r(i),i=1,nr)
  write(ioae) (drdi(i),i=1,nr)
! convention is lmax+1
  write(ioae) lmax+1, nameat, norb, ncore

  write(ioae) (no(i),i=1,norb)
  write(ioae) (lo(i),i=1,norb)
  write(ioae) (iso(i)*PFIVE,i=1,norb)
  write(ioae) (zo(i),i=1,norb)
  write(ioae) znuc, zel
  write(ioae) (cdc(i),i=1,nr)

  do l = 0,lmax
    write(ioae) (vionic(i,l,id),i=1,nr)
  enddo
  do l = 0,lmax
    write(ioae) (vionic(i,l,-1),i=1,nr)
  enddo

  write(ioae) (vhxc(i,id),i=1,nr)
  write(ioae) (vhxc(i,-1),i=1,nr)
  write(ioae) (ev(i),i=1,norb)

  close (unit=ioae)

  return

end subroutine atom_atm_datout
