!>  reads calculation parameters, by default from atom.dat (set by ioread)
!>
!>  \author       Sverre Froyen, Norm Troullier, Carlos Balbas, Jose Luis Martins
!>  \version      6.013
!>  \date         4 August 2021
!>  \copyright    GNU Public License v2

subroutine atom_atm_read_input(nameat, ctype, ititle, kerker, icorr, ispp,  &
    znuc, rmax, aa, bb, ncore, nval, ni, li, zd, zu, evd,                   &
    ioread, mxdorb)


! njtj ***  modifications  ***
!   The input and output variables passed have been changed.
!   There are five new pseudopotential generation options
!   The input variables znuc,zsh,rsh,rmax,aa,bb are
!   compared to a small positive value - eliminates
!   floating point comparisions errors(zero is
!   not always zero).
! njtj ***  modifications  ***

! lcb
!   modified for GGA by Carlos Balbas,   January 97.
! lcb

! cleanup and new interface, July 2019. JLM
! minor cleanup, new name, 22 June 2021. JLM
! jlm  version 6.011
! jlm adapted from atm_input.f90.  zsh, rsh deprecated


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

!  data for orbitals, 1s,2s,2p,3s,3p,3d,4s,4p,4d,5s,5p,4f,5d,6s,6p,5f,6d,7s,7p

  integer, parameter    ::  ncmax = 19                                   !  maximum number of core orbitals
  integer, parameter    ::  nc(ncmax) = (/1,2,2,3,3,3,4,4,4,5,5,4,5,6,6,5,6,7,7/)
  integer, parameter    ::  lc(ncmax) = (/0,0,1,0,1,2,0,1,2,0,1,3,2,0,1,3,2,0,1/)


! input

  integer, intent(in)               ::  ioread                           !<  default tape for reading

  integer, intent(in)               ::  mxdorb                           !<  dimension of the number of orbitals

! output

  character(len=2), intent(out)     ::  nameat                           !<  chemical symbol of the atom

  character(len=2), intent(out)     ::  ctype                            !<  type of calculation flag
  character(len=10), intent(out)    ::  ititle(5)                        !<  title of calculation

  character(len=3), intent(out)     ::  kerker                           !<  type of pseudopotential flag
  character(len=2), intent(out)     ::  icorr                            !<  correlation type
  character(len=1), intent(out)     ::  ispp                             !<  spin polarization  ' ', 's', 'r'

  real(REAL64), intent(out)         ::  znuc                             !<  nuclear charge
  real(REAL64), intent(out)         ::  rmax                             !<  maximum mesh radius
  real(REAL64), intent(out)         ::  aa                               !<  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64), intent(out)         ::  bb                               !<  a = exp(-aa)/znuc, b = 1/bb

  integer, intent(out)              ::  ncore                            !<  number of orbitals treated as core
  integer, intent(out)              ::  nval                             !<  number of valence orbitals

  integer, intent(out)              ::  ni(mxdorb)                       !<  valence principal quantum number
  integer, intent(out)              ::  li(mxdorb)                       !<  valence angular quantum number
  real(REAL64), intent(out)         ::  zd(mxdorb), zu(mxdorb)           !<  occupation of valence down and up orbitals
  real(REAL64), intent(out)         ::  evd(mxdorb)                      !<  default eigenvalue

! local

  real(REAL64)      ::  zsh                                              !  shell charge (deprecated)
  real(REAL64)      ::  rsh                                              !  shell radius (deprecated)
  integer           ::  ios

! counters

  integer     ::  i


!  read the type of calculation and title card

  ctype = '  '
  do i = 1,5
    ititle(i) = '          '
  enddo

  read(ioread,'(3x,a2,5a10)',iostat = ios) ctype,ititle

!  if ctype = ' ' , no more data, program ends

  if (ctype == '  ' .or. ios /= 0) then

    ctype = '  '

    return

  endif

  if (ctype /= 'ae') then
    read(ioread,*) kerker
  endif


! njtj  ***  major modification end  ***

!   read element name and correlation type
!   ispp = ' ' - nonspin calculation
!   ispp = s  - spin polarized calculation
!   ispp = r  - relativistic calculation

  read(ioread,'(3x,a2,3x,a2,a1)') nameat, icorr, ispp

! njtj   ***  major modification start  ***
!   Floating point comparison error modification.
!   Read the atomic number (nuclear charge),
!   shell charge and radius (added to the nuclear potential),
!   and radial grid parameters.

  read(ioread,*) znuc,zsh,rsh,rmax,aa,bb

! read the number of core and valence orbitals

  read(ioread,*) ncore,nval

! reads valence information

  if (nval /= 0) then

    do i = 1,nval

      read(ioread,'(2i5,4f10.3)') ni(i),li(i),zd(i),zu(i),evd(i)

    enddo

  endif

  return

end subroutine atom_atm_read_input
