!>  reads calculation parameters, by default from atom.dat (set by ioread)
!>
!>  \author       Sverre Froyen, Jose Luis Martins
!>  \version      6.0.2
!>  \date         4 August 2021
!>  \copyright    GNU Public License v2

subroutine atom_kb_test_read_input(nameat, ititle, icorr, ispp,          &
    nval, ni, li, zd, zu, evd,                                           &
    ioread, mxdorb)



! cleanup and new interface, July 2019. JLM
! minor cleanup, new name, 22 June 2021. JLM
! jlm adapted from atm_input.f90.  zsh, rsh deprecated


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  ioread                           !<  default tape for reading

  integer, intent(in)               ::  mxdorb                           !<  dimension of the number of orbitals

! output

  character(len=2), intent(out)     ::  nameat                           !<  chemical symbol of the atom

  character(len=10), intent(out)    ::  ititle(5)                        !<  title of calculation

  character(len=2), intent(out)     ::  icorr                            !<  correlation type
  character(len=1), intent(out)     ::  ispp                             !<  spin polarization  ' ', 's', 'r'

  integer, intent(out)              ::  nval                             !<  number of valence orbitals

  integer, intent(out)              ::  ni(mxdorb)                       !<  valence principal quantum number
  integer, intent(out)              ::  li(mxdorb)                       !<  valence angular quantum number
  real(REAL64), intent(out)         ::  zd(mxdorb), zu(mxdorb)           !<  occupation of valence down and up orbitals
  real(REAL64), intent(out)         ::  evd(mxdorb)                      !<  default eigenvalue

! local

  integer              ::  ios
  integer              ::  ncore                                         !  number of orbitals treated as core
  character(len=2)     ::  ctype                                         !  type of calculation flag

! parameters

  real(REAL64), parameter    ::  UM = 1.0_REAL64, ZERO = 0.0_REAL64

! counters

  integer     ::  i


  evd = ZERO

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
    read(ioread,*)
  endif


! njtj  ***  major modification end  ***

!   read element name and correlation type
!   ispp = ' ' - nonspin calculation
!   ispp = s  - spin polarized calculation
!   ispp = r  - relativistic calculation

  read(ioread,'(3x,a2,3x,a2,a1)') nameat, icorr, ispp

  read(ioread,*)

! read the number of core and valence orbitals

  read(ioread,*) ncore,nval

! reads valence information

  if (nval > 0) then

    do i = 1,nval

      read(ioread,'(2i5,4f10.3)') ni(i),li(i),zd(i),zu(i),evd(i)

    enddo

  endif

  return

end subroutine atom_kb_test_read_input
