!>  reads by default from atom.dat and obtain the required array dimensions
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.8
!>  \date         22 June 2021. 19 May 2022.
!>  \copyright    GNU Public License v2

subroutine atom_kb_test_input_size(ioread, iopsd, filepsd, nr, norb, lmax)

! mxdl. 18 September 2021. JLM
! psdtitle, 19 May 2022. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  ioread                           !<  default tape for reading

  integer, intent(in)               ::  iopsd                            !<  default tape for KB pseudopotential
  character(len=*), intent(in)      ::  filepsd                          !<  name of default tape for reading KB pseudopotential

! output

  integer, intent(out)              ::  norb                             !<  number of orbitals
  integer, intent(out)              ::  nr                               !<  number of radial points
  integer, intent(out)              ::  lmax                             !<  maximum angular value

! local variables

  integer                           ::  itype                            !  type of calculation (-1 signals end of calculation)
  character(len=2)                  ::  nameat                           !  chemical symbol of the atom
  character(len=2)                  ::  icorr                            !  correlation type
  character(len=1)                  ::  ispp                             !  spin polarization  ' ', 's', 'r'

  real(REAL64)                      ::  znuc                             !  nuclear charge

  integer                           ::  ncore                            !  number of orbitals treated as core

! local

  character(len=2)           ::  ctype                                   !  type of calculation flag
  integer                    ::  nval                                    !  number of valence orbitals
  integer                    ::  ni, li                                  !  principal and angular quantum number
  integer                    ::  iz

  character(len=2)           ::  icorrt, namet
  character(len=3)           ::  irel
  character(len=4)           ::  nicore
  character(len=10)          ::  iray(6), psdtitle(20)

  integer                    ::  nrm
  integer                    ::  npotd, npotu

! counters

  integer     ::  i

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  EPS = 0.00001_REAL64

  read(ioread,'(3x,a2)') ctype

!  if ctype = ' ' , no more data, program ends

  if (ctype == 'ae') then
    itype=0
  elseif (ctype == 'pg') then
    itype=1
  elseif (ctype == 'pe') then
    itype=2
  elseif (ctype == 'ph') then
    itype=3
  elseif (ctype == 'pt') then
    itype=4
  elseif (ctype == 'pm') then
    itype=5
  else
    itype=-1

    return

  endif


! whatever the ctype it is a test!!! (itype = 4)

  itype = 4

  read(ioread,*)

  read(ioread,'(3x,a2,3x,a2,a1)') nameat,icorr,ispp


  read(ioread,*) znuc

  if (abs(znuc) < EPS) then
    call atom_p_tbl_charge(nameat,iz)
    znuc = ONE*iz
  endif

  open(unit=iopsd, file=trim(filepsd), form='unformatted', status='unknown')

  read(iopsd) namet, icorrt, irel ,nicore, (iray(i),i=1,6),              &
      (psdtitle(i),i=1,7), npotd, npotu, nrm

  nr = nrm + 1

  close(unit=iopsd)

! read the number of core and valence orbitals


  read(ioread,*) ncore, nval


! compute occupation numbers and orbital quantum numbers for the valence

  lmax = 0
  norb = 0

  if (nval > 0) then

    do i = 1,nval

      read(ioread,'(2i5)') ni,li

      if(ispp == ' ') then
        norb = norb + 1
      elseif (ispp == 'r' .and. li == 0) then
        norb = norb + 1
      else
        norb = norb + 2
      endif

      if (lmax < li) lmax = li

    enddo

  else

    write(6,*) '  Stopped in atom_kb_test_input_size. nval = ',nval

    STOP

  endif

  return

end subroutine atom_kb_test_input_size
