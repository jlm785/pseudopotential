!>  Calls the atomic program
!>
!>  \author       Jose Luis Martins
!>  \version      6.013
!>  \date         14 July 2021
!>  \copyright    GNU Public License v2

program atom

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! tape stuff

  integer                           ::  iowrite                          !  default output tape

  integer                           ::  ioread                           !  default tape for reading
  character(len=8)                  ::  filein                           !  name of default tape for reading

  integer                           ::  ioae                             !  default tape for writing
  character(len=12)                 ::  fileae                           !  name of default tape for all-electron results

  integer                           ::  iopsd                            !  default tape for pseudopotential in old format
  character(len=10)                 ::  filepsd                          !  name of default tape for reading pseudopotential in old format

! local variables

  logical                           ::  lex
  character(len=2)                  ::  nameat                           !  chemical symbol of the atom


  iowrite = 6

  iopsd = 2
  filepsd = 'pseudo.dat'

  ioread = 3
  filein = 'atom.dat'

  ioae = 7
  fileae = 'datafile.dat'


  if(iowrite /= 6) then
    open(unit=iowrite,file='atom.out',form='FORMATTED',status='UNKNOWN')
  endif

  inquire(file=filein, exist = lex)

  if(lex) then
    write(6,*)
    write(6,*) '  using file ',filein
    write(6,*)
  else

    write(6,*)
    write(6,*)
    write(6,*) '  Enter chemical symbol'
    write(6,*)

    nameat = '  '
    read(5,*) nameat
    call atom_atm_write_sub(ioread, filein, nameat)
    write(6,*)
  endif

! runs the calculation

  call atom_atm_sub(iowrite, ioread, filein, ioae, fileae, iopsd, filepsd)

  stop

end program atom
