!>  Calls the atomic program
!>
!>  \author       Jose Luis Martins
!>  \version      6.2.0
!>  \date         14 July 2021. 14 November 2026.
!>  \copyright    GNU Public License v2

program atom

! choice of tested pseudopotential. 12 November 2026. JLM+claude
! choice of exchange and correlation functional. 14 November 2026. JLM+claude

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
  character(len=4)                  ::  vtype                            !  label of tested pseudopotential, 'NONE' old tables
  character(len=2)                  ::  icorr                            !  exchange and correlation functional (ca, pw, pb)


  iowrite = 6

  iopsd = 13
  filepsd = 'pseudo.dat'

  ioread = 11
  filein = 'atom.dat'

  ioae = 12
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
    call atom_atm_tested_choose(nameat, vtype)
    call atom_atm_xc_choose(icorr)
    call atom_atm_write_sub(ioread, filein, nameat, vtype, icorr)
    write(6,*)
  endif

! runs the calculation

  call atom_atm_sub(iowrite, ioread, filein, ioae, fileae, iopsd, filepsd)

  stop

end program atom
