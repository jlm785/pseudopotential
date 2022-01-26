!>  Calls the atomic program
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.6
!>  \date         28 October 2021
!>  \copyright    GNU Public License v2

program atom_test

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! tape stuff

  integer                           ::  iowrite                          !  default output tape

  integer                           ::  ioread                           !  default tape for reading
  character(len=8)                  ::  filein                           !  name of default tape for reading

  integer                           ::  iopsd                            !  default tape for pseudopotential in old format
  character(len=10)                 ::  filepsd                          !  name of default tape for reading pseudopotential in old format

  integer                           ::  iopsdkb                          !  default tape for KB pseudopotential in real space
  character(len=15)                 ::  filepsdkb                        !  name of default tape for writeing KB pseudopotential in real space


  iowrite = 6

  iopsd = 2
  filepsd = 'pseudo.dat'

  ioread = 3
  filein = 'atom.dat'

  iopsdkb = 10
  filepsdkb = 'pseudokb.dat   '

  if(iowrite /= 6) then
    open(unit=iowrite,file='atom.out',form='FORMATTED',status='UNKNOWN')
  endif

! runs the test

  call atom_atm_test(iowrite, ioread, filein,                            &
         iopsd, filepsd, iopsdkb, filepsdkb)

  stop

end program atom_test
