program kb_test

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! tape stuff

  integer                           ::  iowrite                          !  default output tape

  integer                           ::  ioread                           !  default tape for reading
  character(len=8)                  ::  filein                           !  name of default tape for reading

  integer                           ::  iopsd                            !  default tape for pseudopotential in old format
  character(len=12)                 ::  filepsd                          !  name of default tape for reading pseudopotential in old format


  iowrite = 6

  ioread = 3
  filein = 'atom.dat'

  iopsd = 2
  filepsd = 'pseudokb.dat'

  if(iowrite /= 6) then
    open(unit=iowrite,file='atom.out',form='FORMATTED',status='UNKNOWN')
  endif

! runs the calculation

  call atom_kb_test_sub(iowrite, ioread, filein, iopsd, filepsd)

end program kb_test
