program plot_show

  implicit none

  integer                           ::  iowrite                          !  default tape for writing

  integer                           ::  ioplot                           !  tape to read input
  character(len=10)                 ::  fname                            !  name of file with plot information

  integer                           ::  iost                             !  data files use tapes iost+1 to iost+4
  integer                           ::  iocomm                           !  tape number for gnuplot command file
  integer                           ::  iotmp                            !  tape number for temporary files

  logical                           ::  lint                             !  interactive run

  lint = .TRUE.
!  lint = .FALSE.

  iowrite = 6

  ioplot = 8
  iost = 30
  iocomm = 21
  iotmp = 23

! asking for outputfile to be used

  fname = 'plot.dat  '

  if(lint) then
    write(6,*)
    write(6,*) ' Which file shall be used to read the plot information?'
    write(6,*) ' Please give the name. (format = a10)'
    write(6,*) ' Press <ENTER> to choose the default <plot.dat>'
    write(6,*)

    read(5,'(a10)') fname
    if(fname == ' ') fname = 'plot.dat  '

  endif

  call atom_plot_sub(iowrite, ioplot, fname, iost, iocomm, iotmp, lint)

  stop

end program plot_show


