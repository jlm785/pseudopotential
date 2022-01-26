program plot_kb_show

  implicit none

  integer                           ::  iowrite                          !  default tape for writing

  integer                           ::  ioplot                           !  tape to read input
  character(len=15)                 ::  fname                            !  name of file with plot information

  integer                           ::  iost                             !  data files use tapes iost+1 to iost+4
  integer                           ::  iocomm                           !  tape number for gnuplot command file
  integer                           ::  iotmp                            !  tape number for temporary files

  logical                           ::  lint                             !  interactive run

  lint = .TRUE.
!  lint = .FALSE.

  iowrite = 6

  ioplot = 8
  iost = 10
  iocomm = 21
  iotmp = 23

! asking for outputfile to be used

  fname = 'plot_kb.dat    '

  if(lint) then
    write(6,*)
    write(6,*) ' Which file shall be used to read the plot information?'
    write(6,*) ' Please give the name. (format = a15)'
    write(6,*) ' Press <ENTER> to choose the default <plot_kb.dat>'
    write(6,*)

    read(5,'(a15)') fname
    if(fname == ' ') fname = 'plot_kb.dat    '

  endif

  call atom_plot_kb_sub(iowrite, ioplot, fname, iost, iocomm, iotmp, lint)

  stop

end program plot_kb_show


