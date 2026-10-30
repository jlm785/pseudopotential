!>  Checks (once) if gnuplot is available in the PATH.
!>  If not, writes a warning the first time it is called.
!>
!>  \author       claude
!>  \version      6.1.0
!>  \date         29 October 2026, 30 October 2026
!>  \copyright    GNU Public License v2

logical function atom_plot_gnuplot_ok(iotmp)

! written 29 October 2026. claude
! call system replaced by the standard execute_command_line. 30 October 2026. JLM+claude

  implicit none

! input

  integer, intent(in)               ::  iotmp                            !<  tape number

! local variables

  logical, save       ::  lfirst = .TRUE.
  logical, save       ::  lok = .FALSE.

  character(len=12)   ::  cpath
  integer             ::  ios
  integer             ::  iexit, icmd                      !  exit and command status of execute_command_line (not used)

  if(lfirst) then

    lfirst = .FALSE.

    call execute_command_line('command -v gnuplot > tmp_mmga_c 2> /dev/null', &
        exitstat=iexit, cmdstat=icmd)

    open(unit=iotmp, file='tmp_mmga_c', status='unknown', form='formatted')

    cpath = ' '
    read(iotmp,'(a12)',iostat=ios) cpath

    close(iotmp, status='delete')

    lok = (ios == 0 .and. len_trim(cpath) > 0)

    if(.not. lok) then
      write(6,*)
      write(6,*) '  W A R N I N G :  gnuplot was not found.'
      write(6,*) '  The plots are not displayed and the PDF files are not written.'
      write(6,*) '  The gnuplot command files (command_*.gp) and the data files'
      write(6,*) '  are written.  Tools/atom_plot_mpl.py draws the figures from'
      write(6,*) '  the fig_*.dat files with matplotlib.'
      write(6,*)
    endif

  endif

  atom_plot_gnuplot_ok = lok

  return

end function atom_plot_gnuplot_ok
