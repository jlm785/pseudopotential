!>  Manages the temporary files (tmp_mmga_*) associated with gnuplot
!>
!>  \author       Manuel Maria Alemany, Jose Luis Martins
!>  \version      6.1.0
!>  \date         January 2000, 25 October 2026, 30 October 2026
!>  \copyright    GNU Public License v2

subroutine atom_plot_gnuplot(iotmp, sfname, task)

! adapted from outdis. 9 August 2021. JLM
! removes fig_*.dat. 25 October 2026. claude
! kill: own processes (pgrep), no crash if none, pid longer than 6. 25 October 2026. claude
! run: nothing (one warning) if gnuplot is not installed. 29 October 2026. JLM+claude
! call system replaced by the standard execute_command_line. 30 October 2026. JLM+claude

  implicit none

! input

  integer, intent(in)               ::  iotmp                            !<  tape number
  character(len=4), intent(in)      ::  sfname                           !<  stem of file name
  character(len=5), intent(in)      ::  task                             !<  task to perform

! local variables

  character(len=15)   :: fname

  character(len=12)   ::  c_pid
  integer             ::  ios
  integer             ::  iexit, icmd                      !  exit and command status of execute_command_line (not used)

! external function

  logical, external   ::  atom_plot_gnuplot_ok

  fname = 'command_'//sfname//'.gp'

  if(trim(task) == 'run') then

    if(.not. atom_plot_gnuplot_ok(iotmp)) return

    open(unit=iotmp, file='tmp_mmga_x', status='unknown', form='formatted')
    write(iotmp,'("gnuplot  ",a16," 1> /dev/null 2>&1 &")') fname
    close(iotmp)

    call execute_command_line('chmod +x ./tmp_mmga_x',                     &
        exitstat=iexit, cmdstat=icmd)
    call execute_command_line('./tmp_mmga_x',                              &
        exitstat=iexit, cmdstat=icmd)

  elseif(trim(task) == 'kill') then

!   most recent gnuplot of the user, nothing to do if there is none

    call execute_command_line('pgrep -n -u "$(id -u)" gnuplot > tmp_mmga_t 2> /dev/null', &
        exitstat=iexit, cmdstat=icmd)

    open(unit=iotmp, file='tmp_mmga_t', status='unknown', form='formatted')

    c_pid = ' '
    read(iotmp,*,iostat=ios) c_pid

    close(iotmp)

    if(ios == 0 .and. len_trim(c_pid) > 0) then

      open(unit=iotmp, file='tmp_mmga_x', status='unknown', form='formatted')

      write(iotmp,'("kill -9 ",a," 1> /dev/null 2>&1  ")') trim(c_pid)

      close(iotmp)

      call execute_command_line('chmod +x ./tmp_mmga_x',                   &
          exitstat=iexit, cmdstat=icmd)
      call execute_command_line('./tmp_mmga_x',                            &
          exitstat=iexit, cmdstat=icmd)

    endif

  elseif(trim(task) == 'clean') then

    call execute_command_line('rm ./tmp_mmga_x ./tmp_mmga_t 2> /dev/null', &
        exitstat=iexit, cmdstat=icmd)

  elseif(trim(task) == 'files') then

    call execute_command_line('rm ./command_*.gp 2> /dev/null',            &
        exitstat=iexit, cmdstat=icmd)
    call execute_command_line('rm ./*t.gp ./log*.gp 2> /dev/null',         &
        exitstat=iexit, cmdstat=icmd)
    call execute_command_line('rm ./kb*.gp ./bas*.gp 2> /dev/null',        &
        exitstat=iexit, cmdstat=icmd)
    call execute_command_line('rm ./fig_*.dat 2> /dev/null',               &
        exitstat=iexit, cmdstat=icmd)

  else

    write(6,*)
    write(6,*) '  W A R N I N G   in atom_plot_gnuplot:  unknown task ',task
    write(6,*) '  Doing nothing'
    write(6,*)

  endif

  return

end subroutine atom_plot_gnuplot
