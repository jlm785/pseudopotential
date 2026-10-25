!>  Manages the temporary files (tmp_mmga_*) associated with gnuplot
!>
!>  \author       Manuel Maria Alemany, Jose Luis Martins
!>  \version      6.1.0
!>  \date         January 2000, 25 October 2026
!>  \copyright    GNU Public License v2

subroutine atom_plot_gnuplot(iotmp, sfname, task)

! adapted from outdis. 9 August 2021. JLM
! removes fig_*.dat. 25 October 2026. claude
! kill: own processes (pgrep), no crash if none, pid longer than 6. 25 October 2026. claude

  implicit none

! input

  integer, intent(in)               ::  iotmp                            !<  tape number
  character(len=4), intent(in)      ::  sfname                           !<  stem of file name
  character(len=5), intent(in)      ::  task                             !<  task to perform

! local variables

  character(len=15)   :: fname

  character(len=12)   ::  c_pid
  integer             ::  ios

  fname = 'command_'//sfname//'.gp'

  if(trim(task) == 'run') then

    open(unit=iotmp, file='tmp_mmga_x', status='unknown', form='formatted')
    write(iotmp,'("gnuplot  ",a16," 1> /dev/null 2>&1 &")') fname
    close(iotmp)

    call system ('chmod +x ./tmp_mmga_x')
    call system ('./tmp_mmga_x')

  elseif(trim(task) == 'kill') then

!   most recent gnuplot of the user, nothing to do if there is none

    call system('pgrep -n -u "$(id -u)" gnuplot > tmp_mmga_t 2> /dev/null')

    open(unit=iotmp, file='tmp_mmga_t', status='unknown', form='formatted')

    c_pid = ' '
    read(iotmp,*,iostat=ios) c_pid

    close(iotmp)

    if(ios == 0 .and. len_trim(c_pid) > 0) then

      open(unit=iotmp, file='tmp_mmga_x', status='unknown', form='formatted')

      write(iotmp,'("kill -9 ",a," 1> /dev/null 2>&1  ")') trim(c_pid)

      close(iotmp)

      call system ('chmod +x ./tmp_mmga_x')
      call system ('./tmp_mmga_x')

    endif

  elseif(trim(task) == 'clean') then

    call system('rm ./tmp_mmga_x ./tmp_mmga_t 2> /dev/null')

  elseif(trim(task) == 'files') then

    call system('rm ./command_*.gp 2> /dev/null')
    call system('rm ./*t.gp ./log*.gp 2> /dev/null')
    call system('rm ./kb*.gp ./bas*.gp 2> /dev/null')
    call system('rm ./fig_*.dat 2> /dev/null')

  else

    write(6,*)
    write(6,*) '  W A R N I N G   in atom_plot_gnuplot:  unknown task ',task
    write(6,*) '  Doing nothing'
    write(6,*)

  endif

  return

end subroutine atom_plot_gnuplot
