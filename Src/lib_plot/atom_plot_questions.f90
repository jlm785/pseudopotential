!>  Asks the questions about what plots to show
!>
!>  \author       Manuel Maria Alemany, Jose Luis Martins
!>  \version      6.013
!>  \date         January 2000, 6 August 2021
!>  \copyright    GNU Public License v2

subroutine atom_plot_questions(task, x_min, x_max, y_min, y_max, lcentr, lnew)

! adapted from outdis. 9 August 2021. JLM

  implicit none


  integer, parameter          :: REAL64 = selected_real_kind(12)

!input

  character(len=5), intent(in)      ::  task                             !<  select question

! modified

  real(REAL64), intent(inout)       ::  x_min, x_max                     !<  x range
  real(REAL64), intent(inout)       ::  y_min, y_max                     !<  y range

! output

  integer, intent(out)              ::  lcentr                           !<  angular momentum
  logical, intent(out)              ::  lnew                             !<  true if the ranges changed

! local variables

  character(len=1)    ::  opt, opt1, opt2
  integer             ::  loop, loop1, loop2

  integer             ::  ioerr

  if(trim(task) == 'range') then

    do loop1 = 1,10

      write(6,*)
      write(6,*) 'Do you want to change the axis ranges? (y/n)'
      write(6,*)
      read(5,*) opt1

      if(opt1 == 'y' .or. opt1 == 'Y') then
        do loop2 = 1,10
          write(6,*)
          write(6,*) 'x-range (1), y-range (2) or both (3)'
          write(6,*)
          read(5,*) opt2
          if(opt2 == '1') then
            write(6,*) 'Values of xmin xmax'
            read(5,*,iostat=ioerr) x_min,x_max
            if(ioerr == 0) exit
          elseif(opt2 == '2') then
            write(6,*) 'Values of ymin ymax'
            read(5,*,iostat=ioerr) y_min,y_max
            if(ioerr == 0) exit
          elseif(opt2 == '3') then
            write(6,*) 'Values of xmin xmax'
            read(5,*,iostat=ioerr) x_min,x_max
            if(ioerr == 0) then
              write(6,*) 'Values of ymin ymax'
              read(5,*,iostat=ioerr) y_min,y_max
              if(ioerr == 0) exit
            endif
          else
            write(6,*) ' Valid input is two real numbers'
            ioerr = 1
            if(loop2 == 10) write(6,*) ' You had 10 chances...'
          endif
        enddo
        if(ioerr == 0) then
          lnew = .TRUE.
          exit
        else
          lnew = .FALSE.
        endif
      elseif(opt1 == 'n' .or. opt1 == 'N') then
        lnew = .FALSE.
        exit
      else
        write(6,*) ' Valid input is y or n'
        lnew = .FALSE.
        if(loop1 == 10) write(6,*) ' You had 10 chances...'
      endif
    enddo

  elseif(trim(task) == 'pdf') then

    do loop = 1,10

      write(6,*)
      write(6,*) 'Do you want to save the plot in a PDF file? (y/n)'
      write(6,*)
      read(5,*) opt
      if(opt == 'y' .or. opt == 'Y') then
        lnew = .TRUE.
        exit
      elseif(opt == 'n' .or. opt == 'N') then
        lnew = .FALSE.
        exit
      else
        write(6,*) ' Valid input is y or n'
        lnew = .FALSE.
        if(loop == 10) write(6,*) ' You had 10 chances...'
      endif
    enddo

   elseif(trim(task) == 'clean') then

    do loop = 1,10

      write(6,*)
      write(6,*) 'Do you want to remove the *.gp gnuplot files? (y/n)'
      write(6,*)
      read(5,*) opt
      if(opt == 'y' .or. opt == 'Y') then
        lnew = .TRUE.
        exit
      elseif(opt == 'n' .or. opt == 'N') then
        lnew = .FALSE.
        exit
      else
        write(6,*) ' Valid input is y or n'
        lnew = .FALSE.
        if(loop == 10) write(6,*) ' You had 10 chances...'
      endif
    enddo

  elseif(trim(task) == 'centr') then

    do loop1 = 1,10

      write(6,*)
      write(6,*) 'Do you want to plot a centrifugal potential? (y/n)'
      write(6,*)
      read(5,*) opt
      if(opt == 'y' .or. opt == 'Y') then

        do loop2 = 1,10
          write(6,*)
          write(6,*) 'For which angular momentum? (0/1/2/...)'
          write(6,*)
          read(5,*,iostat=ioerr) lcentr
          if(ioerr == 0 .and. lcentr > -1 .and. lcentr < 8) then
            exit
          else
            write(6,*) ' Valid input is an integer between 0 and 8'
            ioerr = 1
            if(loop2 == 10) write(6,*) ' You had 10 chances...'
          endif
        enddo
        if(ioerr == 0) then
          lnew = .TRUE.
          exit
        else
          lnew = .FALSE.
        endif
      elseif(opt == 'n' .or. opt == 'N') then
        lnew = .FALSE.
        exit
      else
        write(6,*) ' Valid input is y or n'
        lnew = .FALSE.
        if(loop1 == 10) write(6,*) ' You had 10 chances...'
      endif
    enddo

  else

    write(6,*)
    write(6,*) '  W A R N I N G   in atom_plot_questions:  unknown question ',task
    write(6,*) '  Doing nothing'
    write(6,*)

    lnew = .FALSE.

  endif

  return

end subroutine atom_plot_questions
