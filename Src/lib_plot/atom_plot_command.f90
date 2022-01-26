!>  Writes the files with gnuplot commands
!>
!>  \author       Peter Schuster, Manuel Maria Alemany, Jose Luis Martins
!>  \version      6.0.4
!>  \date         April 1993, January 2000, 9 August 2021, 9 October 2021.
!>  \copyright    GNU Public License v2

subroutine atom_plot_command(iotape, sfname, nlabel, zion, lcentr, lpdf, &
        x_min, x_max, y_min, y_max, iw, label, ip, label2)

! extracted from outdis. 9 August 2021. JLM
! kinetic, 9 October 2021. JLM
! offset removed. 18 October 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  iotape                           !<  tape number
  character(len=4), intent(in)      ::  sfname                           !<  stem of file name

  integer, intent(in)               ::  nlabel                           !<  dimension of labels

  real(REAL64), intent(in)          ::  zion                             !<  ionic charge or log derivative radius

  logical, intent(in)               ::  lpdf                             !<  outputs a pdf file

  integer, intent(in)               ::  lcentr                           !<  plot centrifugal potential for l = lcentr

  real(REAL64), intent(in)          ::  x_min, x_max                     !<  x range
  real(REAL64), intent(in)          ::  y_min, y_max                     !<  y range

  integer, intent(in)               ::  iw                               !<  number of all-electron wave-functions
  integer, intent(in)               ::  ip                               !<  number of pseudo wave-functions

  character(len=12), intent(in)     ::  label(nlabel)                    !<  first label of potentials
  character(len=12), intent(in)     ::  label2(nlabel)                   !<  second label

! local variables

  character(len=15)      ::  fname

! counter

  integer         ::  i

  fname = 'command_'//sfname//'.gp'

  open(unit=iotape, file=trim(fname), status='unknown', form='formatted')

  if(lpdf) then
    write(iotape,'("set term pdf color font ''Helvetica,12''")')
    write(iotape,'("set output ''",a4,".pdf''")') sfname
  else
    write(iotape,'("set term wxt persist")')
  endif

  if(sfname /= 'kbek') then

    write(iotape,'("unset label")')
    write(iotape,'("set key")')

    write(iotape,'("set format x ''%.1f''")')
    write(iotape,'("set format y ''%.1f''")')

    write(iotape,'("set yrange [",f5.1,":",f5.1,"]")')  y_min,y_max
    write(iotape,'("set xrange [",f5.1,":",f5.1,"]")')  x_min,x_max

  endif

  if(sfname == 'wfct') then

    write(iotape,'("set title ''Wave Functions''")')

    write(iotape,'("set xlabel ''r [a.u.]''")')
    write(iotape,'("set ylabel ''rR(r)'' ")')


    write(iotape,'("plot \")')

    do i = 2,2*iw,2
      write(iotape,'("''",a4,".gp'' using 1:",i2," ti ''",a12,"''  w li ls ",i2," dashtype 2, \")') sfname, i, label(i/2), i/2
    enddo

    do i = 3,2*ip+1,2
      write(iotape,'("''",a4,".gp'' using 1:",i2," ti ''",a12,"'' w li ls ",i2,", \")') sfname, i, label2((i-1)/2), (i-1)/2
    enddo

  elseif(sfname == 'fwft') then

    write(iotape,'("set title ''Wave Function Transforms''")')

!    write(iotape,'("unset key")')
    write(iotape,'("set format y ''%.2f''")')
    write(iotape,'("set xlabel ''q [1/a.u.]''")')
    write(iotape,'("set ylabel ''R(q)'' ")')

    write(iotape,'("plot \")')

    do i = 2,iw+1
      write(iotape,'("''",a4,".gp'' using 1:",i2," ti ''",a12,"'' with lines ls ",i2,", \")') sfname, i, label(i-1),i-1
    enddo

  elseif(sfname == 'ppot') then

    write(iotape,'("set title ''Pseudopotentials''")')

    write(iotape,'("set xlabel ''r [a.u.]''")')
    write(iotape,'("set ylabel ''V(r) [Ry]'' ")')

    write(iotape,'("plot \")')

    do i = 2,iw+1
       write(iotape,'("''",a4,".gp'' using 1:",i2," ti ''",a12,"'' w li ls ",i2,", \")') &
            sfname, i, label(i-1), i-1
    enddo

    if(lcentr > 0) then
      write(iotape,'(i3,"/(x*x)  ti ''cetrif. pot.'' w li dashtype 2, \")')  &
                 lcentr*(lcentr+1)
    endif

    write(iotape,'("-",f5.2,"/x"," ti ''Z - ion'' w li ls 5, \")') zion*2.0

  elseif(sfname == 'fppt') then

    write(iotape,'("set title ''Fourier Transforms of Potential''")')

    write(iotape,'("set size 2.8/5.,3/3.")')

    write(iotape,'("set xlabel ''q [1/a.u.]''")')
    write(iotape,'("set ylabel ''q^2/4piZion V(q) [Ry]'' ")')

    write(iotape,'("plot \")')
    do i = 2,iw+1
      write(iotape,'("''",a4,".gp'' using 1:",i2," ti ''",a12,"'' w li ls ",i2,", \")') &
             sfname, i, label(i-1), i-1
    enddo

  elseif(sfname(1:3) == 'log' .or. sfname(1:2) == 'ln') then

    write(iotape,'("set title ''Logarithmic Derivatives at",f5.2," (a.u.)''")') zion

    write(iotape,'("set xlabel ''Energy (Ry)''")')
    write(iotape,'("set ylabel ''Log. Derivatives (1/a.u.)'' ")')

    write(iotape,'("plot \")')
    write(iotape,'("''log_true.gp'' using 1:",i2," ti ''AE ",            &
       &    a12,"'' w li ls ",i2,", \")')   ip+1, label(ip), 2
    write(iotape,'("''log_ppot.gp'' using 1:",i2," ti ''SL pseudo ",     &
       &    a12,"'' w li ls ",i2,", \")')   ip+1, label(ip), 3
    if(iw == 1) then
      write(iotape,'("''log_ppkb.gp'' using 1:",i2," ti ''KB pseudo ",   &
       &    a12,"'' w li ls ",i2,", \")')   ip+1, label(ip), 4
    endif

  elseif(sfname(1:3) == 'bas') then

    write(iotape,'("set title ''Wave and Basis Functions for l = ",a1,"''")') sfname(4:4)

    write(iotape,'("set xlabel ''r [a.u.]''")')
    write(iotape,'("set ylabel ''rR(r)'' ")')


    write(iotape,'("plot \")')

    if(iw /= 0) then
      write(iotape,'("''",a4,".gp'' using 1:2 ti ''",a12,"''  w li ls 2 dashtype 2, \")') sfname, label(1)
    endif

    do i = 1,ip
      write(iotape,'("''",a4,".gp'' using 1:",i1," ti ''",a12,"'' w li ls ",i2,", \")') sfname, i+iw+1, label2(i), i
    enddo

  elseif(sfname(1:4) == 'kbop') then

    write(iotape,'("set title '' Kleinman-Bylander Operators '' ")')

    write(iotape,'("set xlabel ''r [a.u.]''")')
    write(iotape,'("set ylabel ''P(r)'' ")')

    write(iotape,'("plot \")')
    do i = 2,iw+1
      write(iotape,'("''",a4,".gp'' using 1:",i2," ti ''",a12,"'' w li ls ",i2,", \")') &
             sfname, i, label(i-1), i-1
    enddo

  elseif(sfname(1:4) == 'kbtr') then

    write(iotape,'("set title '' Fourier Transform Kleinman-Bylander Operators '' ")')

    write(iotape,'("set xlabel ''q [1 / a.u.]''")')
    write(iotape,'("set ylabel ''P(q)'' ")')

    write(iotape,'("plot \")')
    do i = 2,iw+1
      write(iotape,'("''",a4,".gp'' using 1:",i2," ti ''",a12,"'' w li ls ",i2,", \")') &
             sfname, i, label(i-1), i-1
    enddo

  elseif(sfname == 'kbek') then

    write(iotape,'("unset label")')
    write(iotape,'("set key")')

    write(iotape,'("set format x ''%.1f''")')
    write(iotape,'("set format y ''%.1e''")')

    write(iotape,'("set yrange [",f18.12,":",f10.2,"]")')  y_min,y_max
    write(iotape,'("set xrange [",f5.1,":",f5.1,"]")')  x_min,x_max

    write(iotape,'("set title ''Integral of Kinetic energy (difference)''")')

    write(iotape,'("set xlabel ''q [1/a.u.]''")')
    write(iotape,'("set ylabel ''E_kin(inf) - E_kin(q) [Ry]'' ")')
    write(iotape,'("set log y")')

    write(iotape,'("plot \")')
    write(iotape,'("''",a4,".gp'' using 1:2 ti ''",a12,"'' w li ls 1, \")') &
             sfname, label(1)

  elseif(sfname == 'kbvk') then

    write(iotape,'("set title ''Choice of Local Pseudopotential''")')

    write(iotape,'("set xlabel ''r [a.u.]''")')
    write(iotape,'("set ylabel ''V_local(r) [Ry]'' ")')

    write(iotape,'("plot \")')
    write(iotape,'("''",a4,".gp'' using 1:2 ti ''",a12,"'' w li ls 1, \")') &
             sfname, label(1)

    else

    write(6,*) '  stopped in atom_plot_command'
    write(6,*) '  unknown type of file:  ',sfname

    STOP

  endif

  if(sfname /= 'kbek') then

    write(iotape,'("0 w li ls 8  notitle")')

  else

    write(iotape,'("0.0004 w li ls 6  notitle, \")')
  write(iotape,'("0.005 w li ls 6  notitle")')

  endif

  close(iotape)

  return

end subroutine atom_plot_command
