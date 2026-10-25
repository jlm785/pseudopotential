!>  Writes a plain text description of a figure (fig_sfname.dat)
!>  to be plotted with matplotlib by Tools/atom_plot_mpl.py.
!>  It mirrors the gnuplot commands of atom_plot_command.
!>
!>  \author       claude+JlMartins
!>  \version      6.1.0
!>  \date         25 October 2026.
!>  \copyright    GNU Public License v2

subroutine atom_plot_spec(iotape, sfname, nlabel, zion, lcentr,         &
        x_min, x_max, y_min, y_max, iw, label, ip, label2)

! Written by claude from atom_plot_command. 25 October 2026. claude

! Same arguments as atom_plot_command (without lpdf).
! Only wfct, ppot, bas* and log*/ln* are implemented, other figures
! are silently skipped.
!
! Format of the file, one item per line, "key = value", fields of
! curves separated by "|":
!   title  = text
!   xlabel = text
!   ylabel = text
!   xrange = xmin xmax
!   yrange = ymin ymax
!   curve  = file | x column | y column | label | solid/dashed | color index
!   coulomb     = c   | label | style | color index      (y = -c / x)
!   centrifugal = c   | label | style | color index      (y = c / x**2)
!   hline  = y
! Lines starting with # are comments.
! Color index 0 is used for reference curves (Coulomb tail, pseudo
! wave-function in the basis plots).

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  iotape                           !<  tape number
  character(len=4), intent(in)      ::  sfname                           !<  stem of file name

  integer, intent(in)               ::  nlabel                           !<  dimension of labels

  real(REAL64), intent(in)          ::  zion                             !<  ionic charge or log derivative radius

  integer, intent(in)               ::  lcentr                           !<  plot centrifugal potential for l = lcentr

  real(REAL64), intent(in)          ::  x_min, x_max                     !<  x range
  real(REAL64), intent(in)          ::  y_min, y_max                     !<  y range

  integer, intent(in)               ::  iw                               !<  number of all-electron wave-functions
  integer, intent(in)               ::  ip                               !<  number of pseudo wave-functions

  character(len=12), intent(in)     ::  label(nlabel)                    !<  first label of potentials
  character(len=12), intent(in)     ::  label2(nlabel)                   !<  second label

! local variables

  character(len=12)      ::  fname
  logical                ::  lknown

! counter

  integer         ::  i

  lknown = sfname == 'wfct' .or. sfname == 'ppot' .or.                   &
           sfname(1:3) == 'bas' .or. sfname(1:3) == 'log' .or.           &
           sfname(1:2) == 'ln'

  if(.not. lknown) return

  fname = 'fig_'//sfname//'.dat'

  open(unit=iotape, file=trim(fname), status='unknown', form='formatted')

  write(iotape,'("# figure ",a4," for Tools/atom_plot_mpl.py")') sfname
  write(iotape,'("xrange = ",2es16.8)') x_min, x_max
  write(iotape,'("yrange = ",2es16.8)') y_min, y_max

  if(sfname == 'wfct') then

    write(iotape,'("title = Wave Functions")')
    write(iotape,'("xlabel = r [a.u.]")')
    write(iotape,'("ylabel = rR(r)")')

    do i = 1,iw
      write(iotape,'("curve = ",a4,".gp | 1 | ",i3," | ",a," | dashed | ",i3)') &
          sfname, 2*i, trim(label(i)), i
    enddo
    do i = 1,ip
      write(iotape,'("curve = ",a4,".gp | 1 | ",i3," | ",a," | solid | ",i3)')  &
          sfname, 2*i+1, trim(label2(i)), i
    enddo

  elseif(sfname == 'ppot') then

    write(iotape,'("title = Pseudopotentials")')
    write(iotape,'("xlabel = r [a.u.]")')
    write(iotape,'("ylabel = V(r) [Ry]")')

    do i = 1,iw
      write(iotape,'("curve = ",a4,".gp | 1 | ",i3," | ",a," | solid | ",i3)')  &
          sfname, i+1, trim(label(i)), i
    enddo

    if(lcentr > 0) then
      write(iotape,'("centrifugal = ",i4," | centrif. pot. | dashed | 0")')    &
          lcentr*(lcentr+1)
    endif

    write(iotape,'("coulomb = ",es16.8," | Z - ion | dashed | 0")') 2*zion

  elseif(sfname(1:3) == 'bas') then

    write(iotape,'("title = Wave and Basis Functions for l = ",a1)') sfname(4:4)
    write(iotape,'("xlabel = r [a.u.]")')
    write(iotape,'("ylabel = rR(r)")')

    if(iw /= 0) then
      write(iotape,'("curve = ",a4,".gp | 1 | 2 | ",a," | dashed | 0")')        &
          sfname, trim(label(1))
    endif
    do i = 1,ip
      write(iotape,'("curve = ",a4,".gp | 1 | ",i3," | ",a," | solid | ",i3)')  &
          sfname, i+iw+1, trim(label2(i)), i
    enddo

  else

!   log derivatives, zion is the radius

    write(iotape,'("title = Logarithmic Derivatives at ",f5.2," (a.u.)")') zion
    write(iotape,'("xlabel = Energy (Ry)")')
    write(iotape,'("ylabel = Log. Derivatives (1/a.u.)")')

    write(iotape,'("curve = log_true.gp | 1 | ",i3," | AE ",a," | solid | 2")')        &
        ip+1, trim(label(ip))
    write(iotape,'("curve = log_ppot.gp | 1 | ",i3," | SL pseudo ",a," | solid | 3")') &
        ip+1, trim(label(ip))
    if(iw == 1) then
      write(iotape,'("curve = log_ppkb.gp | 1 | ",i3," | KB pseudo ",a," | dashed | 4")') &
          ip+1, trim(label(ip))
    endif

  endif

  write(iotape,'("hline = 0.0")')

  close(iotape)

  return

end subroutine atom_plot_spec
