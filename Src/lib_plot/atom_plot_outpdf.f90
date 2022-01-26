!>  Writes the command-file to be accepted by gnuplot
!>  and gets the plots in pdf format
!>
!>  \author       Peter Schuster, Manuel Maria Alemany, Jose Luis Martins
!>  \version      6.013
!>  \date         April 1993, January 2000, 30 August 2021
!>  \copyright    GNU Public License v2

subroutine atom_plot_outpdf(iocomm, iotmp, numw, nump, zion,             &
      iw, ip, iv, ifp, ifw,                                              &
      labelw, labelp, labelv, labelf, labelfw, yminwp, ymaxwp, yminv)

! adapted from old code
! originally written by Peter Schuster, April 1993
! interative ploting added by Manuel Maria Alemany, January 2000
! removed H formats 22 January 2008. JLM
! modified June 24 2012 for new gnuplot conventions. JLM
! converted to fortran 90, 6 August 2021


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

!  input

  integer, intent(in)               ::  iocomm                           !<  tape number for gnuplot command file
  integer, intent(in)               ::  iotmp                            !<  tape number for temporary files

  integer, intent(in)               ::  numw                             !<  dimension of number of true wave-functions
  integer, intent(in)               ::  nump                             !<  dimension of number of pseudo-wave-functions

  real(REAL64), intent(in)          ::  zion                             !<  ionic charge

  integer, intent(inout)            ::  iw                               !<  number of true wave-functions
  character(len=12), intent(in)     ::  labelw(numw)                     !<  label of true wave-functions

  integer, intent(in)               ::  ip                               !<  number of pseudo-wave-functions
  character(len=12), intent(in)     ::  labelp(numw)                     !<  label of pseudo-wave-functions

  integer, intent(in)               ::  iv                               !<  number of pseudo-potentials
  character(len=12), intent(in)     ::  labelv(nump)                     !<  label of pseudo-potentials

  integer, intent(in)               ::  ifp                              !<  number of pseudo-potentials transforms
  character(len=12), intent(in)     ::  labelf(nump)                     !<  label of pseudo-potentials transforms

  integer , intent(in)              ::  ifw                              !<  number of pseudo-wave-function transforms
  character(len=12), intent(in)     ::  labelfw(nump)                    !<  label of pseudo-wave-function transforms

  real(REAL64), intent(in)          ::  yminwp, ymaxwp, yminv            !<  min and max

! local variables

  real(REAL64)          ::  x_min, x_max                                 !  x range
  real(REAL64)          ::  y_min, y_max                                 !  y range

  character(len=4)      ::  sfname

! constants

  real(REAL64), parameter   ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  real(REAL64), parameter   ::  HALF = 0.5_REAL64


! create first sheet: true and pseudo wave-functions

  x_min = ZERO
  x_max = 5*ONE + HALF

  y_min = HALF*floor(2*yminwp)
  y_max = HALF*(floor(2*ymaxwp)+1)

  sfname = 'wfct'

  call atom_plot_command(iocomm, sfname, numw, zion, 0, .TRUE.,          &
      x_min, x_max, y_min, y_max, iw, labelw, ip, labelp)

  call atom_plot_gnuplot(iotmp, sfname, 'run  ')

! create second sheet: pseudo wave-function transforms

  x_min = ZERO
  x_max = 10*ONE

  y_min = -3*HALF
  y_max =  3*HALF

  sfname = 'fwft'

  call atom_plot_command(iocomm, sfname, nump, zion, 0, .TRUE.,          &
      x_min, x_max, y_min, y_max, ifw, labelfw, ip, labelfw)

  call atom_plot_gnuplot(iotmp, sfname, 'run  ')

! create third sheet: pseudopotentials

  x_min = ZERO
  x_max = 4*ONE

  y_min = 5*(floor(yminv/5.0))*ONE
  y_min = max(y_min,-50*ONE)
  y_max = 5*ONE

  sfname = 'ppot'

  call atom_plot_command(iocomm, sfname, nump, zion, 0, .TRUE.,          &
      x_min, x_max, y_min, y_max, iv, labelv, ip, labelv)

  call atom_plot_gnuplot(iotmp, sfname, 'run  ')

!     create fourth sheet: fourier transforms of pseudopotentials

  x_min = ZERO
  x_max = 12*ONE

  y_min = -3*HALF
  y_max =  3*HALF

  sfname = 'fppt'

  call atom_plot_command(iocomm, sfname, nump, zion, 0, .TRUE.,          &
      x_min, x_max, y_min, y_max, ifp, labelf, ip, labelw)

  call atom_plot_gnuplot(iotmp, sfname, 'run  ')

! cleanup

  call atom_plot_gnuplot(iotmp, sfname, 'clean')

!  call atom_plot_gnuplot(iotmp, sfname, 'files')

  return

end subroutine atom_plot_outpdf
