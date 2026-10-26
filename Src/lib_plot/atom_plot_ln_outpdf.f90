!>  Uses gnuplot to create pdf plots of the log derivatives.
!>  Prepares files for matplotlib
!>
!>  \author       Norm Troullier, Manuel Maria Alemany, J.L.Martins
!>  \version      6.1.0
!>  \date         1990s, 30 August 2021, 25 October 2026.
!>  \copyright    GNU Public License v2

subroutine atom_plot_ln_outpdf(ioplot, iotmp, lkb,                       &
       numtotal, evlow, evhi, rpoint,                                    &
       c_val, nsc, lo, iso,  ehist, dlog_t, dlog_p, dlog_k,              &
       iowrite, mxdsc, mxdpts)


! adapted from the old program August 2021. JLM
! mxdsc. 25 September 2021. JLM
! label format, figure description for matplotlib. 25 October 2026. claude
! files written by atom_plot_ln_write, failed points as NaN. 25 October 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing
  integer, intent(in)               ::  ioplot                           !<  default tape for plots
  integer, intent(in)               ::  iotmp                            !<  default tape for temporary files

  logical, intent(in)               ::  lkb                              !<  KB pseudopotential or semi-local

  integer, intent(in)               ::  mxdsc                            !<  dimension of maximum scattering channels
  integer, intent(in)               ::  mxdpts                           !<  dimension of number of histogram points

  integer, intent(in)               ::  numtotal                         !<  number of energy in histogram
  real(REAL64), intent(in)          ::  evlow, evhi                      !<  energy range for plot
  real(REAL64), intent(in)          ::  rpoint                           !<  radius for calculating log derivatives

  character(len=4), intent(in)      ::  c_val(mxdsc)                     !<  valence identifiers
  integer, intent(in)               ::  nsc                              !<  number of scattering channels
  integer, intent(in)               ::  lo(mxdsc)                        !<  angular quantum number l
  integer, intent(in)               ::  iso(mxdsc)                       !<  2*spin or 2*(j-l)

  real(REAL64), intent(in)          ::  ehist(mxdpts)                    !<  energy
  real(REAL64), intent(in)          ::  dlog_t(mxdpts,mxdsc)             !<  log derivatives all-electron
  real(REAL64), intent(in)          ::  dlog_p(mxdpts,mxdsc)             !<  log derivatives semi-local pseudopotential
  real(REAL64), intent(in)          ::  dlog_k(mxdpts,mxdsc)             !<  log derivatives non-local KB pseudopotential

! local variables

  character(len=15)                ::  fname

  real(REAL64)                     ::  x_min, x_max, y_min, y_max

  integer                          ::  iw

  character(len=12), allocatable   ::  label(:)
  character(len=12), allocatable   ::  label2(:)

! counters

  integer                  :: i


  allocate(label(nsc), label2(nsc))

  call atom_plot_ln_write(ioplot, 'log_true.gp', c_val, nsc, numtotal,   &
      ehist, dlog_t, mxdsc, mxdpts)


  call atom_plot_ln_write(ioplot, 'log_ppot.gp', c_val, nsc, numtotal,   &
      ehist, dlog_p, mxdsc, mxdpts)

  if(lkb) then

    call atom_plot_ln_write(ioplot, 'log_ppkb.gp', c_val, nsc, numtotal, &
        ehist, dlog_k, mxdsc, mxdpts)
    iw = 1
  else
    iw = 0
  endif

  do i = 1,nsc
    if(iso(i) == 0) then
      write(label(i),'("l=",i2)') lo(i)
    else
      write(label(i),'("l=",i2,", s=",f4.1)') lo(i), 0.5*iso(i)
    endif
  enddo

  do i=1,nsc

    write(iowrite,*)
    write(iowrite,*)

    fname='command_'//c_val(i)//'.gp'

    open (unit=ioplot,file=fname,status='unknown',form='formatted')

    x_min = evlow
    x_max = evhi
    y_min =-5.
    y_max = 5.


    call  atom_plot_command(ioplot, c_val(i), nsc, rpoint, 0, .TRUE.,    &
          x_min, x_max, y_min, y_max, iw, label, i, label2)

    call  atom_plot_spec(ioplot, c_val(i), nsc, rpoint, 0,               &
          x_min, x_max, y_min, y_max, iw, label, i, label2)

    call atom_plot_gnuplot(iotmp,c_val(i),'run  ')

  enddo

! cleanup

  call atom_plot_gnuplot(iotmp, c_val(1), 'clean')

!  call atom_plot_gnuplot(iotmp, c_val(1), 'files')

  deallocate(label, label2)

  return

end subroutine atom_plot_ln_outpdf
