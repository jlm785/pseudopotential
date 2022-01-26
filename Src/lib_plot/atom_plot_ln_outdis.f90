!>  Generates the log derivative plots and uses gnuplot to show it on the screen
!>
!>  \author       Norm Troullier, Manuel Maria Alemany, J.L.Martins
!>  \version      6.0.3
!>  \date         1990s, 23 August 2021. 25 September 2021.
!>  \copyright    GNU Public License v2

subroutine atom_plot_ln_outdis(ioplot, iotmp, lkb,                       &
       numtotal, evlow, evhi, rpoint,                                    &
       c_val, nsc, lo, iso, ehist, dlog_t, dlog_p, dlog_k,               &
       iowrite, mxdsc, mxdpts)


! adapted from the old program August 2021. JLM
! mxdsc. 25 September 2021. JLM

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

  logical                          ::  lnew

  integer                          ::  iw
  integer                          ::  lcentr

  character(len=12), allocatable   ::  label(:)
  character(len=12), allocatable   ::  label2(:)


! counters

  integer                  :: i, k


  allocate(label(nsc), label2(nsc))

  write(iowrite,*)
  write(iowrite,*) '... WRITING TRUE DERIVATIVES IN log_true.gp ...'
  write(iowrite,*)

  open (unit=ioplot,file='log_true.gp',status='unknown',form='formatted')

  write(ioplot,'("#    E(Ry)",10x,20(a4,11x))') (c_val(i),i=1,nsc)
  do i = 1,numtotal
    write(ioplot,'(2x,e13.6,20(2x,e13.6))') ehist(i),(dlog_t(i,k),k=1,nsc)
  enddo
  close(ioplot)

  write(iowrite,*)
  write(iowrite,*) '... WRITING SEMI-LOCAL PSEUDO DERIVATIVES IN log_ppot.gp ...'
  write(iowrite,*)

  open (unit=ioplot,file='log_ppot.gp',status='unknown',form='formatted')

  write(ioplot,'("#    E(Ry)",10x,20(a4,11x))') (c_val(i),i=1,nsc)
  do i = 1,numtotal
    write(ioplot,'(2x,e13.6,20(2x,e13.6))') ehist(i),(dlog_p(i,k),k=1,nsc)
  enddo
  close(ioplot)

  if(lkb) then
    write(iowrite,*)
    write(iowrite,*) '... WRITING NON-LOCAL KB PSEUDO DERIVATIVES IN log_ppkb.gp ...'
    write(iowrite,*)

    open (unit=ioplot,file='log_ppkb.gp',status='unknown',form='formatted')

    write(ioplot,'("#    E(Ry)",10x,20(a4,11x))') (c_val(i),i=1,nsc)
    do i = 1,numtotal
      write(ioplot,'(2x,e13.6,20(2x,e13.6))') ehist(i),(dlog_k(i,k),k=1,nsc)
    enddo
    close(ioplot)
    iw = 1
  else
    iw = 0
  endif

  write(iowrite,*)
  write(iowrite,*)
  write(iowrite,*) '****************************************************'
  write(iowrite,*) '*** IF YOU PROCEED OTHER GNUPLOTS WILL BE KILLED ***'
  write(iowrite,*) '****************************************************'

  do i = 1,nsc
    if(iso(i) == 0) then
      write(label(i),'("8x,l=",i2)') lo(i)
    else
      write(label(i),'("l=",i2,", s=",f4.1)') lo(i), 0.5*iso(i)
    endif
  enddo

  do i = 1,nsc

    write(iowrite,*)
    write(iowrite,*)

    fname='command_'//c_val(i)//'.gp'

    open (unit=ioplot,file=fname,status='unknown',form='formatted')

    x_min = evlow
    x_max = evhi
    y_min =-20.
    y_max = 20.

    call  atom_plot_command(ioplot, c_val(i), nsc, rpoint, 0, .FALSE., &
        x_min, x_max, y_min, y_max, iw, label, i, label2)

    call atom_plot_gnuplot(iotmp,c_val(i),'run  ')

    call atom_plot_questions('range', x_min, x_max, y_min, y_max, lcentr, lnew)

    if(lnew) then

      call atom_plot_gnuplot(iotmp,c_val(i),'kill ')

      call  atom_plot_command(ioplot, c_val(i), nsc, rpoint, 0, .FALSE., &
          x_min, x_max, y_min, y_max, iw, label, i, label2)

      call atom_plot_gnuplot(iotmp,c_val(i),'run  ')

    endif

    call atom_plot_questions('pdf  ', x_min, x_max, y_min, y_max, lcentr, lnew)

    call atom_plot_gnuplot(iotmp,c_val(i),'kill ')

    if(lnew) then

      call  atom_plot_command(ioplot, c_val(i), nsc, rpoint, 0, .TRUE., &
          x_min, x_max, y_min, y_max, iw, label, i, label2)

      call atom_plot_gnuplot(iotmp,c_val(i),'run  ')

    endif

  enddo

! cleanup

  call atom_plot_gnuplot(iotmp, c_val(1), 'clean')

  call atom_plot_questions('clean',x_min, x_max, y_min, y_max, lcentr, lnew)

  if(lnew) then
    call atom_plot_gnuplot(iotmp, c_val(1), 'files')
  endif

  deallocate(label, label2)

  return

end subroutine atom_plot_ln_outdis
