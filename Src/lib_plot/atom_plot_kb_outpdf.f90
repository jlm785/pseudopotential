!>  Writes the files with data for subsequent plotting
!>  of wave-functions and potentials with gnuplot
!>
!>  \author       Peter Schuster, Manuel Maria Alemany, Jose Luis Martins
!>  \version      6.0.5
!>  \date         4 September 2021, 9 October 2021.
!>  \copyright    GNU Public License v2

subroutine atom_plot_kb_outpdf(iocomm, iotmp, numw, numb, nump, nftp,    &
     iw, labelw, lo_w,                                                   &
     ib, labelb, lo_b,                                                   &
     ip, labelp,                                                         &
     iq, labelq,                                                         &
     labelek,                                                            &
     labelv,                                                             &
     yminwp, ymaxwp, yminp, ymaxp, yminq, ymaxq,                         &
     xmaxek, yminek, ymaxek, yminv, ymaxv)

! adapted from all-electron+pseudopotential plot code.
! kinetic, 9 October 2021. JLM
! nump can be zero. 18 October 2021. JLM
! printing, 21 October 2021. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  iocomm                           !<  tape number for gnuplot command file
  integer, intent(in)               ::  iotmp                            !<  tape number for temporary files

  integer, intent(in)               ::  numw                             !<  dimension of number of wave-functions
  integer, intent(in)               ::  numb                             !<  dimension of number of basis-functions
  integer, intent(in)               ::  nump                             !<  dimension of projectors
  integer, intent(in)               ::  nftp                             !<  dimension of number of Fourier transform points

  integer, intent(in)               ::  iw                               !<  number of wave-functions
  character(len=12), intent(in)     ::  labelw(numw)                     !<  label of wave-functions
  integer, intent(in)               ::  lo_w(numw)                       !<  angular momentum of wave-function

  integer, intent(in)               ::  ib                               !<  number of basis-functions
  character(len=12), intent(in)     ::  labelb(numb)                     !<  label of basis-functions
  integer, intent(in)               ::  lo_b(numb)                       !<  angular momentum of wave-function

  integer, intent(in)               ::  ip                               !<  number of projectors
  character(len=12), intent(in)     ::  labelp(nump)                     !<  label of projectors

  integer, intent(in)               ::  iq                               !<  number of projectors
  character(len=12), intent(in)     ::  labelq(nump)                     !<  label of projectors

  character(len=12), intent(in)     ::  labelek                          !<  label of kinetic energy

  character(len=12), intent(in)     ::  labelv                           !<  label of local potential

  real(REAL64), intent(in)          ::  yminwp, ymaxwp                   !<  min and max
  real(REAL64), intent(in)          ::  yminp, ymaxp                     !<  min and max
  real(REAL64), intent(in)          ::  yminq, ymaxq                     !<  min and max
  real(REAL64), intent(in)          ::  xmaxek, yminek, ymaxek           !<  min and max
  real(REAL64), intent(in)          ::  yminv, ymaxv                     !<  min and max

! local variables

  character(len=12), allocatable    ::  label_bl(:), label_wl(:)

  integer, allocatable              ::  ipoint_b(:)
  integer                           ::  ipoint_w

  real(REAL64)          ::  x_min, x_max                                 !  x range
  real(REAL64)          ::  y_min, y_max                                 !  y range

  integer               ::  lmax
  integer               ::  nbas, nw, nlabel

  character(len=4)      ::  sfname
  real(REAL64)          ::  zion                             !<  ionic charge

! parameters

  real(REAL64), parameter   ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  real(REAL64), parameter   ::  HALF = 0.5_REAL64

! counters

  integer     ::  i, l


  x_min = ZERO
  x_max = 10*ONE

  y_min = HALF*floor(2*yminwp)
  y_max = HALF*(floor(2*ymaxwp)+1)

  lmax = 0
  do i = 1,numw
    if(lmax < lo_w(i)) lmax = lo_w(i)
  enddo
  do i = 1,numb
    if(lmax < lo_b(i)) lmax = lo_b(i)
  enddo

! loop over angular momentum of basis functions

  do l = 0,lmax

    write(sfname,'("bas",i1)') l

    ipoint_w = 0
    do i = 1,iw
      if(lo_w(i) == l) ipoint_w = i
    enddo

    nbas = 0
    do i = 1,ib
      if(lo_b(i) == l) nbas = nbas + 1
    enddo
    allocate(ipoint_b(nbas))
    nbas = 0
    do i = 1,ib
      if(lo_b(i) == l) then
        nbas = nbas + 1
        ipoint_b(nbas) = i
      endif
    enddo

    if(ipoint_w == 0) then
      nw = 0
    else
      nw = 1
    endif

    nlabel = max(nbas,nw)

    if(nlabel > 0) then

      allocate(label_wl(nlabel))
      allocate(label_bl(nlabel))
      do i = 1,nlabel
        label_wl(i) = '            '
        label_bl(i) = '            '
      enddo

      if(ipoint_w /= 0) then
        nw = 1
        label_wl(1) = labelw(ipoint_w)
      else
        nw = 0
      endif

      if(nbas > 0) then
        do i = 1,nbas
          label_bl(i) = labelb(ipoint_b(i))
        enddo
      endif

      call atom_plot_command(iocomm, sfname, nlabel, zion, 0, .TRUE.,    &
         x_min, x_max, y_min, y_max, nw, label_wl, nbas, label_bl)

      call atom_plot_gnuplot(iotmp, sfname, 'run  ')

      deallocate(label_wl)
      deallocate(label_bl)

    endif

    deallocate(ipoint_b)

  enddo


! Local potetntial

  x_min = ZERO
  x_max = 5*ONE

  y_min = yminv
  y_max = ymaxv

  sfname = 'kbvk'

  call atom_plot_command(iocomm, sfname, 1, zion, 0, .TRUE.,             &
         x_min, x_max, y_min, y_max, 1, labelv, 0, labelb)

  call atom_plot_gnuplot(iotmp, sfname, 'run  ')

! projectors real space

  if(nump > 0) then

    x_min = ZERO
    x_max = 5*ONE

    y_min = HALF*floor(2*yminp)
    y_max = HALF*(floor(2*ymaxp)+1)

    sfname = 'kbop'

    call atom_plot_command(iocomm, sfname, nump, zion, 0, .TRUE.,          &
           x_min, x_max, y_min, y_max, ip, labelp, 0, labelb)

    call atom_plot_gnuplot(iotmp, sfname, 'run  ')

  endif

! projetors reciprocal space

  if(nump > 0) then

    x_min = ZERO
    x_max = 12*ONE

    y_min = HALF*floor(2*yminq)
    y_max = HALF*(floor(2*ymaxq)+1)

    sfname = 'kbtr'

    call atom_plot_command(iocomm, sfname, nftp, zion, 0, .TRUE.,          &
           x_min, x_max, y_min, y_max, iq, labelq, 0, labelb)

    call atom_plot_gnuplot(iotmp, sfname, 'run  ')

  endif

! integral kinetic energy

  x_min = ZERO
  x_max = xmaxek

  y_min = yminek
  y_max = ymaxek

  sfname = 'kbek'

  call atom_plot_command(iocomm, sfname, 1, zion, 0, .TRUE.,             &
         x_min, x_max, y_min, y_max, 1, labelek, 0, labelb)

  call atom_plot_gnuplot(iotmp, sfname, 'run  ')


! cleanup

  call atom_plot_gnuplot(iotmp, sfname, 'clean')


  return

end subroutine atom_plot_kb_outpdf
