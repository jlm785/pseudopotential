!>  Swows the comparison of log derivatives of all-electron and pseudo potentials.
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.3
!>  \date         September 2021
!>  \copyright    GNU Public License v2

program plot_ln_show

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter                ::  lc = 4                           !  maximum valence angular momentum


  integer                           ::  iowrite                          !  default output

  integer                           ::  iodata                           !  tape number for datafile
  character(len=15)                 ::  datafile                         !  file name for all electron results

  integer                           ::  iopsd                            !  tape number for pseudo
  character(len=15)                 ::  pseudo                           !  file name for semi-local pseudopotential results

  integer                           ::  iokb                             !  tape number for pseudokb
  character(len=15)                 ::  pseudokb                         !  file name for non-local KB pseudopotential results

  integer                           ::  ioplot                           !  tape number for plots
  integer                           ::  iotmp                            !  tape number for temporary command files

! other variables

  logical                           ::  lkb                              !  If KB results are included in the plots
  real(REAL64)                      ::  rpoint                           !  radius where log derivatives are calculated

  logical                           ::  lint                             !  interactive run

  real(REAL64)                      ::  rc_tab(0:lc)                     !  tabulated core radii
  character(len=30)                 ::  status                           !  quality of psedopotential

  character(len=2)                  ::  nameat

! parameters

  real(REAL64), parameter    ::  ONE = 1.0_REAL64

! counters

  integer     ::  l


  lint = .TRUE.
!  lint = .FALSE.

  iowrite = 6

  iodata = 7
  datafile = 'datafile.dat   '

  iopsd = 8
  pseudo = 'pseudo.dat     '
  iokb = 9
  pseudokb = 'pseudokb.dat   '

  ioplot = 18
  iotmp = 19

  call atom_psd_print_info(iowrite, iodata, datafile, nameat)

  if(lint) then

    call atom_plot_ln_setup(datafile, pseudo, pseudokb, lkb, rpoint)

  else

    call atom_p_tbl_psd_tm2(nameat, rc_tab, status)

    rpoint = ONE
    do l = 0,lc
      if(1.2*rc_tab(l) > rpoint) rpoint = 1.2*rc_tab(l)
    enddo

  endif

  call atom_plot_ln_sub(iowrite, iodata, datafile, iopsd, pseudo,        &
         iokb, pseudokb, ioplot, iotmp, lkb, rpoint, lint)

  stop

end program plot_ln_show
