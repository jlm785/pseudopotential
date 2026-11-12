!>  Swows the comparison of log derivatives of all-electron and pseudo potentials.
!>
!>  \author       Jose Luis Martins
!>  \version      6.2.0
!>  \date         September 2021. 25 May 2022. 12 November 2026.
!>  \copyright    GNU Public License v2

program plot_ln_show


! Initializes lkb. 25 May 2022. JLM
! tapes >= 10 consistent with the other programs. 7 November 2026. claude+JLM
! rpoint from the core radii used (atom.dat) as in atom_all. 12 November 2026. claude

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

  integer                           ::  ioread                           !  tape number for atom.dat
  character(len=15)                 ::  filein                           !  file name for the input of atom and psd_gen

  integer                           ::  ioplot                           !  tape number for plots
  integer                           ::  iotmp                            !  tape number for temporary command files

  logical, external                 ::  atom_plot_gnuplot_ok             !  checks if gnuplot is installed

! other variables

  logical                           ::  lkb                              !  If KB results are included in the plots
  real(REAL64)                      ::  rpoint                           !  radius where log derivatives are calculated

  logical                           ::  lint                             !  interactive run

  real(REAL64)                      ::  rc(0:lc)                         !  core radii used to generate the pseudopotential
  real(REAL64)                      ::  rc_tab(0:lc)                     !  tabulated core radii
  character(len=30)                 ::  status                           !  quality of psedopotential
  real(REAL64)                      ::  cfac, rcfac, cfactau             !  partial core parameters (not used)
  logical                           ::  lex, lres

  character(len=2)                  ::  nameat

! parameters

  real(REAL64), parameter    ::  ONE = 1.0_REAL64

! counters

  integer     ::  l


  lkb = .TRUE.
  lint = .TRUE.
!  lint = .FALSE.

  iowrite = 6

  iodata = 12
  datafile = 'datafile.dat   '

  iopsd = 13
  pseudo = 'pseudo.dat     '
  iokb = 16
  pseudokb = 'pseudokb.dat   '

  ioread = 11
  filein = 'atom.dat       '

  ioplot = 30
  iotmp = 23

  if(.not. atom_plot_gnuplot_ok(iotmp)) then
    write(6,*) '  plot_ln_show needs gnuplot.  Stopping.'
    write(6,*)

    stop

  endif

  call atom_psd_print_info(iowrite, iodata, datafile, nameat)

  if(lint) then

    call atom_plot_ln_setup(datafile, pseudo, pseudokb, lkb, rpoint)

  else

!   core radii used by psd_gen (or atom_all): those of atom.dat,
!   or from the table if atom.dat does not exist or has no core radii
!   only used if the hard coded lint is changed!

    lres = .FALSE.
    inquire(file = trim(filein), exist = lex)
    if(lex) then
      call atom_psd_read_atom(ioread, filein, nameat, rc,                &
          cfac, rcfac, cfactau, lres)
    endif

    if(.not. lres) then
      call atom_p_tbl_psd_tm2(nameat, rc_tab, status)
      do l = 0,lc
        rc(l) = rc_tab(l)
      enddo
    endif

    rpoint = ONE
    do l = 0,lc
      if(1.2*rc(l) > rpoint) rpoint = 1.2*rc(l)
    enddo

  endif

  call atom_plot_ln_sub(iowrite, iodata, datafile, iopsd, pseudo,        &
         iokb, pseudokb, ioplot, iotmp, lkb, rpoint, lint)

  stop

end program plot_ln_show
