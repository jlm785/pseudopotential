!>  Calls the whole atomic pseudo-potential test plot chain
!>
!>  \author       Jose Luis Martins
!>  \version      6.2.0
!>  \date         21 October 2021. 26 January 2026. 12 November 2026.
!>  \copyright    GNU Public License v2

program atom_all


! Initializes lkb. 25 May 2022. JLM
! Several real space pseudopotential formats, 26 January 2026. JLM
! up to quadruple zeta basis, 24 October 2026. JLM+claude
! tbasis, lmax_bas allocated with mxdset, 24 October 2026. JLM+claude
! tapes >= 10 consistent with the other programs. 7 November 2026. claude+JLM
! tested pseudopotentials, arguments identified by content,
! rpoint from the core radii before rounding to the grid. 12 November 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter                ::  lc = 4                           !  maximum valence angular momentum

! tape stuff

  integer                           ::  iowrite                          !  default output tape

  integer                           ::  ioread                           !  default tape for reading configuration file
  character(len=8)                  ::  filein                           !  name of default tape for reading configuration file

  integer                           ::  ioae                             !  default tape for writing all-electron results
  character(len=12)                 ::  fileae                           !  name of default tape for all-electron results

  integer                           ::  iopsd                            !  default tape for pseudopotential in old format
  character(len=10)                 ::  filepsd                          !  name of default tape for reading pseudopotential in old format

  integer                           ::  ioreal                           !  default tape for pseudopotential in real space (parsec) format
  character(len=7)                  ::  filereal                         !  name of default tape for reading/writing pseudopotential in real space (parsec) format
  character(len=10)                 ::  sfileparsec                      !  suffix of default tape for writing pseudopotential in parsec format
  character(len=4)                  ::  sfilesiesta                      !  suffix of default tape for writing pseudopotential in siesta format


  integer                           ::  ioplot                           !  default tape for plot file
  character(len=8)                  ::  fileplot                         !  name of default tape for plot file

  integer                           ::  iokb                             !  default tape for pseudopotential in KB format
  character(len=12)                 ::  sfilekb                          !  suffix for default tape for writing pseudopotential in KB format
  character(len=5)                  ::  sfilepsp8                        !  suffix for default tape for writing pseudopotential in psp8 abinit format

  integer                           ::  ioupf                            !  default tape for pseudopotential in UPF format
  character(len=6)                  ::  sfileupf                         !  suffix for default tape for writing pseudopotential in UPF format

  integer                           ::  iopsdkb                          !  default tape for KB pseudopotential in real space
  character(len=15)                 ::  filepsdkb                        !  name of default tape for writeing KB pseudopotential in real space

  integer                           ::  ioplotkb                         !  tape for later plotting
  character(len=15)                 ::  fileplotkb                       !  name of file with plot information

  integer                           ::  iost                             !  data files use tapes iost+1 to iost+4
  integer                           ::  iocomm                           !  tape number for gnuplot command file
  integer                           ::  iotmp                            !  tape number for temporary files

! core radii

  real(REAL64)                      ::  rc(0:lc)                         !  core radius r_c(l)
  real(REAL64)                      ::  cfac                             !  criteria for pseudo-core charge
  real(REAL64)                      ::  rcfac                            !  pseudo-core radius
  real(REAL64)                      ::  cfactau                          !  core correction factor for kinetic energy density

  real(REAL64)                      ::  rc_tab(0:lc)                     !  tabulated core radii
  real(REAL64)                      ::  rz(0:lc)                         !  zero radius l
  real(REAL64)                      ::  rx(0:lc)                         !  extrema radius l
  integer                           ::  lpmax                            !  maximum value of l in pseudo

! KB variables

  integer                           ::  llocal                           !  angular momentum for local potential (negative: maximum of l-dependent)
  integer                           ::  nql                              !  number of points for the Fourier grid of local potential and densities
  real(REAL64)                      ::  delql                            !  spacing of points in Fourier grid

  integer                           ::  lmax_pot                         !  maximum angular momentum in potential

  integer                           ::  mxdl                             !  dimension for angular momentum
  integer                           ::  mxdnr                            !  dimension for radial grid
  integer                           ::  mxdset                           !  dimension for number of atomic basis sets

  integer                           ::  n_bsets                          !  number of atomic basis sets

  character(len=3), allocatable     ::  tbasis(:)                        !  type of basis
  integer, allocatable              ::  lmax_bas(:)                      !  maximum angular momentum in basis

  integer, allocatable              ::  n_bas(:,:)                       !  basis functions for angular momentum l
  real(REAL64), allocatable         ::  r_bas(:,:,:)                     !  cutoff for the basis (up to quadruple zeta), 0: energy derivative of previous
  integer, allocatable              ::  nz_bas(:,:,:)                    !  number of non-trivial zeroes in basis function
  real(REAL64), allocatable         ::  r_siesta(:)                      !  cutoff using the SIESTA recipe
  real(REAL64), allocatable         ::  r_99(:)                          !  radius with 99% of charge


! other variables

  logical                           ::  lex
  character(len=2)                  ::  nameat                           !  chemical symbol of the atom

  integer                           ::  narg

  integer                           ::  ifcore                           !  0: no partial core correction; 1 : xc partial core correction; 2 : not recommended
  integer                           ::  ifcore_atm, ifcore_dat, ifcore_tab

  logical                           ::  lkb                              !  If KB results are included in the plots
  real(REAL64)                      ::  rpoint                           !  radius where log derivatives are calculated

  logical                           ::  lint                             !  interactive run

  integer                           ::  iz                               !  atomic number

  character(len=2)                  ::  firstarg                         !  first argument
  character(len=3)                  ::  tbasisin                         !  type of basis in program line
  character(len=3)                  ::  pccin                            !  type of partial core correction in program line
  character(len=4)                  ::  vtypein                          !  label of tested pseudopotential, '    ' default, 'NONE' old tables
  character(len=4)                  ::  argin                            !  argument in program line
  integer                           ::  larg                             !  length of argument
  integer                           ::  iarg                             !  counter

! parameters

  real(REAL64), parameter    ::  ONE = 1.0_REAL64

! counters

  integer     ::  l


! interactivity
!  lint = .TRUE.

  lint = .FALSE.

! only warnings go to default output
! iowrite = 6
  iowrite = 10

  iopsd = 13
  filepsd = 'pseudo.dat'

  ioread = 11
  filein = 'atom.dat'

  ioae = 12
  fileae = 'datafile.dat'

  ioreal = 14
  filereal = 'psd.pot'
  sfileparsec = '_POTRE.DAT'
  sfilesiesta = '.psf'

  ioplot = 15
  fileplot = 'plot.dat'

  iokb = 18
  sfilekb = '_POTKB_F.DAT'

  ioupf = 19
  sfileupf = 'TM.UPF'
  sfilepsp8 = '.psp8'

  iopsdkb = 16
  filepsdkb = 'pseudokb.dat   '

  ioplotkb = 17
  fileplotkb = 'plot_kb.dat    '

  iost = 30
  iocomm = 21
  iotmp = 23


  if(iowrite /= 6) then
    open(unit=iowrite,file='atom.out',form='FORMATTED',status='UNKNOWN')
  endif

! checks if the program is run with line arguments

  narg = command_argument_count()

  tbasisin = '   '
  pccin = '   '

  if(narg > 0) then

    call get_command_argument(1,firstarg)

    if(firstarg == '-h') then

      write(6,*)
      write(6,*) '  Run the code without arguments or'
      write(6,*) '  with a chemical symbol as first argument, or'
      write(6,*) '  with -x  as first argument to enter an interactive mode.'
      write(6,*)
      write(6,*) '  The other arguments, in any order, can be:'
      write(6,*)
      write(6,*) '  yes/no for with/without partial core corrections'
      write(6,*) '  for exchange and correlation.'
      write(6,*)
      write(6,*) '  SZ/DZ/DZP/... for choice of basis set.'
      write(6,*)
      write(6,*) '  the label of a tested pseudopotential (HARD, SOFT,...),'
      write(6,*) '  or NONE for the old tables (educated guess).'
      write(6,*) '  The default is the first tested pseudopotential.'
      write(6,*)
      write(6,*) '  If the arguments are not present default values are used.'
      write(6,*)

      stop

    endif

    if(firstarg == '-x') then

      lint = .TRUE.

      write(6,*)
      write(6,*)
      write(6,*) '  Enter chemical symbol'
      write(6,*)

      nameat = '  '
      read(5,*) nameat

      call atom_p_tbl_charge(nameat,iz)

      if(iz > 150) then

        write(6,*)
        write(6,*) '  Unknown element, STOPPING'
        write(6,*)

        stop

      endif

      call atom_atm_tested_choose(nameat, vtypein)

      call atom_atm_write_sub(ioread, filein, nameat, vtypein)
      write(6,*)

    else

      nameat = firstarg

      call atom_p_tbl_charge(nameat,iz)

      if(iz > 150) then

        write(6,*)
        write(6,*) '  Unknown element, STOPPING'
        write(6,*)

        stop

      endif

      inquire(file=filein, exist = lex)

      if(lex) then
        write(6,*)
        write(6,*) '  Cannot overwrite ',filein
        write(6,*) '  rm ',filein
        write(6,*) '  before re-running the code'
        write(6,*)

        stop

      endif

!     more then one argument, identified by content

      vtypein = '    '

      do iarg = 2,narg

        call get_command_argument(iarg, argin, larg)

        if(larg > 4) then
          write(6,*)
          write(6,'("  Argument ",i2," is too long (more than 4 ",       &
               &   "characters), STOPPING")') iarg
          write(6,*)

          stop

        endif

        call chrcap(argin, 4)

        if(argin == 'YES ') then
          pccin = 'yes'
        elseif(argin == 'NO  ') then
          pccin = 'no '
        elseif(argin(2:2) == 'Z') then
          tbasisin = argin(1:3)
        else
          vtypein = argin
        endif

      enddo

      if(iz < 3) pccin = '   '

      call atom_atm_write_sub(ioread, filein, nameat, vtypein)
      write(6,*)

    endif

  else

!   no line arguments when the code was run

    inquire(file=filein, exist = lex)

    if(lex) then
      write(6,*)
      write(6,*) '  using file ',filein
      write(6,*)
    else

      write(6,*)
      write(6,*)
      write(6,*) '  Enter chemical symbol'
      write(6,*)

      nameat = '  '
      read(5,*) nameat

      call atom_p_tbl_charge(nameat,iz)

      if(iz > 150) then

        write(6,*)
        write(6,*) '  Unknown element, STOPPING'
        write(6,*)

        stop

      endif

      call atom_atm_tested_choose(nameat, vtypein)

      call atom_atm_write_sub(ioread, filein, nameat, vtypein)
      write(6,*)

    endif

  endif

! runs the all-ellecron calculation

  call atom_atm_sub(iowrite, ioread, filein, ioae, fileae, iopsd, filepsd)

! gets default core radii, wave-function zero and extrema

  call atom_psd_core_radii(rc, cfac, rcfac, cfactau, rc_tab,             &
         rz, rx, lpmax,                                                  &
         iowrite, ioae, ioread, filein, fileae)

  if(lint) then

    call atom_psd_change_radii(lpmax, rc, rc_tab, rz, rx)

  endif

  call atom_psd_ifcore(ioread, filein, ioae, fileae,                     &
           ifcore_atm, ifcore_dat, ifcore_tab)

  if(lint) then

    call atom_psd_change_ifcore(ifcore_atm, ifcore_dat, ifcore_tab,      &
           ifcore)

  else

    if(ifcore_atm >= 0) then
      ifcore = ifcore_atm
    else
      ifcore = ifcore_tab
    endif

    if(pccin == 'yes') ifcore = 1
    if(pccin == 'no ') ifcore = 0

  endif

  if(lint) then

    call atom_psd_change_cfac(ifcore, cfac, rcfac, cfactau)

  endif

! default radius for the logarithmic derivatives, from the core radii
! before atom_psd_sub rounds them to the grid

  rpoint = ONE
  do l = 0,lc
    if(1.2*rc(l) > rpoint) rpoint = 1.2*rc(l)
  enddo

! constructs the pseudopotential

  call atom_psd_sub(rc, ifcore, cfac, rcfac, cfactau, .FALSE.,           &
         iowrite, ioae, fileae,                                          &
         iopsd, filepsd, ioreal, filereal, sfileparsec, sfilesiesta,     &
         ioplot, fileplot)

! plots all-electron pseudo comparison

  call atom_plot_sub(iowrite, ioplot, fileplot, iost, iocomm, iotmp, lint)


! KB transformation

  llocal = -1
  nql = 4000
  delql = 0.01_REAL64

! default basis sets

  if(tbasisin(2:2) == 'Z') then
    n_bsets = 1
  else
    n_bsets = 2
  endif

  mxdset = n_bsets + 2

  allocate(tbasis(mxdset))
  allocate(lmax_bas(mxdset))

  tbasis(:) = '   '
  lmax_bas(:) = 0

  if(tbasisin(2:2) == 'Z') then
    tbasis(1) = tbasisin
  else
    tbasis(1) = 'DZP'
    tbasis(2) = 'SZ'
  endif


  call atom_kb_psd_in_real_size(ioreal, filereal, lmax_pot, mxdnr)

! give enough headroom for larger basis sets, increase if needed

  mxdl = lmax_pot + 2

  allocate(n_bas(0:mxdl,mxdset))
  allocate(r_bas(4,0:mxdl,mxdset))
  allocate(nz_bas(4,0:mxdl,mxdset))
  allocate(r_siesta(0:mxdl))
  allocate(r_99(0:mxdl))

  call atom_kb_basis_cutoff(n_bsets, tbasis, lmax_pot, nameat,           &
      lmax_bas, n_bas, r_bas, nz_bas, r_siesta, r_99,                    &
      iowrite, ioreal, filereal,                                         &
      mxdl, mxdset, mxdnr)

! llocal of the tested pseudopotential or from table

  call atom_kb_llocal_default(nameat, llocal, iowrite, ioreal, filereal)

  if(llocal > lmax_pot) llocal = lmax_pot

  if(lint) then

    call atom_kb_setup_pseudo(llocal, lmax_pot, nql, delql)

    call atom_kb_setup_basis(n_bsets, tbasis,                            &
        lmax_bas, n_bas, r_bas, nz_bas, r_siesta, r_99,                  &
        mxdl, mxdset)

  endif


  call atom_kb_sub(llocal, nql, delql, nql, delql,                       &
      n_bsets, lmax_bas, n_bas, r_bas, nz_bas,                           &
      iowrite, ioreal, filereal, iokb, sfilekb,                          &
      ioupf, sfileupf, sfilepsp8,                                        &
      iopsdkb, filepsdkb, ioplotkb, fileplotkb,                          &
      mxdnr, mxdl, mxdset)

! plots with KB projectors

  call atom_plot_kb_sub(iowrite, ioplotkb, fileplotkb, iost, iocomm, iotmp, lint)

! plots of the logarithmic derivatives

  call atom_psd_print_info(iowrite, ioae, fileae, nameat)

  lkb = .TRUE.

  if(lint) then

    write(6,*)' Enter distance (in a.u.) at which the'
    write(6,*)' logaritmic derivatives are calculated'
    read(5,*) rpoint
    write(6,*)
    write(6,*)

  endif

  call atom_plot_ln_sub(iowrite, ioae, fileae, iopsd, filepsd,           &
         iopsdkb, filepsdkb, iost, iotmp, lkb, rpoint, lint)


! runs the final tests

  call atom_atm_test(iowrite, ioread, filein,                            &
         iopsd, filepsd, iopsdkb, filepsdkb)

  deallocate(n_bas, r_bas, nz_bas)
  deallocate(r_siesta, r_99)
  deallocate(tbasis)
  deallocate(lmax_bas)

  stop

end program atom_all
