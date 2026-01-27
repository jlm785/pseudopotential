!>  Manages the pseudopotential generation toolchain
!>
!>  \author       Jose Luis Martins
!>  \version      6.1.0
!>  \date         22 June 2021, 26 January 2026.
!>  \copyright    GNU Public License v2

program psd_gen

! names files, bit rot of files. 26 January 2026. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter                ::  lc = 4                           !  maximum valence angular momentum

! tape identification

  integer                           ::  iowrite                          !  default output tape

  integer                           ::  ioread                           !  default tape for reading
  character(len=8)                  ::  filein                           !  name of default tape for reading

  integer                           ::  ioae                             !  default tape for writing
  character(len=12)                 ::  fileae                           !  name of default tape for all-electron results

  integer                           ::  iopsd                            !  default tape for pseudopotential in old format
  character(len=10)                 ::  filepsd                          !  name of default tape for writing pseudopotential in old format

  integer                           ::  ioreal                           !  default tape for pseudopotential in real space (parsec) format
  character(len=7)                  ::  filereal                         !  name of default tape for reading pseudopotential in real space (parsec) format
  character(len=4)                  ::  sfilesiesta                      !<  suffix of default tape for writing pseudopotential in siesta format
  character(len=10)                 ::  sfileparsec                      !<  suffix of default tape for writing pseudopotential in parsec format

  integer                           ::  ioplot                           !  default tape for plot file
  character(len=8)                  ::  fileplot                         !  name of default tape for plot file

! core radii

  real(REAL64)                      ::  rc(0:lc)                         !  core radius r_c(l)
  real(REAL64)                      ::  cfac                             !  criteria for pseudo-core charge
  real(REAL64)                      ::  rcfac                            !  pseudo-core radius

  real(REAL64)                      ::  rc_tab(0:lc)                     !  tabulated core radii
  real(REAL64)                      ::  rz(0:lc)                         !  zero radius l
  real(REAL64)                      ::  rx(0:lc)                         !  extrema radius l
  integer                           ::  lpmax                            !  maximum value of l in pseudo

! other variables

  integer                           ::  ifcore                           !  0: no partial core correction; 1 : xc partial core correction; 2 : not recommended
  integer               ::  ifcore_atm, ifcore_dat, ifcore_tab

  logical                           ::  lint                             !  interactive run


  lint = .TRUE.
!  lint = .FALSE.

  iowrite = 6

  ioread = 3
  filein = 'atom.dat'

  iopsd = 2
  filepsd = 'pseudo.dat'

  ioae = 7
  fileae = 'datafile.dat'

  ioreal = 4
  filereal = 'psd.pot'
  sfileparsec = '_POTRE.DAT'
  sfilesiesta = '.psf'

  ioplot = 8
  fileplot = 'plot.dat'


! gets default core radii, wave-function zero and extrema

  call atom_psd_core_radii(rc, cfac, rcfac, rc_tab, rz, rx, lpmax,       &
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

  endif

! constructs the pseudopotential

  call atom_psd_sub(rc, ifcore, cfac, rcfac, .FALSE.,                    &
         iowrite, ioae, fileae,                                          &
         iopsd, filepsd, ioreal, filereal, sfileparsec, sfilesiesta,     &
         ioplot, fileplot)


  stop

end program psd_gen
