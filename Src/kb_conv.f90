program kb_conv

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


  integer                           ::  iowrite                          !  default output tape

  integer                           ::  ioparsec                         !  default tape for pseudopotential in parsec format
  character(len=7)                  ::  fileparsec                       !  name of default tape for reading pseudopotential in parsec format

  integer                           ::  iokb                             !  default tape for pseudopotential in KB format
  character(len=12)                  ::  sfilekb                         !  suffix for default tape for writing pseudopotential in KB format

  integer                           ::  ioupf                            !  default tape for pseudopotential in UPF format
  character(len=7)                  ::  sfileupf                         !  suffix for default tape for writing pseudopotential in UPF format

  integer                           ::  iopsdkb                          !  default tape for KB pseudopotential in real space
  character(len=15)                 ::  filepsdkb                        !  name of default tape for writeing KB pseudopotential in real space

  integer                           ::  ioplotkb                         !  tape for later plotting
  character(len=15)                 ::  fileplotkb                       !  name of file with plot information

  integer                           ::  llocal                           !  angular momentum for local potential (negative: maximum of l-dependent)
  integer                           ::  nql                              !  number of points for the Fourier grid of local potential and densities
  real(REAL64)                      ::  delql                            !  spacing of points in Fourier grid

  character(len=3)                  ::  tbasis                           !  type of basis

  logical                           ::  lint                             !  interactive run

  integer                           ::  lmax_bas                         !  maximum angular momentum in basis
  integer                           ::  lmax_pot                         !  maximum angular momentum in potential

  integer                           ::  mxdl                             !  dimension for angular momentum
  integer                           ::  mxdnr                            !  dimension for radial grid

  integer, allocatable              ::  n_bas(:)                         !  basis functions for angular momentum l
  real(REAL64), allocatable         ::  r_bas(:,:)                       !  cutoff for the basis (up to triple zeta and l = 4)
  integer, allocatable              ::  nz_bas(:,:)                      !  number of non-trivial zeroes in basis function
  real(REAL64), allocatable         ::  r_siesta(:)                      !  cutoff using the SIESTA recipe
  real(REAL64), allocatable         ::  r_99(:)                          !  radius with 99% of charge

  lint = .TRUE.
!  lint = .FALSE.

  iowrite = 6

  ioparsec = 7
  fileparsec = 'psd.pot'

  iokb = 9
  sfilekb = '_POTKB_F.DAT'

  ioupf = 9
  sfileupf = '_TM.UPF'

  iopsdkb = 10
  filepsdkb = 'pseudokb.dat   '

  ioplotkb = 11
  fileplotkb = 'plot_kb.dat    '

  llocal = -1
  nql = 4000
  delql = 0.01_REAL64
  tbasis = 'DZP'


  call atom_kb_psd_in_parsec_size(ioparsec, fileparsec, lmax_pot, mxdnr)

! give enough headroom for larger basis sets, increase if needed

  mxdl = lmax_pot + 2

  allocate(n_bas(0:mxdl), r_bas(3,0:mxdl), nz_bas(3,0:mxdl))
  allocate(r_siesta(0:mxdl), r_99(0:mxdl))

  call atom_kb_basis_cutoff(tbasis, llocal, lmax_pot,                    &
      lmax_bas, n_bas, r_bas, nz_bas, r_siesta, r_99,                    &
      iowrite, ioparsec, fileparsec, mxdl)

  if(lint) then
    call atom_kb_setup(tbasis, llocal, lmax_pot, nql, delql,             &
        lmax_bas, n_bas, r_bas, nz_bas, r_siesta, r_99,                  &
        mxdl)
  endif

  call atom_kb_sub(llocal, nql, delql, nql, delql,                       &
      lmax_bas, n_bas, r_bas, nz_bas,                                    &
      iowrite, ioparsec, fileparsec, iokb, sfilekb, ioupf, sfileupf,     &
      iopsdkb, filepsdkb, ioplotkb, fileplotkb, mxdnr, mxdl)

  deallocate(n_bas, r_bas, nz_bas)
  deallocate(r_siesta, r_99)

  stop

end program kb_conv


