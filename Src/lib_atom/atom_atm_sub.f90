!>  The recent 6.0X  version of the venerable Froyen/Berkeley/Martins/Troullier/... atomic program
!>
!>  \author       Sverre Froyen, Jose Luis Martins
!>  \version      6.0.5
!>  \date         22 June 2021
!>  \copyright    GNU Public License v2

subroutine atom_atm_sub(iowrite, ioread, filein, ioae, fileae, iopsd, filepsd)

! translated to f90 from atm.f version 5.805
! printing, 21 October 2021. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  ioread                           !<  default tape for reading
  character(len=*), intent(in)      ::  filein                           !<  name of default tape for reading

  integer, intent(in)               ::  ioae                             !<  default tape for writing
  character(len=*), intent(in)      ::  fileae                           !<  name of default tape for all-electron results

  integer, intent(in)               ::  iopsd                            !<  default tape for pseudopotential in parsec format
  character(len=*), intent(in)      ::  filepsd                          !<  name of default tape for reading pseudopotential in parsec format

! array dimensions

  integer                           ::  mxdorb                           !  dimension of the number of orbitals
  integer                           ::  mxdnr                            !  dimension of the number of radial points
  integer                           ::  mxdl                             !  dimension for angular momentum l

! local variables

  character(len=2)                  ::  nameat                           !  chemical symbol of the atom

  character(len=2)                  ::  ctype                            !  type of calculation flag (converted to itype)
  character(len=10)                 ::  ititle(5)                        !  title of calculation

  character(len=3)                  ::  kerker                           !  type of pseudopotential flag

  character(len=2)                  ::  icorr                            !  correlation type
  character(len=1)                  ::  ispp                             !  spin polarization  ' ', 's', 'r'

  real(REAL64)                      ::  znuc                             !  nuclear charge

  real(REAL64)                      ::  rmax                             !  maximum mesh radius
  real(REAL64)                      ::  aa                               !  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64)                      ::  bb                               !  a = exp(-aa)/znuc, b = 1/bb

  integer                           ::  ncore                            !  number of orbitals treated as core
  integer                           ::  nval                             !  number of valence orbitals

  integer, allocatable              ::  ni(:)                            !  valence principal quantum number
  integer, allocatable              ::  li(:)                            !  valence angular quantum number
  real(REAL64), allocatable         ::  zd(:), zu(:)                     !  occupation of valence down and up orbitals
  real(REAL64), allocatable         ::  evd(:)                           !  default eigenvalue

  real(REAL64)                      ::  etotal                           !  total energy

! finds dimensions

  open(unit=ioread, file=trim(filein), status='OLD', form='FORMATTED')

  call atom_atm_input_size(ioread, iopsd, filepsd, mxdnr, mxdorb, mxdl)

  close(unit=ioread)

! gets configuration

  open(unit=ioread, file=trim(filein), status='OLD', form='FORMATTED')

  allocate(ni(mxdorb))
  allocate(li(mxdorb))
  allocate(zu(mxdorb),zd(mxdorb))
  allocate(evd(mxdorb))

  call atom_atm_read_input(nameat, ctype, ititle, kerker, icorr, ispp,   &
      znuc, rmax, aa, bb, ncore, nval, ni, li, zd, zu, evd,              &
      ioread, mxdorb)

  close(unit=ioread)

! calculates the ae structure

  write(6,*)
  write(6,*) '  Starting all-electron calculation'
  write(6,*)

  call atom_atm_scf(etotal,                                              &
      nameat, ctype, ititle, kerker, icorr, ispp,                        &
      znuc, rmax, aa, bb, ncore, nval, ni, li, zd, zu, evd,              &
      iowrite, ioae, fileae, iopsd, filepsd,                             &
      mxdnr, mxdorb, mxdl)

  deallocate(ni)
  deallocate(li)
  deallocate(zu,zd)
  deallocate(evd)

  return

end subroutine atom_atm_sub
