!>  Tests the pseudopotential in the Kleinman-Bylander form
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.3
!>  \date         22 June 2021, 20 September 2021. JLM
!>  \copyright    GNU Public License v2

subroutine atom_kb_test_sub(iowrite, ioread, filein, iopsd, filepsd)



  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  ioread                           !<  default tape for reading
  character(len=*), intent(in)      ::  filein                           !<  name of default tape for reading

  integer, intent(in)               ::  iopsd                            !<  default tape for KB pseudopotential input
  character(len=*), intent(in)      ::  filepsd                          !<  name of default tape for reading KB pseudopotential

! array dimensions

  integer                           ::  mxdorb                           !  dimension of orbitals
  integer                           ::  mxdnr                            !  dimension of radial points
  integer                           ::  mxdl                             !  dimension of angular momentum

! local variables

  integer                           ::  lmax                             !  maximum angular value in potential
  integer                           ::  nr                               !  number of points
  integer                           ::  norb                             !  number of orbitals

  character(len=2)                  ::  nameat                           !  chemical symbol of the atom

  character(len=10)                 ::  ititle(5)                        !  title of calculation

  character(len=2)                  ::  icorr                            !  correlation type
  character(len=1)                  ::  ispp                             !  spin polarization  ' ', 's', 'r'

  integer                           ::  nval                             !  number of valence orbitals

  real(REAL64)                      ::  etotal                           !  total energy

  integer, allocatable              ::  ni(:)                            !  valence principal quantum number
  integer, allocatable              ::  li(:)                            !  valence angular quantum number
  real(REAL64), allocatable         ::  zd(:), zu(:)                     !  occupation of valence down and up orbitals
  real(REAL64), allocatable         ::  evd(:)                           !  default eigenvalue

! finds dimensions

  open(unit=ioread, file=trim(filein), status='OLD', form='FORMATTED')

  call atom_kb_test_input_size(ioread, iopsd, filepsd, nr, norb, lmax)

  close(unit=ioread)

! give enough room for excitations

  mxdl = lmax + 2
  mxdorb = norb + 4
  mxdnr = nr

! gets configuration

  open(unit=ioread, file=trim(filein), status='OLD', form='FORMATTED')

  allocate(ni(mxdorb))
  allocate(li(mxdorb))
  allocate(zu(mxdorb),zd(mxdorb))
  allocate(evd(mxdorb))

  call atom_kb_test_read_input(nameat, ititle, icorr, ispp,              &
      nval, ni, li, zd, zu, evd,                                         &
      ioread, mxdorb)

  close(unit=ioread)

! calculates the ae structure

   call atom_kb_test_scf(etotal,                                         &
      nameat, ititle, icorr, ispp,                                       &
      nval, ni, li, zd, zu, evd,                                         &
      iowrite, iopsd, filepsd,                                           &
      mxdnr, mxdorb, mxdl)

  deallocate(ni)
  deallocate(li)
  deallocate(zu,zd)
  deallocate(evd)

  return

end subroutine atom_kb_test_sub
