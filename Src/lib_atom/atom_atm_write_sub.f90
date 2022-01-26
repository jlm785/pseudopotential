!>  Writes a default input configuration file from the chemical symbol.
!>  The default is Perdew-Zunger non-spin-polarized LDA for the
!>  atomic ground state configuration.
!>  User may edit that file to have the other cases.
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.7
!>  \date         22 June 2021, 21 December 2021
!>  \copyright    GNU Public License v2

subroutine atom_atm_write_sub(ioread, filein, nameat)

! nval2, jhard. 21 December 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  ioread                           !<  default tape for reading
  character(len=*), intent(in)      ::  filein                           !<  name of default tape for reading

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbol of the atom

! local variables

  character(len=10)                 ::  ititle(5)                        !  title of calculation
  character(len=9)                  ::  bdate                            !  date
  character(len=2)                  ::  icorr                            !  correlation type
  character(len=1)                  ::  ispp                             !  spin polarization  ' ', 's', 'r'

  integer                           ::  iz                               !  atomic number

  real(REAL64)                      ::  rmax, aa, bb                     !  logarithmic mesh parameters

  character(len=2)                  ::  pcc

  integer                           ::  ncore                            !  canonical number of core orbitals
  integer                           ::  nval                             !  canonical number of interesting valence orbitals
  integer                           ::  jhard                            !  flag for accuracy/speed compromise
  integer                           ::  no(5), lo(5)                     !  configuration
  real(REAL64)                      ::  zo(5)                            !  occupation
  real(REAL64)                      ::  rc(5)                            !  core radii
  character(len=30)                 ::  status                           !  quality of
! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64

! counters

  integer     ::  n


  jhard = 0
  call atom_p_tbl_psd_pcc(nameat, pcc, jhard)

  call zedate(bdate)

  ititle(1) = '     '//nameat//'   '
  ititle(2) = '   file ge'
  ititle(3) = 'nerated by'
  ititle(4) = ' atm_write'
  ititle(5) = ' '//bdate

  icorr = 'ca'

  call atom_p_tbl_charge(nameat, iz)

  ispp = ' '
  if(iz > 20) ispp = 'r'

  rmax = 120*ONE
  aa = 4*ONE
  bb = 200*ONE

  jhard = 0
  call atom_p_tbl_config(nameat, ncore, nval, no, lo, zo, jhard)

  call atom_p_tbl_psd_tm2(nameat, rc, status)

  open(unit=ioread, file=trim(filein), status='NEW', form='FORMATTED')

  write(ioread,'(3x,a2,5a10)') pcc,ititle
  write(ioread,'(8x,a3)') 'tm2'
  write(ioread,'(" n=",a2," c=",a2,a1)') nameat,icorr,ispp

  write(ioread,'(6f10.1)') ONE*iz, ZERO, ZERO, rmax, aa, bb

  write(ioread,'(1x,i4,1x,i4)') ncore,nval

  do n = 1,nval
    write(ioread,'(1x,i4,1x,i4,1x,f8.2,2x,f8.2)') no(n), lo(n), zo(n), ZERO
  enddo

  write(ioread,'(1x,5(f7.2,2x))') (rc(n),n=1,nval)
  write(ioread,'(50x)')

  close(unit=ioread)

  return

end subroutine atom_atm_write_sub
