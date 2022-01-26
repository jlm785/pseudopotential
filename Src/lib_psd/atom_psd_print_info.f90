!>  Reads the basic information from the  all-electron results file
!>  and prints it.
!>
!>  \author       Sverre Froyen,  Norm Troullier, Jose Luis Martins
!>  \version      6.014
!>  \date         1980s and 1990s, April 2018, 27 August 2021.
!>  \copyright    GNU Public License v2

subroutine atom_psd_print_info(iowrite, ioae, fileae, nameat)

!  Adapted from atom_psd_dat_in


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer , intent(in)              ::  iowrite                          !<  default output tape

  integer, intent(in)               ::  ioae                             !<  default tape for writing
  character(len=*), intent(in)      ::  fileae                           !<  name of default tape for all-electron results

! output

  character(len=2), intent(out)     ::  nameat                           !<  chemical symbot of the atom

! local variables

  integer             ::  ioerror
  logical             ::  lold                                           !  indicates it is an  old format

  character(len=9)     ::  bdate                                         !  date the subroutine was called
  character(len=5)     ::  vers                                          !  version of the code

  integer              ::  itype                                         !  type of calculation (-1 signals end of calculation)
  character(len=2)     ::  icorr                                         !  correlation type
  character(len=1)     ::  ispp                                          !  spin polarization  ' ', 's', 'r'

  real(REAL64)         ::  a
  real(REAL64)         ::  b

  integer              ::  nr

  integer              ::  norb
  integer              ::  ncore

  integer              ::  lmax


! Open and read out data (default: from  datafile.dat).

  open (unit=ioae, file=trim(fileae), status='old', form='unformatted')

  read(ioae, iostat = ioerror ) itype, icorr, ispp, nr, a, b, bdate, vers
  lold = .FALSE.
  if(ioerror /= 0) lold = .TRUE.

  read(ioae)
  read(ioae)

  read(ioae) lmax, nameat, norb, ncore

  close (unit=ioae)

  write(6,*)
  write(6,*)  '  Found a file with all-electron results for ',nameat
  write(6,*)

  if(iowrite /= 6) then
    write(iowrite,*)
    write(iowrite,*)  '  Found a file with all-electron results for ',nameat
    write(iowrite,*)
  endif

  write(iowrite,*)
  if (icorr == 'pb') then
    write(iowrite,*) '  Calculation was done with PBE GGA'
  elseif (icorr== 'pw') then
    write(iowrite,*) '  Calculation was done with Perdew-Wang LDA'
  elseif (icorr == 'ca') then
    write(iowrite,*) '  Calculation was done with Ceperley-Alder LDA'
  else
    write(iowrite,*) '  Unrecognized correlation choice'
  endif
  write(iowrite,*)

  write(iowrite,*)
  if (ispp == 'r') then
    write(iowrite,*) '  Based on a relativistic calculation'
  elseif (ispp == 's') then
    write(iowrite,*) '  Based on a spin-polarized calculation'
  else
    write(iowrite,*) '  Based on a non-spin-polarized and non-relativistic calculation'
  endif
  write(iowrite,*)

  if(.not. lold) then
  write(iowrite,*)
  write(iowrite,*) '  Calculation was run on ',bdate,' with code version ',vers
  write(iowrite,*)
  endif

  return

end subroutine atom_psd_print_info
