!>  Finds the dimensions of arrays to read datafile.dat
!>
!>  \author       Sverre Froyen,  Norm Troullier, Jose Luis Martins
!>  \version      6.0.3
!>  \date         May 1, 1991, April 2018
!>  \copyright    GNU Public License v2

subroutine atom_psd_dat_in_size(ioae, fileae, lmax, nr, norb)

!  ***********************************************************
!  *                                                         *
!  *  The routine reads data from file 'datafile.dat'        *
!  *  for latter use                                         *
!  *                                                         *
!  *  Reverse of datout.f Version dated May 1, 1991          *
!  *                                                         *
!  ***********************************************************

!  converted to f90, April 2018
!  jlm  version 5.805
!  lmax.  14 September 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


  integer, intent(in)               ::  ioae                             !<  default tape for writing
  character(len=*), intent(in)      ::  fileae                           !<  name of default tape for all-electron results

! output

  integer, intent(out)              ::  nr                               !<  dimension of the number of radial points

  integer, intent(out)              ::  norb                             !<  dimension of the number of orbitals
  integer, intent(out)              ::  lmax                             !<  maximum angular momentum

! local variables

  integer               ::  itype
  character(len=2)      ::  icorr
  character(len=1)      ::  ispp

  character(len=2)      ::  nameat
  integer               ::  lmaxp1

! Open and read out data to current file datafile.dat.

  open (unit=ioae, file=trim(fileae), status='old',form='unformatted')

  read(ioae) itype, icorr, ispp, nr

  read(ioae)
  read(ioae)
  read(ioae) lmaxp1, nameat, norb
! convention of old files
  lmax = lmaxp1-1

  close (unit=ioae)

  return

end subroutine atom_psd_dat_in_size
