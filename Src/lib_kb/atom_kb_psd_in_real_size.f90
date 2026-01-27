!<  Finds the size of the arrays from the parsec output of the atomic program
!>
!>  \author       J.L.Martins
!>  \version      6.1.0
!>  \date         25 May 2012, 17 September 2021.
!>  \copyright    GNU Public License v2

subroutine atom_kb_psd_in_real_size(iotape, fname, lmax, nrmax)

! nrmax->mxdnr. 17 September 2021. JLM
! parsec_size -> real_size in name. 26 January 2026. JLM


  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)

! input:

  integer, intent(in)               ::  iotape                           !<  io tape number
  character(len=*), intent(in)      ::  fname                            !<  file name of the output of the atomic program, parsec style

! output:

  integer, intent(out)              ::  lmax                             !<  maximum number of angular momentum components
  integer, intent(out)              ::  nrmax                            !<  maximum number of radial grid points.

! local variables

  integer                      ::  npotd, npotu, nr, nrm, l
  real(REAL64)                 ::  dummy

! counters

  integer                  :: i, j


  open(unit = iotape, file = trim(fname), status='OLD', form='FORMATTED')

  read(iotape,*)
  read(iotape,*)
  read(iotape,*)
  read(iotape,'(1x,2i3,i5)') npotd, npotu, nrm

  nr = nrm+1
  nrmax = nr

  read(iotape,*)                                  !'Radial grid follows'
  read(iotape,'(4g20.12)') (dummy,j=2,nr)

  lmax = 0
  do i = 1, npotd
    read (iotape,*)                               !'Pseudopotential follows (l on next line)'
    read(iotape,'(1x,i2)') l
    lmax = max(lmax,l)
    read(iotape,'(4g20.12)') (dummy,j=2,nr)
  enddo

  do i = 1, npotu
    read(iotape,*)                                !   'Minor comp. of Pseudop. follows (l on next line)'
    read(iotape,'(1x,i2)') l
    lmax = max(lmax,l)
    read(iotape,'(4g20.12)') (dummy,j=2,nr)
  enddo

  close(unit=iotape)

  return

end subroutine atom_kb_psd_in_real_size
