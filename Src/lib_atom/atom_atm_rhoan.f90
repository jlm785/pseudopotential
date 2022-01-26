!>  Write charge density to some file (iowrite).
!>  Not called by default.
!>
!>  \author       Sverre Froyen, Norm troullier, Jose Luis Martins
!>  \version      6.0.3
!>  \date         1980's 22 June 2021, 12 September 2021. JLM
!>  \copyright    GNU Public License v2

subroutine atom_atm_rhoan(nr, r, cdv, iowrite, mxdnr)


! converted to f90, March 2018
! cleanup and new interface, July 2019. JLM
! cdv, 12 September 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points

  integer, intent(in)               ::  nr                               !<  dimension of the number of radial points

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(in)          ::  cdv(mxdnr,-1:1)                  !<  4*pi*r**2 * valence charge density

! counters

  integer     ::  i

  if(iowrite > 0) then
    do i = 2,nr
      write(iowrite,'(3x,2e20.12)') r(i), cdv(i, 0) / (r(i)*r(i))
    enddo
  endif

  return

end subroutine atom_atm_rhoan
