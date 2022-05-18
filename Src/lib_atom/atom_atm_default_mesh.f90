!>  Returns the default parameters of the quasi-logarithmic mesh
!>  Could be an include file, but this is old fortran...
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.8
!>  \date         21 April 2022
!>  \copyright    GNU Public License v2

subroutine atom_atm_default_mesh(rmax, aa, bb)

! Consistent use of default parameters.
! Bug reported by Raymond Atta-Fynn. 21 April 2022. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! output

  real(REAL64), intent(out)         ::  rmax                             !<  maximum mesh radius
  real(REAL64), intent(out)         ::  aa                               !<  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64), intent(out)         ::  bb                               !<  a = exp(-aa)/znuc, b = 1/bb

 ! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64

  rmax = 120*ONE
  aa = 6*ONE
  bb = 200*ONE

  return

end subroutine atom_atm_default_mesh
