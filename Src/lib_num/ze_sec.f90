!>  Gets cpu time in seconds (sum over cores)
!>  linux gfortran/ifort version
!>
!>  \author       Jose Luis Martins
!>  \version      6.013
!>  \date         1980s, 22 June 2021
!>  \copyright    GNU Public License v2

subroutine zesec(tback)

! Written 21 october 2003
! Modified 11 September 2015. f90. JLM
! copyright INESC-MN/Jose Luis Martins

! version 4.60 pw
! version 6.00 atom

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)
  integer, parameter          :: REAL32 = selected_real_kind(6)

! output

  real(REAL64), intent(out)          ::  tback                           !<  cpu time in seconds

! local variables

  real(REAL32)    :: tin

  call cpu_time(tin)
  tback = tin

  return

end subroutine zesec
