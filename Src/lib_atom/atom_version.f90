!>  Hard written version of the code
!>
!>  \author       Jose Luis Martins
!>  \version      6.2.0
!>  \date         20 August 2020, 8 November 2026
!>  \copyright    GNU Public License v2

subroutine atom_version(vers)

  implicit none

  character(len=5), intent(out)     ::  vers                             !<  version of the code

  vers = '6.2.0'

  return

end subroutine atom_version

