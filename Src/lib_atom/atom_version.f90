!>  Hard written version of the code
!>
!>  \author       Jose Luis Martins
!>  \version      6.2.1
!>  \date         20 August 2020, 17 November 2026
!>  \copyright    GNU Public License v2

subroutine atom_version(vers)

  implicit none

  character(len=5), intent(out)     ::  vers                             !<  version of the code

  vers = '6.2.1'

  return

end subroutine atom_version

