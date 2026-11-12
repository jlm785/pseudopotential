!>  Removes leading blanks and replaces multiple blanks
!>  by a single one in a string.
!>
!>  \author       claude
!>  \version      6.2.0
!>  \date         12 November 2026.
!>  \copyright    GNU Public License v2

subroutine chr_squeeze(string)

! Written 12 November 2026. claude

  implicit none

! input and output

  character(len=*), intent(inout)   ::  string                           !<  string to be squeezed

! local variables

  integer          ::  n

! counters

  integer          ::  i


  string = adjustl(string)

  n = 0
  do i = 1,len(string)
    if(string(i:i) /= ' ' .or. n == 0) then
      n = n + 1
      string(n:n) = string(i:i)
    elseif(string(n:n) /= ' ') then
      n = n + 1
      string(n:n) = string(i:i)
    endif
  enddo

  if(n < len(string)) string(n+1:) = ' '

  return

end subroutine chr_squeeze
