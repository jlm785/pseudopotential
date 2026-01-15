!------------------------------------------------------------!
! This file is distributed as part of the cpw2000 code and   !
! under the terms of the GNU General Public License. See the !
! file `LICENSE' in the root directory of the cpw2000        !
! distribution, or http://www.gnu.org/copyleft/gpl.txt       !
!                                                            !
! The webpage of the cpw2000 code is not yet written         !
!                                                            !
! The cpw2000 code is hosted on GitHub:                      !
!                                                            !
! https://github.com/jlm785/cpw2000                          !
!------------------------------------------------------------!

!>  Accepts a string of nchar characters and replaces
!>  any lowercase letters by uppercase ones.
!>  If n = 0, uses the string length.
!>
!>  \author       JLM adapted from the web
!>  \version      5.12
!>  \date         before August 1997, 22 November 2025.
!>  \copyright    GNU Public License v2

subroutine chrcap(string,nchar)

! Modified documentation August 2019.  JLM
! Indentation case. 22 November 2025. JLM

  implicit none

  integer, intent(in)                ::  nchar                           !<  string length
  character(len=*), intent(inout)    ::  string                          !<  string to be uppercased

  integer       ::  ncopy, i, itemp

  ncopy = nchar
  if(ncopy < 1) ncopy = len(string)
  do i=1,ncopy

    if(lge(string(i:i),'a') .and. lle(string(i:i),'z'))then
      itemp=ichar(string(i:i))+ichar('A')-ichar('a')
      string(i:i)=char(itemp)
      endif
  enddo

  return

end subroutine chrcap

