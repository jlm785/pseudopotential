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

!>    indexes a real array by the heapsort method
!>    adapted from http://rosettacode.org
!>    see also W. H. Preuss et al. Numerical Recipes     

      subroutine sort(n,a,indx)

!     written 24 June 2013. JLM
!     Modified documentation August 2019.  JLM
!     copyright  J.L.Martins, INESC-MN.

!      version 4.94

      implicit none
      integer, parameter  :: REAL64 = selected_real_kind(12)

!     input

      integer, intent(in)        ::  n                                   !<  length of array
      real(REAL64), intent(in)   ::  a(n)                                !<  array to be indexed

!     output

      integer, intent(out)       ::  indx(n)                             !<  index of array a

!     local variables

      integer    ::  iroot,ichild,istart,ibot,indxt,ic

      if(n < 1) return
      
      do ichild=1,n
        indx(ichild) = ichild
      enddo

      if(n == 1) return
      
!     hiring phase

      ibot = n
      
      do istart = n/2,1,-1
        indxt = indx(istart)
        iroot = istart

!       long enough siftdown loop does not exceed ~log(n)/log(2)

        do ic = 1,n+5
          ichild = 2*iroot
          if(ichild <= ibot) then
            if(ichild < ibot) then
              if(a(indx(ichild)) < a(indx(ichild+1)))                    &
     &              ichild = ichild + 1
            endif

            if(a(indxt) < a(indx(ichild))) then
              indx(iroot) = indx(ichild)
              iroot = ichild
            else

              exit

            endif

          else

            exit

          endif
          
        enddo
        
        indx(iroot) = indxt

      enddo

!     retirement and promotion phase

      istart = 1
      
      do ibot = n-1,1,-1
        indxt = indx(ibot+1)
        indx(ibot+1) = indx(1)
        if(ibot == 1) then
            
          exit

        endif

        iroot = istart

!       long enough siftdown loop does not exceed ~log(n)/log(2)  (repeated...)
        
        do ic = 1,n+5
          ichild = 2*iroot
          if(ichild <= ibot) then
            if(ichild < ibot) then
              if(a(indx(ichild)) < a(indx(ichild+1)))                    &
     &              ichild = ichild + 1
            endif

            if(a(indxt) < a(indx(ichild))) then
              indx(iroot) = indx(ichild)
              iroot = ichild
            else

              exit

            endif

          else

            exit

          endif
          
        enddo
        
        indx(iroot) = indxt

      enddo

      indx(1) = indxt

      return

      end subroutine sort
