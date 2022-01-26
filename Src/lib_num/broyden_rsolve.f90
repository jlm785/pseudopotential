!>  Solves a set of upper triangular equations ut.x = b, 
!>  where ut is an upper triangular matrix stored in the following way:
!>  the off-diagonal elements are stored in a, the diagonal elements in d.
!>  b is overwritten with x on output.
!>  Based on numerical recipes.

subroutine broyden_rsolv(a,n,np,d,b)


! Written 1 June 2018.  J.L.Martins
! Documentation 3 November 2020. JLM
! Copyright Jose Luis Martins/INESC-MN

! version 4.98

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  n                                !<  active space dimension
  integer, intent(in)               ::  np                               !<  leading dimension of a

  real(REAL64), intent(in)          ::  a(np,np)                         !<  matrix decomposed
  real(REAL64), intent(in)          ::  d(n)                             !<  matrix decomposed (diagonal)

! input and output

  real(REAL64), intent(inout)       ::  b(n)                             !<  rigth handside in / solution out

! local

  real(REAL64)              ::  xsum

! constant

  real(REAL64), parameter   ::  ZERO = 0.0_REAL64

! counters

  integer         ::  i, j


  b(n) = b(n)/d(n)

  if(n > 1) then
    do i = n-1,1,-1
      xsum = ZERO
      do j = i+1,n
        xsum = xsum + a(i,j)*b(j)
      enddo
      b(i) = (b(i)-xsum)/d(i)
    enddo
  endif

  return
end subroutine broyden_rsolv
