!>  Carries out a Jacobi rotation on rows i and i+1 of r and qt matrices.
!>  Matrix r is triangular.  Based on numerical recipes

subroutine broyden_rotate(r,qt,n,np,i,a,b)

! Written 1 June 2018.  J.L.Martins
! Documentation 3 November 2020. JLM
! Copyright Jose Luis Martins/INESC-MN

! version 4.98

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  n                                !<  active space dimension
  integer, intent(in)               ::  np                               !<  leading dimension of a

  integer, intent(in)               ::  i                                !<  rotate 1,i+1

  real(REAL64), intent(in)          ::  a                                !<  rotation parameter cos(theta) = a / sqrt(a^2+b^2)
  real(REAL64), intent(in)          ::  b                                !<  rotation parameter sin(theta) = b / sqrt(a^2+b^2)

! input and output

  real(REAL64), intent(inout)       ::  qt(np,np)                        !<  q transposed
  real(REAL64), intent(inout)       ::  r(np,np)                         !<  r

! local

  real(REAL64)              ::  costh, sinth
  real(REAL64)              ::  fact, w, y

! constant

  real(REAL64), parameter   ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64

! counters

  integer         ::  j

  if(a == ZERO)then
    costh = ZERO
    sinth = sign(UM,b)
  else if(abs(a) > abs(b))then
    fact = b/a
    costh = sign(UM/sqrt(UM+fact*fact),a)
    sinth = fact*costh
  else
    fact = a/b
    sinth = sign(UM/sqrt(UM+fact*fact),b)
    costh = fact*sinth
  endif

  do j = i,n
    y = r(i,j)
    w = r(i+1,j)
    r(i,j) = costh*y - sinth*w
    r(i+1,j) = sinth*y + costh*w
  enddo

  do j = 1,n
    y = qt(i,j)
    w = qt(i+1,j)
    qt(i,j) = costh*y - sinth*w
    qt(i+1,j) = sinth*y + costh*w
  enddo

  return
end subroutine broyden_rotate

