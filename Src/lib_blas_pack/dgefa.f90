!> DGEFA factors A=L*U, a real general matrix.

subroutine dgefa ( a, lda, n, ipvt, info )

!*****************************************************************************80
!
!! DGEFA factors a real general matrix.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2001
!    12 October 2021. Documentation, constants, remove cycle. José Luís Martins
!
!  Author:
!
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
!    LINPACK User's Guide,
!    SIAM, 1979,
!    ISBN13: 978-0-898711-72-1,
!    LC: QA214.L56.
!
!  Parameters:
!
!    Input/output, real ( kind = 8 ) A(LDA,N).
!    On intput, the matrix to be factored.
!    On output, an upper triangular matrix and the multipliers used to obtain
!    it.  The factorization can be written A=L*U, where L is a product of
!    permutation and unit lower triangular matrices, and U is upper triangular.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
!
!    Output, integer ( kind = 4 ) INFO, singularity indicator.
!    0, normal value.
!    K, if U(K,K) == 0.  This is not an error condition for this subroutine,
!    but it does indicate that DGESL or DGEDI will divide by zero if called.
!    Use RCOND in DGECO for a reliable indication of singularity.

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input and output

  integer, intent(in)               ::  lda                              !<  the leading dimension of A
  integer, intent(in)               ::  n                                !<  the order of the matrix A

  real(REAL64), intent(inout)       ::  a(lda,n)                         !<  matrix to be factored

  integer, intent(out)              ::  info                             !<  singularity indicator
  integer, intent(out)              ::  ipvt(n)                          !<  pivot indices

! local variables

  real(REAL64)       ::  t

! external functions

  integer, external  ::  idamax

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64

! counters

  integer  ::  j, k, l

!  Gaussian elimination with partial pivoting.

  info = 0

  do k = 1, n - 1

!  Find L = pivot index.

    l = idamax ( n-k+1, a(k,k), 1 ) + k - 1
    ipvt(k) = l

!  Zero pivot implies this column already triangularized.

    if ( a(l,k) == ZERO ) then
      info = k
    else

!    Interchange if necessary.

      if ( l /= k ) then
        t = a(l,k)
        a(l,k) = a(k,k)
        a(k,k) = t
      end if

!    Compute multipliers.

      t = -ONE / a(k,k)
      a(k+1:n,k) = a(k+1:n,k) * t

!    Row elimination with column indexing.

      do j = k + 1, n
        t = a(l,j)
        if ( l /= k ) then
          a(l,j) = a(k,j)
          a(k,j) = t
        end if
        a(k+1:n,j) = a(k+1:n,j) + t * a(k+1:n,k)
      end do

    endif

  end do

  ipvt(n) = n

  if ( a(n,n) == ZERO ) then
    info = n
  end if

  return

end subroutine dgefa
