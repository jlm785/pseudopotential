!> DGESL solves a real general linear system A * X = B.

subroutine dgesl ( a, lda, n, ipvt, b, job )

!*****************************************************************************80
!
!! DGESL solves a real general linear system A * X = B.
!
!  Discussion:
!
!    DGESL can solve either of the systems A * X = B or A' * X = B.
!
!    The system matrix must have been factored by DGECO or DGEFA.
!
!    A division by zero will occur if the input factor contains a
!    zero on the diagonal.  Technically this indicates singularity
!    but it is often caused by improper arguments or improper
!    setting of LDA.  It will not occur if the subroutines are
!    called correctly and if DGECO has set 0.0 < RCOND
!    or DGEFA has set INFO == 0.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2001
!    12 October 2021.  Documentation, constants, remove dot_product.  José Luís Martins
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
!    Input, real ( kind = 8 ) A(LDA,N), the output from DGECO or DGEFA.
!
!    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
!
!    Input, integer ( kind = 4 ) N, the order of the matrix A.
!
!    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from DGECO or DGEFA.
!
!    Input/output, real ( kind = 8 ) B(N).
!    On input, the right hand side vector.
!    On output, the solution vector.
!
!    Input, integer ( kind = 4 ) JOB.
!    0, solve A * X = B;
!    nonzero, solve A' * X = B.
!
  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input and output

  integer, intent(in)               ::  lda                              !<  the leading dimension of A
  integer, intent(in)               ::  n                                !<  the order of the matrix A

  real(REAL64), intent(in)          ::  a(lda,n)                         !<  factorized matrix
  real(REAL64), intent(inout)       ::  b(n)                             !<  right hand side or solution

  integer, intent(in)               ::  ipvt(n)                          !<  pivot indices
  integer, intent(in)               ::  job                              !<  indicates transpose of A

! local variables

  real(REAL64)             ::  t
  real(REAL64), external   ::  ddot

! counters

  integer  ::  k, l

!  Solve A * X = B.

  if ( job == 0 ) then

    do k = 1, n - 1

      l = ipvt(k)
      t = b(l)

      if ( l /= k ) then
        b(l) = b(k)
        b(k) = t
      end if

      b(k+1:n) = b(k+1:n) + t * a(k+1:n,k)

    end do

    do k = n, 1, -1
      b(k) = b(k) / a(k,k)
      t = -b(k)
      b(1:k-1) = b(1:k-1) + t * a(1:k-1,k)
    end do

  else

!  Solve A' * X = B.

    do k = 1, n
      t = ddot(k-1,a(1,k),1,b(1),1)
      b(k) = ( b(k) - t ) / a(k,k)
    end do

    do k = n - 1, 1, -1

      b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
      l = ipvt(k)

      if ( l /= k ) then
        t = b(l)
        b(l) = b(k)
        b(k) = t
      end if

    end do

  end if

  return

end subroutine dgesl

