!>  Constructs the QR decomposition of a SQUARE matrix a.
!   The upper triangular matrix R is stored in thefollowing way:
!   the off-diagonal elements are stored in a, the diagonal elements in d.
!   The matrix Q is represented by a product of Householder matrices
!   Q_1,...,Q_n-1 Q_j = 1 - (1/c_j)|v_j><v_j| where the
!   v_j is stored in the j-th column of a.  Note that v_j(i) = 0 for j < i.

  subroutine broyden_qr_decomp(a,n,np,c,d,sing)

! Written May/June 2018.  J.L.Martins
! Documentation 3 November 2020. JLM
! Copyright Jose Luis Martins/INESC-MN

! version 4.98

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  n                                !<  active space dimension
  integer, intent(in)               ::  np                               !<  leading dimension of a

! output

  real(REAL64), intent(out)         ::  c(n)                             !<  prefactor of the "projector" c = <v|v>/2
  real(REAL64), intent(out)         ::  d(n)                             !<  diagonal of the R matrix

  logical,intent(out)               ::  sing                             !<  if the matrix is singular sing = .TRUE.

! input and output

  real(REAL64), intent(inout)       ::  a(np,np)                         !<  matrix to de decomposed/decomposition

! local

  real(REAL64)              ::  xscale, sigma, xsum, tau

! constant

  real(REAL64), parameter   ::  ZERO = 0.0_REAL64

! counters

  integer         ::  i, j, k

  sing = .FALSE.

  if(n > 1) then

    do k = 1,n-1
      xscale =  ZERO
      do i = k,n
        xscale = max(xscale,abs(a(i,k)))
      enddo

      if(xscale == ZERO)then

        sing = .TRUE.
        c(k)= ZERO
        d(k)= ZERO

      else
        do i = k,n
          a(i,k) = a(i,k)/xscale
        enddo

        xsum= ZERO
        do i = k,n
          xsum = xsum+a(i,k)*a(i,k)
        enddo
        sigma = sign(sqrt(xsum),a(k,k))
        a(k,k) = a(k,k) + sigma
        c(k) = sigma*a(k,k)
        d(k) =-xscale*sigma

        do j=k+1,n
          xsum =  ZERO
          do i = k,n
            xsum = xsum + a(i,k)*a(i,j)
          enddo
          tau = xsum/c(k)
          do i = k,n
            a(i,j) = a(i,j) - tau*a(i,k)
          enddo
        enddo

      endif

    enddo

  else
    c(1) = ZERO
  endif

  d(n) = a(n,n)
  if(d(n) == ZERO) sing=.TRUE.

  return
  end subroutine broyden_qr_decomp

