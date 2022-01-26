!>  QR update

  subroutine broyden_qr_update(r,qt,n,np,u,v)

! Written May/June 2018.  J.L.Martins
! Documentation 3 November 2020. JLM
! Copyright Jose Luis Martins/INESC-MN

! version 4.98

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  n                                !<  active space dimension
  integer, intent(in)               ::  np                               !<  leading dimension of a

  real(REAL64), intent(in)          ::  v(np)                            !<  v vector  

! input and output

  real(REAL64), intent(inout)       ::  qt(np,np)                        !<  q transposed
  real(REAL64), intent(inout)       ::  r(np,np)                         !<  r matrix

  real(REAL64), intent(inout)       ::  u(np)                            !<  u vector    

! constant

  real(REAL64), parameter   ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64

! counters

  integer         ::  i, j, k

  do k = n,1,-1
    j = k
    if(u(k) /= ZERO) exit
  enddo
  k = j

  if(k > 1) then
    do i = k-1,1,-1
      call broyden_rotate(r,qt,n,np,i,u(i),-u(i+1))
      if(u(i) == ZERO)then
        u(i) = abs(u(i+1))
      else if(abs(u(i)) > abs(u(i+1))) then
        u(i) = abs(u(i))*sqrt(UM + (u(i+1)/u(i))*(u(i+1)/u(i)))
      else
        u(i) = abs(u(i+1))*sqrt(UM + (u(i)/u(i+1))*(u(i)/u(i+1)))
      endif
    enddo
  endif

  do j = 1,n
    r(1,j)=r(1,j) + u(1)*v(j)
  enddo

  if(k > 1) then
    do i = 1,k-1
      call broyden_rotate(r,qt,n,np,i,r(i,i),-r(i+1,i))
    enddo
  endif

  return
  end subroutine broyden_qr_update
