!>  calculates a fourier transform of the local potential
!>  it deals separately with the coulomb tail
!>
!>  \author       Sverre Froyen, Norm Troullier, JL Martins
!>  \version      6.013
!>  \date         80s, May 2012
!>  \copyright    GNU Public License v2

subroutine ft_local(nr, nql, r, drdi, qp, zion, vlocal, vql0, vlocft)

! adapted from old code J.L.Martins 20/5/2012
! copyright jlm inesc-mn

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  nr                               !<  number of points in the radial grid
  integer, intent(in)               ::  nql                              !<  number of points for the local Fourier grid

  real(REAL64), intent(in)          ::  r(0:nr)                          !<  radial grid points   r(i) = a*(exp(b*i)-1), i=1,...,nr  
  real(REAL64), intent(in)          ::  drdi(0:nr)                       !<  d r(i) / d i
  real(REAL64), intent(in)          ::  qp(0:nql)                        !<  fourier grid points
  real(REAL64), intent(in)          ::  zion                             !<  ionic charge
  real(REAL64), intent(in)          ::  vlocal(0:nr)                     !<  r*local pseudopotential

! output

  real(REAL64), intent(out)         ::  vlocft(0:nql)                    !<  Fourier transform of local pseudopotential
  real(REAL64), intent(out)         ::  vql0                             !<  zero frequency of the integral without the Coulomb part

! allocated arrays

  real(REAL64), allocatable  ::  w(:), y(:)

! local variables

  real(REAL64)  ::  vt

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  real(REAL64), parameter    ::  PI4 = 16*atan(ONE)
  
! counters

  integer     :: j, k


  allocate(w(0:nr+4),y(0:nr+4))

! Fourier transform local potential, first for q=0.
! Use Bode's rule to integrate.

  vt = 2*zion 

  w(0) = ZERO
  do k = 1,nr
    w(k) = drdi(k)*r(k)*(vlocal(k)+vt)
  enddo
  do k = nr+1,nr+4
    w(k) = ZERO
  enddo

  vql0 = ZERO
  do k = 0,nr,4
    vql0 = vql0 + 7*w(k) + 32*w(k+1) + 12*w(k+2) + 32*w(k+3) + 7*w(k+4)
  enddo

  vql0 = (2 * vql0 * PI4) / 45

! fourier transform local potential, for rest q>0.
! use bode's rule to integrate.

  y(0) = ZERO
  do k = nr+1,nr+4
    y(k) = ZERO
  enddo

  do k = 1,nr
    w(k) = w(k)/r(k)
  enddo

  do j = 1,nql
    vlocft(j) = zero
    do k = 1,nr
      y(k) = sin(qp(j)*r(k))*w(k)
    enddo
    do k = 0,nr,4
      vlocft(j) = vlocft(j) + 7*y(k) + 32*y(k+1) + 12*y(k+2) +           &
                              32*y(k+3) + 7*y(k+4)
    enddo
    vlocft(j) = PI4 * ((2*vlocft(j))/45 - vt/qp(j)) / qp(j)
  enddo

  deallocate(w,y)

  return
  
end  subroutine ft_local
