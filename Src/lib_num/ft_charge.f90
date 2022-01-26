!>  calculates a radial fourier transform.
!>  Four point formula.  4 pi r^2 in the input function
!>
!>  \author       Sverre Froyen, Norm Troullier, JL Martins
!>  \version      6.013
!>  \date         80s, April 2012
!>  \copyright    GNU Public License v2

subroutine ft_charge(nr, nql, r, drdi, qp, fin, fout)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  nr                               !<  number of points in the radial grid
  integer, intent(in)               ::  nql                              !<  number of points for the local Fourier grid

  real(REAL64), intent(in)          ::  r(0:nr)                          !<  radial grid points   r(i) = a*(exp(b*i)-1), i=1,...,nr  
  real(REAL64), intent(in)          ::  drdi(0:nr)                       !<  d r(i) / d i
  real(REAL64), intent(in)          ::  qp(0:nql)                        !<  fourier grid points

  real(REAL64), intent(in)          ::  fin(0:nr)                        !<  4*pi*r**2*function to be transformed

! output

  real(REAL64), intent(out)         ::  fout(0:nql)                      !<  Fourier transform of function fin

! allocated arrays

  real(REAL64), allocatable  ::  w(:), y(:)

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  
! counters

  integer     ::  j, k


  allocate(w(0:nr+4),y(0:nr+4))

  do j = 1,nr
    w(j) = drdi(j)*fin(j)/r(j)
  enddo
  w(0) = ZERO

  do k = 1,nql
    fout(k) = ZERO
    do j = 0,nr
      y(j) = w(j)*sin(qp(k)*r(j))
    enddo
    do j = 0,nr,4
      fout(k) = fout(k) + 7*y(j) + 32*y(j+1) + 12*y(j+2) + 32*y(j+3) + 7*y(j+4)
    enddo
    fout(k) = ((2*fout(k))/45) / qp(k)
  enddo

  deallocate(w,y)

  return

end subroutine ft_charge
