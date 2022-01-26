!>  calculates a fourier transform of a radial function times a
!>  spherical harmonic. Eight point formula
!>
!>  \author       Sverre Froyen, Norm Troullier, JL Martins
!>  \version      6.013
!>  \date         80s, April 2012
!>  \copyright    GNU Public License v2

subroutine ft_precise(nr, nqnl, r, drdi, qp, fin, fout, l)

! adapted from old code J.L.Martins 30/4/2012
! copyright jlm inesc-mn

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  nr                               !<  number of points in the radial grid
  integer, intent(in)               ::  nqnl                             !<  number of points for the non-local Fourier grid

  real(REAL64), intent(in)          ::  r(0:nr)                          !<  radial grid points   r(i) = a*(exp(b*i)-1), i=1,...,nr  
  real(REAL64), intent(in)          ::  drdi(0:nr)                       !<  d r(i) / d i
  real(REAL64), intent(in)          ::  qp(0:nqnl)                       !<  fourier grid points

  real(REAL64), intent(in)          ::  fin(0:nr)                        !<  function to be transformed
  integer, intent(in)               ::  l                                !<  angular momentum

! output

  real(REAL64), intent(out)         :: fout(0:nqnl)                      !<  Fourier transform of function

! allocated arrays

  real(REAL64), allocatable  ::  w(:), y(:)

! local variables

  real(REAL64)  ::  crnorm

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  real(REAL64), parameter    ::  PI4 = 16*atan(ONE)
  
! counters

  integer     :: j, k

! functions

  real(REAL64), external   :: sbessj

  allocate(w(0:nr+7),y(0:nr+7))

  do k = 1,nr
    w(k) = drdi(k)*r(k)*r(k)*fin(k)
  enddo
  y(0) = zero
  do k = nr+1,nr+7
    y(k) = zero
  enddo
  crnorm = (7*sqrt((2*l+1)*PI4))/17280

! due to the high number of oscilations in the integrand,
! an eight point newton-cotes intagration method is used.
! see  abramowitz and stegun eq. 25.4.17

  do j = 0,nqnl

    do k = 1,nr
      y(k) = sbessj(l,qp(j)*r(k))*w(k)
    enddo

    fout(j) = zero
    do k=0,nr,7
      fout(j) = fout(j)+751*(y(k)+y(k+7))+3577*(y(k+1)+y(k+6))+          &
       1323*(y(k+2)+y(k+5))+2989*(y(k+3)+y(k+4))
    enddo
    fout(j) = crnorm * fout(j)

  enddo

  deallocate(w,y)
  
  return

end subroutine ft_precise
