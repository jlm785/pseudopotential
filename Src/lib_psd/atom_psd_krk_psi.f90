!>  Calculates r*psi
!>  rpsi(k) = r(k)^(l+1) * exp(polyr)
!>  polyr = sum_j b(j)*r(k)^j,   j=0,1,2,3,4,... (BUT miss 1 for invertibility, b(1) = 0)
!>  OR  j = 0,2,4,6,8....
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.0.6
!>  \date         1990s, 30 October 2021
!>  \copyright    GNU Public License v2

subroutine atom_psd_krk_psi(n, bkrk, lsq, jrc, r, lo, rpsi, mxdnr)

! generate pseudo wavefunction
! generalization of atom_psd_tm2_psi

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension  of the number of radial points

  integer, intent(in)               ::  n                                !<  order of polynomial
  real(REAL64), intent(in)          ::  bkrk(0:n)                        !<  coefficient of polynomial
  logical, intent(in)               ::  lsq                              !<  polynomial on squares

  integer, intent(in)               ::  jrc                              !<  r(jrc) = r_cut

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid

  integer, intent(in)               ::  lo                               !<  angular momentum l

! input and output

  real(REAL64), intent(inout)       ::  rpsi(mxdnr)                      !<  r*psi, modified up to jrc.

! local variables

  real(REAL64)       ::  r2, polyr

! counters

  integer  ::  k, j


  do k = 1,jrc

    if(lsq) then
      r2 = r(k)*r(k)
    else
      r2 = r(k)
    endif
    polyr = bkrk(n)
    if(n > 0) then
      do j = n-1,0,-1
        polyr = polyr*r2 + bkrk(j)
      enddo
    endif

!   r**(l+1)

    rpsi(k) = r(k)
    if(lo > 0) then
      do j = 1,lo
        rpsi(k) = rpsi(k)*r(k)
      enddo
    endif

    rpsi(k) = rpsi(k) * exp(polyr)

  enddo

  return

end subroutine atom_psd_krk_psi
