!>  Returns the maximum of the exponent polyr of the pseudo
!>  wave-function rpsi(k) = r(k)^(l+1) * exp(polyr) (see atom_psd_krk_psi)
!>  for r(1),...,r(jrc).  Used to reject trial coefficients before
!>  exp(polyr) overflows.
!>
!>  \author       claude+J.L.Martins
!>  \version      6.1.0
!>  \date         26 October 2026.
!>  \copyright    GNU Public License v2

subroutine atom_psd_krk_polymax(n, bkrk, lsq, jrc, r, polymax, mxdnr)

! written from atom_psd_krk_psi. 26 October 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension  of the number of radial points

  integer, intent(in)               ::  n                                !<  order of polynomial
  real(REAL64), intent(in)          ::  bkrk(0:n)                        !<  coefficient of polynomial
  logical, intent(in)               ::  lsq                              !<  polynomial on squares

  integer, intent(in)               ::  jrc                              !<  r(jrc) = r_cut

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid

! output

  real(REAL64), intent(out)         ::  polymax                          !<  maximum of polyr for r(1),...,r(jrc)

! local variables

  real(REAL64)       ::  r2, polyr

! counters

  integer  ::  k, j


  polymax = -huge(polymax)

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

    polymax = max(polymax, polyr)

  enddo

  return

end subroutine atom_psd_krk_polymax
