!>  Generates the screened pseudopotential from the wave-function
!>  rpsi(k) = r(k)^(l+1) * exp(polyr)
!>  polyr = sum_j b(j)*r(k)^j,   j=0,2,3,4,... (missing 1 for invertibility)
!>  OR  j = 0,2,4,6,8....
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.0.6
!>  \date         April 1990, 30 October 2021
!>  \copyright    GNU Public License v2

subroutine atom_psd_krk_v(n, bkrk, lsq, jrc, r, lo, ev, vj, mxdnr)

! generate pseudo-potential
! generalization of atom_psd_tm2_v.

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
  real(REAL64), intent(in)          ::  ev                               !<  orbital energy

! output

  real(REAL64), intent(out)         ::  vj(mxdnr)                        !<  v_pseudo_screened

! local variables

  real(REAL64)       ::  r2, xlambda

! counters

  integer  ::  j, k


  if(n == 0 .or. (n == 1 .and. .not. lsq)) then

    do k = 1,jrc
      vj(k) = ev
    enddo

  else
    if(lsq) then

      do k = 1,jrc

        r2 = r(k)*r(k)

        xlambda = n*bkrk(n)
        if(n > 1) then
          do j = n-1,1,-1
            xlambda = xlambda*r2 + j*bkrk(j)
          enddo
        endif
        xlambda = 2*xlambda

        vj(k) = ev + xlambda * (2*(lo+1) + xlambda*r2)

        xlambda = (2*n)*(2*n-1)*bkrk(n)
        if(n > 1) then
          do j = n-1,1,-1
            xlambda = xlambda*r2 + (2*j)*(2*j-1)*bkrk(j)
          enddo
        endif
        vj(k) = vj(k) + xlambda
      enddo

    else

      do k = 1,jrc

        xlambda = n*bkrk(n)
        if(n > 2) then
          do j = n-1,2,-1
            xlambda = xlambda*r(k) + j*bkrk(j)
          enddo
        endif

        vj(k) = ev + xlambda * (2*(lo+1) + xlambda*r(k)*r(k))

        xlambda = n*(n-1)*bkrk(n)
        if(n > 2) then
          do j = n-1,2,-1
            xlambda = xlambda*r(k) + j*(j-1)*bkrk(j)
          enddo
        endif
        vj(k) = vj(k) + xlambda
      enddo

    endif

  endif

  return

end subroutine atom_psd_krk_v
