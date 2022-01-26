!>  Finds the value of gamma for the v"(0)=0 criteria using bisection.
!>  rpsi(k) = r(k)^(l+1) * exp(polyr)
!>  polyr = delta + gamma^2 + sum_j b(j)*r^(2*j+2)
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.0.6
!>  \date         April 1990, 30 June 2021, 30 October 2021.
!>  \copyright    GNU Public License v2

subroutine atom_psd_tm2_rtbis(x1, x2, iflag,                             &
        jrc, r, drdi, bjin, rcpn, cdrc, lo, rpsi, gama, delta, bj,       &
        mxdnr)

! *************************************************************
! *  njtj                                                     *
! *  Finds the value of gamma for the v"(0)=0 criteria.       *
! *  The method used is bisection.  This routine              *
! *  was taken from Numerical Recipes, page 247.              *
! *  njtj                                                     *
! *************************************************************

! use bjin, instead of recalculating. 1 July 2021. JLM
! lp -> lo, 30 October 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension radial grid

  integer, intent(in)               ::  jrc                              !<  r(jrc) = r_cut

  integer, intent(in)               ::  lo                               !<  angular momentum l
  real(REAL64), intent(in)          ::  cdrc                             !<  int_0^r_cut rho(r) dr

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  real(REAL64), intent(in)          ::  bjin(5)                          !<  right member of the equation
  real(REAL64), intent(in)          ::  rcpn(12)                         !<  r_cut^n

! output

  real(REAL64), intent(out)         ::  gama                             !<  correct value of gamma
  real(REAL64), intent(out)         ::  delta                            !<  correct value of delta
  integer, intent(out)              ::  iflag                            !<  iflag = 0: success
  real(REAL64), intent(out)         ::  bj(5)                            !<  correct polyr coefficients

! input and output

  real(REAL64), intent(inout)       ::  x1, x2                           !<  floor and ceiling values of gamma

  real(REAL64), intent(inout)       ::  rpsi(mxdnr)                      !<  r*psi

! local variables

  real(REAL64)      ::  f, fmid                                           !  v''(0)
  real(REAL64)      ::  xmid, dx

  integer    ::  ifl

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  integer, parameter         ::  NTRY = 80
  real(REAL64), parameter    ::  XACC = 1.0E-10_REAL64

! counters

  integer  ::  j


  call atom_psd_tm2_gamfn(x1, f, ifl,                                    &
      jrc, r, drdi, bjin, rcpn, cdrc, lo, rpsi, delta, bj,               &
      mxdnr)
  if(ifl /= 0) then
    iflag = 3
  else

    call atom_psd_tm2_gamfn(x2, fmid, ifl,                               &
        jrc, r, drdi, bjin, rcpn, cdrc, lo, rpsi, delta, bj,             &
        mxdnr)
    if(ifl /= 0) then
      iflag = 4
    else

      if(f*fmid >= ZERO) then
        write(6,*)
        write(6,*)
        write(6,'("  error in bisection method, psd_tm2_rtbis")')
        write(6,'("  root must be bracketed")')
        iflag = 2
      else

        iflag = 1
        if(f < ZERO)then
          gama = x1
          dx = x2-x1
        else
          gama = x2
          dx = x1-x2
        endif

        do j = 1,NTRY
          dx = dx/2
          xmid = gama + dx
          call atom_psd_tm2_gamfn(xmid, fmid, ifl,                       &
              jrc, r, drdi, bjin, rcpn, cdrc, lo, rpsi, delta, bj,       &
              mxdnr)

          if(fmid < ZERO) gama = xmid
          if(abs(dx) < XACC .or. fmid == ZERO) then
            iflag = 0

            exit

          endif
        enddo

      endif
    endif
  endif

  if(iflag /=0) then
    write(6,*)
    write(6,*)
    write(6,'("  error in bisection method, psd_tm2_rtbis")')
    write(6,'("  too many bisections used")')
    iflag = 1
  endif

  return

end subroutine atom_psd_tm2_rtbis
