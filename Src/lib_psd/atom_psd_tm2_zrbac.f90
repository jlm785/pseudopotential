!>  Brackets the value of gamma
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.0.6
!>  \date         April 1990, 30 June 2021, 30 October 2021.
!>  \copyright    GNU Public License v2

subroutine atom_psd_tm2_zrbac(x1, x2, iflag,                             &
        jrc, r, drdi, bjin, rcpn, cdrc, lo, rpsi, delta, bj,             &
        mxdnr)

! **********************************************************
! *  njtj
! *    Routine brackets the root of the given function.
! *    Taken from Numerical Recipes page 245.
! *  njtj
! **********************************************************
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

  real(REAL64), intent(out)         ::  delta                            !<  value of delta
  integer, intent(out)              ::  iflag                            !<  iflag = 0: success
  real(REAL64), intent(out)         ::  bj(5)                            !<  polyr coefficients

! input and output

  real(REAL64), intent(inout)       ::  x1, x2                           !<  values of gamma

  real(REAL64), intent(inout)       ::  rpsi(mxdnr)                      !<  r*psi

! local variables

  real(REAL64)      ::  f1, f2                                           !  v''(0)

  integer    ::  ifl

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  real(REAL64), parameter    ::  FACTOR = 1.6
  integer, parameter         ::  NTRY = 50

! counters

  integer  ::  j

  if(x1 == x2) then
    iflag = 2
  else
    call atom_psd_tm2_gamfn(x1, f1, ifl,                                 &
        jrc, r, drdi, bjin, rcpn, cdrc, lo, rpsi, delta, bj,             &
        jrc)
    if(ifl /= 0) then
      iflag = 3
    else

      call atom_psd_tm2_gamfn(x2, f2, ifl,                               &
          jrc, r, drdi, bjin, rcpn, cdrc, lo, rpsi, delta, bj,           &
          jrc)
      if(ifl /= 0) then
        iflag = 4
      else

        iflag = 1
        do j = 1,ntry

          if(f1*f2 < ZERO) then
            iflag = 0

            exit

          endif

          if(abs(f1) < abs(f2)) then
            x1 = x1 + factor*(x1-x2)
            call atom_psd_tm2_gamfn(x1, f1, ifl,                         &
               jrc, r, drdi, bjin, rcpn, cdrc, lo, rpsi, delta, bj,      &
               mxdnr)
          else
            x2 = x2 + factor*(x2-x1)
            call atom_psd_tm2_gamfn(x2, f2, ifl,                          &
                jrc, r, drdi, bjin, rcpn, cdrc, lo, rpsi, delta, bj,      &
                mxdnr)
          endif

        enddo

      endif
    endif
  endif

! failure

  if(iflag /= 0) then
    write(6,*)
    write(6,*)
    write(6,'("  error in atom_psd_tm2_zrbac - cannot bracket gamma.  l = ",i3)') lo
  endif

  return

end subroutine atom_psd_tm2_zrbac
