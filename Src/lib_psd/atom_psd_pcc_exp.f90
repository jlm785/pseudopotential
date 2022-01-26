!>  constructs the partial core correction. core charge,
!>  first and second derivatives are conserved at the
!>  point icore. interpolation is done with Lagrange formula.
!>  the core function is exp(ac+bc*r**2+cc*r**4)
!>
!>  \author       Manuel Maria Alemany
!>  \version      6.013
!>  \date         January 2000, 5 July 2021
!>  \copyright    GNU Public License v2

subroutine atom_psd_pcc_exp(nr, icore, ac, bc, cc, r, cdc)

! Converted to fortran90. 5 July 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

!  integer, intent(in)               ::  mxdnr                            !<  dimension  of the number of radial points
  integer, intent(in)               ::  nr                               !<  number of radial points

  real(REAL64), intent(in)          ::  r(nr)                         !<  radial grid

  integer, intent(in)               ::  icore                            !<  partial core correction for r < r(icore)

! output

  real(REAL64), intent(out)         ::  ac, bc, cc                       !<  core function is exp(ac+bc*r**2+cc*r**4)

! input and output

  real(REAL64), intent(inout)       ::  cdc(nr)                          !<  core charge density

! local variables

  integer, parameter       ::  nn = 5
  real(REAL64)             ::  dr_k(-nn:nn), ddr_k(-nn:nn)
  real(REAL64)             ::  cdc_sca(-nn:nn)
  integer                  ::  in1, in2
  real(REAL64)             ::  f1, f2
  real(REAL64)             ::  dr, ddr
  real(REAL64)             ::  cdcp,cdcpp,ddcdc

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64

! counters

  integer     ::  i, in, jn, kn, ln


  if ( (icore - nn) < 1 .or. (icore + nn) > nr ) then
    write(6,*)
    write(6,*)
    write(6,*)  ' error in psd_pcc_exp : '
    write(6,*)  ' derivatives at r(icore) not calculated'
    write(6,*)  ' icore, nn, nr = ', icore, nn, nr

    stop

  endif

  do i = -nn,nn
    cdc_sca(i)= log(cdc(i+icore)) - 2*log(r(i+icore))
  enddo

  in1 =-nn
  in2 = nn

  do in = in1,in2
    if (in == 0) then
      dr_k(in) = ZERO
      do jn = in1,in2
        if (jn /= 0) dr_k(in) = dr_k(in) + ONE/(0 - jn)
      enddo
    else
      f1 = ONE
      f2 = ONE
      do jn = in1,in2
        if (jn /= in .and. jn /= 0) f1 = f1 * (0  - jn)
        if (jn /= in)               f2 = f2 * (in - jn)
      enddo
      dr_k(in) = f1 / f2
    endif
  enddo

  dr = ZERO
  do in = in1,in2
    dr = dr + dr_k(in) * r(icore+in)
  enddo

  cdcp = ZERO
  do in  = in1,in2
    cdcp = cdcp + dr_k(in) * cdc_sca(in)
  enddo
  cdcp = cdcp / dr

  do in = in1,in2
    ddr_k(in)= ZERO
    if (in == 0) then
      do jn = in1,in2
        if (jn /= 0) then
          do kn = jn+1,in2
            if (kn /= 0) ddr_k(in) = ddr_k(in) + (2*ONE) / ((0 - jn)*(0 - kn))
          enddo
        endif
      enddo
    else
      f2 = ONE
      do jn = in1,in2
        if (jn /= in) then
          f2 = f2 * (in - jn)
          do kn = jn+1,in2
            if (kn /= in) then
              f1 = 2*ONE
              do ln = in1,in2
                if (ln /= in .and. ln /= jn .and. ln /= kn) f1 = f1 * (0 - ln)
              enddo
              ddr_k(in) = ddr_k(in) + f1
            endif
          enddo
        endif
      enddo
      ddr_k(in) = ddr_k(in) / f2
    endif
  enddo

  ddr   = ZERO
  ddcdc = ZERO
  do in = in1,in2
    ddr   = ddr   + ddr_k(in) * r(icore+in)
    ddcdc = ddcdc + ddr_k(in) * cdc_sca(in)
  enddo

  cdcpp = (ddcdc - ddr * cdcp) / dr**2

  cc = r(icore) * cdcpp - cdcp
  cc = cc / (8 * r(icore)**3)

  bc = cdcp - 4 * cc * r(icore)**3
  bc = bc / (2 * r(icore))

  ac = cdc_sca(0) - bc*r(icore)**2 - cc*r(icore)**4

  do i = 1,icore
      cdc(i)= r(i)*r(i) *             &
              exp( (cc*r(i)*r(i) + bc) * r(i)*r(i) + ac )
  enddo

  return

end subroutine atom_psd_pcc_exp
