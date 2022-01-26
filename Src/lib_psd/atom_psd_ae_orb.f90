!>  Calculates the all-electron wave-functions, and respective screened
!>  potential, for preparation of the pseudopotential construction.
!>
!>  \author       Sverre Froyen,  Norm Troullier, Jose Luis Martins
!>  \version      6.0.6
!>  \date         1980s and 1990s, 30 June 2021, 31 October 2021.
!>  \copyright    GNU Public License v2

subroutine atom_psd_ae_orb(ispp, nr, r, drdi, d2rodr, v, ar, br,         &
     no, lo, iso, znuc, ev, flgev,                                       &
     vionic, vhxc,                                                       &
     iowrite, mxdnr)

! extracted from pseudo2. 30 June 2021. JLM
! so-> iso, vionic, vhxc, 15, 23 September 2021.
! rextr = 0. 28 October 2021. JLM
! atom_atm_zero_extr. 31 october 2021. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  iowrite                          !<  default tape number

  integer, intent(in)               ::  mxdnr                            !<  dimension  of the number of radial points

  character(len=1), intent(in)      ::  ispp                             !<  spin polarization  ' ', 's', 'r'

  integer, intent(in)               ::  nr                               !<  number of radial points

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i
  real(REAL64), intent(in)          ::  d2rodr(mxdnr)                    !<  (d^2 r(i) / d i^2) / (d r(i) / d i)

  integer, intent(in)               ::  no                               !<  principal quantum number n
  integer, intent(in)               ::  lo                               !<  angular quantum number l
  integer, intent(in)               ::  iso                              !<  2*spin or 2*(j-l)
  real(REAL64), intent(in)          ::  znuc                             !<  nuclear charge

  logical, intent(in)               ::  flgev                            !<  if true use fixed energy

  real(REAL64), intent(in)          ::  vionic(mxdnr)                    !<  r*ionic potential in Rydberg

  real(REAL64), intent(in)          ::  vhxc(mxdnr)                      !<  screening potential (Ry)

! output

  real(REAL64), intent(out)         ::  ar(mxdnr)                        !<  r*psi or major component
  real(REAL64), intent(out)         ::  br(mxdnr)                        !<  d ar / d r or minor component

  real(REAL64), intent(out)         ::  v(mxdnr)                         !<  effective potential

! input and output

  real(REAL64), intent(inout)       ::  ev                               !<  orbital energy, guess on input, accurate on output

! local variables

  real(REAL64)       ::  vzero

  integer            ::  iflag

  real(REAL64)       ::  revi
  integer            ::  nrevi

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  TOL = 1.0E-10_REAL64


! counters

  integer  ::  j


  revi = 10.0

  do j = 1,nr
    ar(j) = ZERO
  enddo

! set up the potential

  do j = 2,nr
    v(j) = vionic(j)/r(j) + vhxc(j)
  enddo
  vzero = vhxc(1)

  if (ispp /= 'r') then
    do j = 2,nr
      v(j) = v(j) + ((lo+1)*lo) / (r(j)*r(j))
    enddo
  endif

! solves Schrodinger or Dirac equation

  iflag = 0
  if (ispp /= 'r') then

    if(flgev) then
      call atom_atm_difnrl_one(nr, r, drdi, d2rodr, v, ar, br,           &
          lo, ev, revi, nrevi, iflag,                                    &
          iowrite, mxdnr)
    else
      call atom_atm_difnrl(nr, r, drdi, d2rodr, v, ar, br,               &
          no, lo, ev, iflag, TOL,                                        &
          iowrite, mxdnr)
    endif

  else
    if(flgev) then
      call atom_atm_difrel_one(nr, r, drdi, v, ar, br,                   &
          lo, iso, znuc, vzero, ev, revi, nrevi, iflag,                  &
          iowrite, mxdnr)
    else
      call atom_atm_difrel(nr, r, drdi, v, ar, br,                       &
          no, lo, iso, znuc, vzero, ev, iflag, TOL,                      &
          iowrite, mxdnr)
    endif

  endif

  if(iflag /= 0) then
    write(iowrite,*)
    write(iowrite,'("  psd_ae_orb: wave-function not found, iflag = ",   &
              &  i3)') iflag
    write(iowrite,*)

    if(iflag > 4) stop

  endif

  return

end subroutine atom_psd_ae_orb
