!>  Integrates the radial Schrödinger equation
!>  from origin up to revi for the input energy evi.
!>  It outputs the radial wave-function rpsi (r u(r))
!>  and the derivative drpsidr = d(rpsi)/dr.
!>
!>  \author       Sverre Froyen, Norm Troullier, Jose Luis Martins
!>  \version      6.013
!>  \date         28 June 2021
!>  \copyright    GNU Public License v2

subroutine atom_atm_difnrl_one(nr, r, drdi, d2rodr, v, rpsi, drpsidr,    &
      l, evi, revi, nrevi, iflag,                                        &
      iowrite, mxdnr)

! adapted from atm_difnrl, 28 June 2021. JLM
! jlm  version 6.013


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points

  integer, intent(in)               ::  nr                               !<  dimension of the number of radial points

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i
  real(REAL64), intent(in)          ::  d2rodr(mxdnr)                    !<  (d^2 r(i) / d i^2) / (d r / d i)

  real(REAL64), intent(in)          ::  v(mxdnr)                         !<  effective potential (RYDBERG)

  integer, intent(in)               ::  l                                !<  angular quantum number l

  real(REAL64), intent(in)          ::  evi                              !<  orbital energy
  real(REAL64), intent(in)          ::  revi                             !<  wave-function calculated up to revi

! output

  real(REAL64), intent(out)         ::  rpsi(mxdnr)                      !<  radial wave-function u = rR(r)  (integral rpsi**2 = 1)
  real(REAL64), intent(out)         ::  drpsidr(mxdnr)                   !<  d rpsi / d r.
!                                                                            In most of the code it is d rpsi / d i, and it is only converted before returning

  integer, intent(out)              ::  iflag                            !<  iflag = 0: success; iflag = 1: failed to converge, iflag > 3: major error

  integer, intent(out)              ::  nrevi                            !<  wave-functions was calculated up to r(nrevi)

! allocatable work arrays

  real(REAL64), allocatable         ::  drdisq(:)                        !  drdi*drdi

! local variables

  real(REAL64)              ::  aa, bb
  real(REAL64)              ::  factor
  integer                   ::  ll

  integer                   ::  joflow
  integer                   ::  jbig
  real(REAL64)              ::  temp

  real(REAL64)              ::  w2, w5, rv0, rv0pr

  real(REAL64)              ::  vrpsi0, vzero

!  real(REAL64)              ::  rpsimax

  real(REAL64)              ::  zeff                                     !  effective coulomb charge

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  EPS = 1.0E-8_REAL64
  real(REAL64), parameter    ::  BIG = 0.01_REAL64*huge(ONE)

! counters

  integer  ::  j


  iflag = 0

  allocate(drdisq(mxdnr))

!  lp = l+1

! initialize psi, initialize startup of outward integration.

  do j = 1,nr
    rpsi(j) = ZERO
  enddo
  do j = 1,nr
    drpsidr(j) = ZERO
  enddo
  if (l == 0) then
    drpsidr(1) = drdi(1)
  endif

  do j = 1,nr
    drdisq(j) = drdi(j)*drdi(j)
  enddo


! determine effective charge and vzero for startup of
! outward integration
! rpsi = r**(l+1) * (1 + aa r + bb r**2 + ... )
! aa = -znuc / (l+1)     bb = (-2 znuc aa + v(0) - e)/(4 l + 6)

  w2 = r(2) *  v(2) - (l*(l+1))/r(2)
  w5 = r(5) *  v(5) - (l*(l+1))/r(5)
  rv0pr = (w5-w2)/(r(5)-r(2))
  rv0 = w2 - rv0pr*r(2)

  zeff = -rv0/2

  if(abs(zeff - ONE*nint(zeff)) < EPS) zeff = ONE*nint(zeff)

  aa = -zeff/(l+1)
  vzero = -2*zeff*aa + rv0pr

  vrpsi0 = ZERO
  if (l == 0) vrpsi0 = -2*zeff
  if (l == 1) vrpsi0 = 2*ONE

! finds r(i) ~ revi

  nrevi = nr
  do j = nr,2,-1

      if (r(j) < revi) exit

      nrevi = j

  enddo

! avoids overflow

  do j = nrevi,2,-1
    temp = v(j) - evi
    if (temp < ZERO) temp = ZERO
    joflow = j

    if (r(j)*sqrt(temp) < 40.0) exit

  enddo
  nrevi = joflow

! avoids uncontroled oscillation

  do j = 2,nrevi
    temp = evi - v(j)
    if (temp < ZERO) temp = ZERO
    joflow = j

    if (drdi(j)*sqrt(temp) > 0.7) exit

  enddo
  nrevi = joflow

  if (nrevi <= 6) then
    write(iowrite,'(//,"error in atm_difnrl_one - nrevi = ",i10,          &
       &  ",   l = ",i3,",  evi = ",g12.3)') nrevi, l, evi
    iflag = 5

    return

  endif

! outward integration from 1 to nrevi
! startup

  bb = (vzero-evi)/(4*l+6)
  do j = 2,5
    rpsi(j) = r(j)**(l+1) * (ONE+(aa+bb*r(j))*r(j))
    drpsidr(j) = drdi(j) * r(j)**l * ((l+1)+(aa*(l+2)+bb*(l+3)*r(j))*r(j))
  enddo

  call atom_atm_integ_nr(6, nrevi, drdisq, d2rodr,                       &
      v, evi, rpsi, drpsidr, vrpsi0,                                     &
      iowrite, mxdnr)

! gets a normalization of whatever has been found
! both of the following tricks to avoid overflow work

!   rpsimax = ZERO
!   do j = 2,nrevi
!     rpsimax = max(rpsimax,abs(rpsi(j)))
!   enddo
!
!   do j = 2,nrevi
!     rpsi(j) = rpsi(j) / rpsimax
!     drpsidr(j) = drpsidr(j) / rpsimax
!   enddo

  jbig = 2
  factor = ZERO
  ll = 4
  do j = 2,nrevi
    factor = factor + ll*rpsi(j)*rpsi(j)*drdi(j)

    if(abs(factor) > BIG) exit

    ll = 6 - ll
    jbig = j
  enddo
  factor = factor / 3

  if(jbig /= nrevi) then
    do j = jbig+1,nrevi
      rpsi(j) = ZERO
      drpsidr(j) = ZERO
    enddo
    nrevi = jbig
    iflag = 3
  endif

! normalize wavefunction and change drpsidr from d(rpsi)/dj to d(rpsi)/dr

  if(factor > ZERO) then
    factor = ONE / sqrt(factor)
    do j = 1,nrevi
      rpsi(j) = factor*rpsi(j)
      drpsidr(j) = factor*drpsidr(j) / drdi(j)
    enddo
  else

!   paranoid case
    write(iowrite,*)
    write(iowrite,'("   atm_difnrl_one:  bad normalization",g14.6)') factor

    iflag = 7

  endif

  deallocate(drdisq)

  return
  end subroutine atom_atm_difnrl_one
