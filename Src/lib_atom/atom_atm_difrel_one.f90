!>  Integrates the relativistic 2 component "Dirac" equation
!>  from origin up to revi for the input energy evi.
!>  Outputs the major and minor component
!>  of the wavefunction, ar and br. (r u(r))
!>
!>  \author       Sverre Froyen, Norm Troullier, Jose Luis Martins
!>  \version      6.0.3
!>  \date         28 June 2021,12 September 2021.
!>  \copyright    GNU Public License v2

subroutine atom_atm_difrel_one(nr, r, drdi, v, ar, br,                   &
          l, iso, znuc, vzero, evi, revi, nrevi, iflag,                  &
          iowrite, mxdnr)

! adapted from atm_difrel. 28 June 2021. JLM
! so->iso, 12 September 2021. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points

  integer, intent(in)               ::  nr                               !<  dimension of the number of radial points
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  real(REAL64), intent(in)          ::  v(mxdnr)                         !<  effective potential (RYDBERG)

  integer, intent(in)               ::  l                                !<  angular quantum number l
  integer, intent(in)               ::  iso                              !<  spin quantum number

  real(REAL64), intent(in)          ::  znuc                             !<  nuclear charge
  real(REAL64), intent(in)          ::  vzero                            !<  effective potential in Rydberg

  real(REAL64), intent(in)          ::  evi                              !<  orbital energy
  real(REAL64), intent(in)          ::  revi                             !<  wave-function should be calculated up to revi


! output

  real(REAL64), intent(out)         ::  ar(mxdnr)                        !<  major component radial wave-function u = rR(r)
  real(REAL64), intent(out)         ::  br(mxdnr)                        !<  minor component (integral ar**2 + br**2 = 1)

  integer, intent(out)              ::  iflag                            !<  iflag = 0: success; iflag = 1: failed to converge; iflag > 3: major error

  integer, intent(out)              ::  nrevi                            !<  wave-functions was calculated up to r(nrevi)

! local variables

  real(REAL64)              ::  s                                        !  r*psi ~ r**s
  integer                   ::  ka                                       !  -l  or  l+1
  real(REAL64)              ::  az                                       !  alfa*znuc

  real(REAL64)              ::  b0, b1, b2
  real(REAL64)              ::  a1, a2
  real(REAL64)              ::  factor
  integer                   ::  ll
  integer                   ::  joflow
  integer                   ::  jbig
  real(REAL64)              ::  temp

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  BIG = 0.01*HUGE(ONE)

  real(REAL64), parameter    ::  AI = 137.0359991_REAL64
  real(REAL64), parameter    ::  AI2 = 4*AI*AI

! counters

  integer  ::  j

  iflag = 0

  az = znuc / AI
  if (iso == -1 .and. l /= 0) then
    ka=-l
  else
    ka = l+1
  endif

! determine effective charge and vzero for startup of
! outward integration
! ar = r**s * (1  + a1 r + a2 r**2 + ... )
! br = r**s * (b0 + b1 r + b2 r**2 + ... )
! s = sqrt (ka**2 - az**2)    b0 = - az / (s + ka)
! an = (az (v0 - e) a(n-1) - (s + n + ka) (v0 - e - 4*ai**2) b(n-1))
!   / (2 n ai (2 s + n))
! bn = ((v0 - e) a(n-1) - 2 znuc an ) / ( 2 ai (s + n + ka))

  s = sqrt(ka*ka-az*az)

  if (ka > 0) then
    b0 = -az/(s+ka)
  else
    b0 = (s-ka)/az
  endif

  do j = 1,nr
    ar(j) = ZERO
    br(j) = ZERO
  enddo


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
    write(iowrite,'(//,"error in atm_difrel_one - nrevi = ",i10,          &
       &  ",   l = ",i3,", j-l = ",f8.2,",  evi = ",g12.3)') &
             nrevi, l, iso*0.5, evi
    iflag = 5

    return

  endif

!  Outward integration from 1 to nrevi

    a1 = (az*(vzero-evi) - (s+1+ka)*(vzero-evi-AI2)*b0) / (2*AI*(2*s+1))
    b1 = ((vzero-evi)-2*znuc*a1) / (2*AI*(s+1+ka))
    a2 = (az*(vzero-evi)*a1 - (s+2+ka)*(vzero-evi-AI2)*b1) / (4*AI*(2*s+2))
    b2 = ((vzero-evi)*a1-2*znuc*a2) / (2*AI*(s+2+ka))
    do j = 2,5
      ar(j) = r(j)**s * (1 +(a1+a2*r(j))*r(j))
      br(j) = r(j)**s * (b0+(b1+b2*r(j))*r(j))
    enddo
    ar(1) = ZERO
    br(1) = ZERO

   call atom_atm_integ_rel(6, nrevi, r, drdi, v, evi, ka, ar, br, temp,  &
         iowrite, mxdnr)

! Normalize the wavefunction of whatever has been found
! Uses trick to avoid overflow

  jbig = 2
  factor = ZERO
  ll = 4
  do j = 2,nrevi
    factor = factor + ll * (ar(j)*ar(j)+br(j)*br(j)) * drdi(j)

    if(abs(factor) > BIG) exit

    ll = 6 - ll
    jbig = j
  enddo
  factor = factor / 3

  if(jbig /= nrevi) then
    do j = jbig+1,nrevi
      ar(j) = ZERO
      br(j) = ZERO
    enddo
    nrevi = jbig
    iflag = 3
  endif

  if(factor > ZERO) then

    factor = ONE / sqrt(factor)
    do j = 1,nrevi
      ar(j) = factor*ar(j)
      br(j) = factor*br(j)
    enddo

  else

!   paranoid case
    write(iowrite,*)
    write(iowrite,'("   atm_difrel_one:  bad normalization",g14.6)') factor
    iflag = 7

  endif

  return

end subroutine atom_atm_difrel_one
