!>  Integrates the relativistic 2 component "Dirac" equation.
!>  Finds the eigenvalue ev, the major and minor component
!>  of the wavefunction, ar and br.
!>  It uses an intial guess for the eigenvalues from dsolv1.
!>
!>  \author       Sverre Froyen, Norm Troullier, Jose Luis Martins
!>  \version      6.0.3
!>  \date         1980s,  12 September 2021.
!>  \copyright    GNU Public License v2

subroutine atom_atm_difrel_bound(nr, r, drdi, v, ar, br,                 &
          n, l, iso, znuc, vzero, ev, iflag, tol,                        &
          iowrite, mxdnr)

! converted to f90, March 2018
! cleanup and new interface, July 2019. JLM
! use of atm_integ_rel, 23 June 2021. JLM
! remove fixed evi. 28 June 2021. JLM
! so->iso, 12 September 2021. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points

  integer, intent(in)               ::  nr                               !<  number of radial points
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  real(REAL64), intent(in)          ::  v(mxdnr)                         !<  effective potential (RYDBERG)

  integer, intent(in)               ::  n                                !<  principal quantum number n
  integer, intent(in)               ::  l                                !<  angular quantum number l
  integer, intent(in)               ::  iso                              !<  2*spin or 2*(j-l)

  real(REAL64), intent(in)          ::  znuc                             !<  nuclear charge
  real(REAL64), intent(in)          ::  vzero                            !<  effective potential in Rydberg

  real(REAL64), intent(in)          ::  tol                              !<  precision of the solution tol ~1.0e-10

! output

  real(REAL64), intent(out)         ::  ar(mxdnr)                        !<  major component radial wave-function u = rR(r)
  real(REAL64), intent(out)         ::  br(mxdnr)                        !<  minor component (integral ar**2 + br**2 = 1)

  integer, intent(out)              ::  iflag                            !<  iflag = 0: success; iflag = 1: failed to converge; iflag > 3: major error

! input and output

  real(REAL64), intent(inout)       ::  ev                               !<  orbital energy, guess on input, accurate on output

! local variables

  integer                   ::  itmax, icmax

  integer                   ::  juflow

  integer                   ::  nctp, ninf                               !  classical turning point and practical infinity


  real(REAL64)              ::  s                                        !  r*psi ~ r**s
  integer                   ::  ka                                       !  -l  or  l+1
  real(REAL64)              ::  az                                       !  alfa*znuc

  real(REAL64)              ::  fanctp

  real(REAL64)              ::  alf
  real(REAL64)              ::  a1, a2
  real(REAL64)              ::  arin, arout, arpin, arpout
  real(REAL64)              ::  b0, b1, b2
  real(REAL64)              ::  dev, dev0
  real(REAL64)              ::  emin, emax
  real(REAL64)              ::  einf
  real(REAL64)              ::  evold
  real(REAL64)              ::  factor
  integer                   ::  ll

  integer                   ::  nodes
  real(REAL64)              ::  temp

  logical                   ::  lmany, lfew, lold                        !  repeated too many/few nodes/old had right number of nodes

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  ETOL = 1.0E-7_REAL64

  real(REAL64), parameter    ::  AI = 137.035999177_REAL64
  real(REAL64), parameter    ::  AI2 = 4*AI*AI

!------Machine dependent parameter-
!------Require exp(-2*expzer) to be within the range of the machine

  real(REAL64), parameter    ::  EXPZER = log10(huge(2*ONE))

! counters

  integer  ::  j, it, icount

  itmax = 300
  icmax = 100                                                            !  maximum iteration steps for ninf, nctp

  iflag = 1

  lmany = .TRUE.
  lfew = .TRUE.
  lold = .FALSE.

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
! an = (az (v0 - e) a(n-1) - (s + n + ka) (v0 - e - 4 ai**2) b(n-1))
!   / (2 n ai (2 s + n))
! bn = ((v0 - e) a(n-1) - 2 znuc an ) / (2 ai (s + n + ka))

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

! set the underflow trap, error from Berkeley version,
! fixed by Troy Barbee, sqrt(expzer) should be expzer/2,
! 4/17/90.

  juflow=1
  do j = 2,nr
    juflow = j
    if (s*abs(log(r(j))) < EXPZER/2) exit
  enddo

  einf = v(nr) + l*(l+1) / (r(nr)*r(nr))

  emax = einf
  emax = emax + abs(emax) / 2

  if(l == 0) then
    emin = ZERO
    do j = 2,nr
      temp = v(j) + 2*znuc/r(j)
      if(emin > temp) emin = temp
    enddo
!   should be safe...
    emin = emin - 2*znuc*znuc/(n*ONE) - ONE
  else
    emin = v(nr) + l*(l+1) / (r(nr)*r(nr))
    do j = 2,nr
      temp = v(j) + l*(l+1) / (r(j)*r(j))
      if(emin > temp) emin = temp
    enddo
  endif

  if(emin >= emax) then
    iflag = 5
    write(iowrite,*)
    write(iowrite,'(/,"   potential has no bound states")')

    return

  endif

  if (ev >= emax) ev = (emax+emin)/2
  if (ev <= emin) ev = (emax+emin)/2

! begin loop over trial energy

  do it = 1,itmax

    if (ev > ZERO) then
      write(iowrite,*)
      write(iowrite,'(" error in atm_difrel:   ev is greater than ",     &
          &  "v(infinty)  n = ",i2,",  l = ",i3,", j-l = ",f7.2," ev =", &
          &  e18.10)') n, l, iso*0.5, ev
      write(iowrite,*)
      iflag = 6

      exit

    endif

!   Find practical infinity ninf and classical turning
!   point nctp for orbital.


    do icount = 1,icmax

      do j = nr,2,-1
        temp = v(j) - ev
        if (temp < ZERO) temp = ZERO
        ninf = j

        if (r(j)*sqrt(temp) < EXPZER) exit

      enddo

      nctp = ninf
      do j = 2,ninf
        if (v(j) < ev) nctp = j
      enddo
      if (ev >= -ETOL*100) nctp = ninf - 5
      if (ev >= -ETOL) ev = ZERO

      if(nctp < 10) then

        ev = 0.9*ev

      elseif(nctp > ninf-5) then

        nctp = ninf-5
        do j = ninf-5,ninf
          if(ev > v(j)) ev = v(j)
        enddo

        exit

      else

        exit

      endif

    enddo

    if (nctp <= 6) then
      write(iowrite,*)
      write(iowrite,'("  error in atm_difrel: - cannot find the")')
      write(iowrite,'("  classical turning point for orbital  n = ",     &
           &  i2,",  l = ",i3,", j-l = ",f7.2)') n, l , iso*0.5
      write(iowrite,*)
      iflag = 3

      exit

    endif


!  Outward integration from 1 to nctp, startup.

    a1 = (az*(vzero-ev) - (s+1+ka)*(vzero-ev-AI2)*b0) / (2*AI*(2*s+1))
    b1 = ((vzero-ev)-2*znuc*a1) / (2*AI*(s+1+ka))
    a2 = (az*(vzero-ev)*a1 - (s+2+ka)*(vzero-ev-AI2)*b1) / (4*AI*(2*s+2))
    b2 = ((vzero-ev)*a1-2*znuc*a2) / (2*AI*(s+2+ka))
    do j = 2,5
      ar(j) = r(j)**s * (1 +(a1+a2*r(j))*r(j))
      br(j) = r(j)**s * (b0+(b1+b2*r(j))*r(j))
    enddo
    ar(1) = ZERO
    br(1) = ZERO

   call atom_atm_integ_rel(6, nctp, r, drdi, v, ev, ka, ar, br, fanctp,  &
         iowrite, mxdnr)

    nodes = 0
    do j = max(6,juflow+1),nctp
      if(ar(j)*ar(j-1) < ZERO) nodes = nodes + 1
    enddo

    arout = ar(nctp)
    arpout = fanctp

!   End outward integration.
!   If number of nodes correct, start inward integration
!   else modify energy stepwise and try again.

    if (nodes /= n-l-1) then
      if (nodes < n-l-1) then

!       too few nodes increase ev

        if (emin < ev) emin = ev
        if(lfew .or. .not. lold) then
          ev = ev - (ev-emax)/2
        else
          ev = (8*evold + 2*ev) / 10
        endif
        lfew = .TRUE.

      else

!       too many nodes decrease ev

        if (ev < emax) emax = ev
        if(lmany .or. .not. lold) then
          ev = ev - (ev-emin)/2
        else
          ev = (8*evold + 2*ev) / 10
        endif
        lmany = .TRUE.

      endif

      lold = .FALSE.

    else

      lmany = .FALSE.
      lfew = .FALSE.

!     Inward integration from ninf to nctp startup.

      do j = ninf,ninf-4,-1
        alf = v(j) - ev
        if (alf < ZERO) alf = ZERO
        alf = sqrt(alf)
        ar(j) = exp(-alf*r(j))
        br(j) = 2*AI*(alf+ka/r(j))*ar(j) / (v(j)-ev-AI2)
      enddo

   call atom_atm_integ_rel(ninf-5, nctp, r, drdi, v, ev, ka, ar, br, fanctp,  &
         iowrite, mxdnr)

      arin = ar(nctp)
      arpin = fanctp

!     End inward integration
!     Rescale ar and br outside nctp to match ar(nctp) from
!     outward integration.

      factor = arout/arin
      do j = nctp,ninf
        ar(j) = factor * ar(j)
        br(j) = factor * br(j)
      enddo
      arpin = factor * arpin

!     Find the normalizing factor.

      factor = ZERO
      ll = 4
      do j = 2,ninf
        factor = factor + ll*(ar(j)*ar(j)+br(j)*br(j))*drdi(j)
        ll = 6 - ll
      enddo
      factor = factor / 3

!     Modify the eigenvalue ev.

      dev0 = arout * (arpout-arpin) / (factor * drdi(nctp))
      if (5*abs(dev0) > (einf-ev)) then
        dev = sign((einf-ev),dev0)/5
      else
        dev = dev0
      endif
      evold = ev
      lold = .TRUE.

      if(dev > ZERO) then
        emin = ev - ETOL
      else
        emax = ev + ETOL
      endif

      ev = ev + dev
      if (ev > emax) ev = (evold + emax) / 2
      if (ev < emin) ev = (evold + emin) / 2

      if (abs(dev0) < tol*(ONE-ev)) then

        iflag = 0
        exit

      endif

    endif

  enddo                                                                  !  loop over it


  if(iflag > 0) then
    write(iowrite,*)
    write(iowrite,'("  atm_difrel: not converged after",i6," iterations")') it
    write(iowrite,*)
    write(iowrite,'("  n = ",i3,",  l = ",i3,", j-l = ",f7.2,",  ev =",  &
          & e18.10,",  nodes =",i2)') n, l, iso*0.5, ev, nodes

  endif

! Normalize the wavefunction (even if convergence is bad)

  factor = ZERO
  ll = 4
  do j = 2,ninf
    factor = factor + ll * (ar(j)*ar(j)+br(j)*br(j)) * drdi(j)
    ll = 6 - ll
  enddo
  factor = factor / 3


  if(factor > ZERO) then

    factor = ONE / sqrt(factor)
    do j = 1,ninf
      ar(j) = factor*ar(j)
      br(j) = factor*br(j)
    enddo

  else

!   paranoid case
    write(iowrite,*)
    write(iowrite,'("   atm_difrel:  bad normalization",g14.6)') factor
    iflag = 7

  endif

  return

end subroutine atom_atm_difrel_bound
