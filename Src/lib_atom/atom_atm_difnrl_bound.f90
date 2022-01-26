!>  Finds the eigenvalue ev, the wavefunction rpsi (r u(r))
!>  and the derivative drpsidr = d(rpsi)/dr,
!>  For the boundary condition psi("infinity") -> 0.
!>  It uses an input guess for the eigenvalue (from atm_dsolv1 for example),
!>  and searches for bound states.
!>  It is only reliable for bound states  ev < v(nr)
!>
!>  \author       Sverre Froyen, Norm Troullier, Jose Luis Martins
!>  \version      6.0.2
!>  \date         1980s, 22 June 2021, 9 September 2021.
!>  \copyright    GNU Public License v2

subroutine atom_atm_difnrl_bound(nr, r, drdi, d2rodr, v, rpsi, drpsidr,  &
      n, l, ev, iflag, tol,                                              &
      iowrite, mxdnr)

! converted to f90, March 2018
! cleanup and new interface, July 2019. JLM
! jlm  version 6.011
! jlm  replaced ev by einf - ef in final test.  dev0. September 2021. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points

  integer, intent(in)               ::  nr                               !<  number of radial points

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i
  real(REAL64), intent(in)          ::  d2rodr(mxdnr)                    !<  (d^2 r(i) / d i^2) / (d r(i) / d i)

  real(REAL64), intent(in)          ::  v(mxdnr)                         !<  effective potential (RYDBERG) includes centrifugal potential

  integer, intent(in)               ::  n                                !<  principal quantum number n
  integer, intent(in)               ::  l                                !<  angular quantum number l

  real(REAL64), intent(in)          ::  tol                              !<  precision of the solution tol ~1.0e-10

! output

  real(REAL64), intent(out)         ::  rpsi(mxdnr)                      !<  r*radial wave function
  real(REAL64), intent(out)         ::  drpsidr(mxdnr)                   !<  d rpsi / d i

  integer, intent(out)              ::  iflag                            !<  iflag = 0:success; iflag = 1: failed to converge; iflag > 3: major error

! input and output

  real(REAL64), intent(inout)       ::  ev                               !<  orbital energy, guess on input, accurate on output

! allocatable work arrays

  real(REAL64), allocatable         ::  drdisq(:)                        !  drdi*drdi

! local variables

  integer                   ::  itmax, icmax                             !  maximum ite

  integer                   ::  juflow

  integer                   ::  nctp, ninf                               !  classical t

  real(REAL64)              ::  aa, bb
  real(REAL64)              ::  alf
  real(REAL64)              ::  arctp, brctp
  real(REAL64)              ::  dev, dev0
  real(REAL64)              ::  emin, emax
  real(REAL64)              ::  evold
  real(REAL64)              ::  factor
  integer                   ::  ll

  integer                   ::  nodes
  real(REAL64)              ::  temp, einf
  real(REAL64)              ::  w2, w5, rv0, rv0pr

  real(REAL64)              ::  vrpsi0, vzero

  real(REAL64)              ::  zeff                                     !  effective coulomb charge

  logical                   ::  lmany, lfew, lold                        !  repeated too many/few nodes/old had right number of nodes

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64

  real(REAL64), parameter    ::  ETOL = 1.0E-07_REAL64


  real(REAL64), parameter    ::  EXPZER = log10(huge(2*ONE))

! counters

  integer  ::  j, it, icount


  allocate(drdisq(mxdnr))

  itmax = 300                                                       !  maximum iteratio
  icmax = 100                                                       !  maximum iteratio

  iflag = 1

  lmany = .FALSE.
  lfew = .FALSE.
  lold = .FALSE.

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

! set underflow trap, error from Berkeley version,
! fixed by Troy Barbee sqrt(expzer) should be expzer/2
! 4/17/90
! modified 9/9/2021

  juflow = 1
  do j = 2,nr
    juflow = j
    if ((l+1)*abs(log(r(j))) < EXPZER/2) exit
  enddo

! njtj  *** end major modification  ***

! determine effective charge and vzero for startup of
! outward integration
! rpsi = r**(l+1) * (1 + aa r + bb r**2 + ... )
! aa = -znuc / (l+1)     bb = (-2 znuc aa + v(0) - e)/(4 l + 6)

  w2 = r(2) *  v(2) - (l*(l+1))/r(2)
  w5 = r(5) *  v(5) - (l*(l+1))/r(5)
  rv0pr = (w5-w2)/(r(5)-r(2))
  rv0 = w2 - rv0pr*r(2)

  zeff = -rv0/2

  if(abs(zeff - ONE*nint(zeff)) < 10*ETOL) zeff = ONE*nint(zeff)

  aa = -zeff/(l+1)
  vzero = -2*zeff*aa + rv0pr

  vrpsi0 = ZERO
  if (l == 0) vrpsi0 = -2*zeff
  if (l == 1) vrpsi0 = 2*ONE

! finds bounds for the energy, emin,emax

  einf = v(nr)
!  einf = max(v(nr),ZERO)
!  if(abs(einf) < 0.01) einf = ZERO

  emax = einf
  emax = emax + abs(emax) / 2

  if(l == 0 .and. abs(zeff) > tol) then
    emin = ZERO
    do j = 2,nr
      temp = v(j) + 2*zeff/r(j)
      if(temp < emin) emin = temp
    enddo
!   should be safe...
    emin = emin - 2*zeff*zeff/(n*n*ONE) - ONE
  else
    emin = v(nr)
    do j = 2,nr
      if(emin > v(j)) emin = v(j)
    enddo
  endif

! paranoid check

  if(emin >= emax) then
    iflag = 5
    write(iowrite,*)
    write(iowrite,'(/,"   potential has no bound states")')

    return

  endif

! guess energy is not a good estimate

  if (ev >= emax) ev = (emax+emin)/2
  if (ev <= emin) ev = (emax+emin)/2


! begin loop over trial energy

  do it = 1,itmax

    if (it > itmax-2) then
      write(iowrite,'(" n =",i3," l =",i3," ev =",e18.10," nodes =",i2)')  &
           n,l,ev,nodes
    endif

    if (ev > einf) then

      iflag = 4
      write(iowrite,*)
      write(iowrite,*)
      write(iowrite,'(" error in atm_difnrl:   ev = ",g18.10,            &
            & " is greater than v(infinty)",g18.10,"for n = ",i3,        &
            & "  l = ",i3)') ev,einf,n,l
      write(iowrite,*)

      exit

    endif

!   find practical infinity ninf and
!   classical turning point nctp for orbital

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
      if (ev > einf - ETOL*10) nctp = ninf - 5
      if (ev > einf - ETOL) ev = einf

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
      write(iowrite,'(//,"error in difnrl - cannot find the classical ", &
           & /," turning point for n = ",i3,"  l = ",i3)') n,l
      iflag = 3

      exit

    endif

!   outward integration from 1 to nctp
!   startup

    bb = (vzero-ev)/(4*l+6)
    do j = 2,5
      rpsi(j) = r(j)**(l+1) * (ONE+(aa+bb*r(j))*r(j))
      drpsidr(j) = drdi(j) * r(j)**l * ((l+1)+(aa*(l+2)+bb*(l+3)*r(j))*r(j))
    enddo

    call atom_atm_integ_nr(6, nctp, drdisq, d2rodr,                      &
        v, ev, rpsi, drpsidr, vrpsi0,                                    &
        iowrite, mxdnr)

!   count nodes - if no underflow

    nodes = 0
    do j = max(6,juflow+1),nctp
      if(rpsi(j)*rpsi(j-1) < ZERO) nodes = nodes + 1
    enddo

!   njtj  ***  end major modification  ***

    arctp = rpsi(nctp)
    brctp = drpsidr(nctp)

!   end outward integration

!   if number of nodes correct, start inward integration
!   else modify energy stepwise and try again

    if (nodes /= n-l-1) then

!     incorrect number of nodes

      if (nodes < n-l-1) then

!       too few nodes; increase ev

        if (emin < ev) emin = ev
        if(lfew .or. .not. lold) then
          ev = ev - (ev-emax)/2
        else
          ev = (8*evold + 2*ev) / 10
        endif
        lfew = .TRUE.

      else

!       too many nodes; decrease ev

        if (emax > ev) emax = ev
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

!     correct number of nodes
!     inward integration from ninf to nctp
!     startup

      do j = ninf,ninf-4,-1
        alf = v(j) - ev
        if (alf < ZERO) alf = ZERO
        alf = sqrt(alf)
        rpsi(j) = exp(-alf*r(j))
        drpsidr(j) = -drdi(j)*alf*rpsi(j)
      enddo

      call atom_atm_integ_nr(ninf-5, nctp, drdisq, d2rodr,               &
          v, ev, rpsi, drpsidr, vrpsi0,                                  &
          iowrite, mxdnr)

!     rescale rpsi and drpsidr outside nctp to match rpsi(nctp) from
!     outward integration

      factor = arctp/rpsi(nctp)
      do j = nctp,ninf
        rpsi(j) = factor * rpsi(j)
        drpsidr(j) = factor * drpsidr(j)
      enddo

!     find normalizing factor

      factor = ZERO
      ll = 4
      do j = 2,ninf
        factor = factor + ll*rpsi(j)*rpsi(j)*drdi(j)
        ll = 6 - ll
      enddo
      factor = factor / 3

!     modify eigenvalue ev

      dev0 = arctp * (brctp-drpsidr(nctp)) / (factor * drdi(nctp))
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

!     convergence with dev0

      if (abs(dev0) < tol*(ONE-ev)) then

        iflag = 0

        exit

      endif

    endif

  enddo                                                             !  loop over it


  if(iflag > 0) then
    write(iowrite,*)
    write(iowrite,'("   atm_difnrl:  solution was not found in ",i6,     &
            & " iterations ")') itmax
    write(iowrite,'(" n =",i3," l =",i3," ev =",e18.10," nodes =",i2)')  &
          n,l,ev,nodes
    write(iowrite,*)

!   gets a normalization of whatever has been found (probably scattering state)

    factor = ZERO
    ll = 4
    do j = 2,ninf
      factor = factor + ll*rpsi(j)*rpsi(j)*drdi(j)
      ll = 6 - ll
    enddo
    factor = factor / 3

  endif

! normalize wavefunction and change drpsidr from d(rpsi)/dj to d(rpsi)/dr

  if(factor > ZERO) then
    factor = ONE / sqrt(factor)
    do j = 1,ninf
      rpsi(j) = factor*rpsi(j)
      drpsidr(j) = factor*drpsidr(j) / drdi(j)
    enddo
  else

!   paranoid case
    write(iowrite,*)
    write(iowrite,'("   atm_difnrl:  bad normalization",g14.6)') factor
    iflag = 7

  endif

  deallocate(drdisq)

  return

  end subroutine atom_atm_difnrl_bound
