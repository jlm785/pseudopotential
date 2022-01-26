!>  dsolv2 finds the non-relativistic wave-function using
!>  difnrl to integrate the Scroedinger equation or
!>  the relativistic wave-function using
!>  difrel to integrate the Dirac equation.
!>  The energy level from the previous iteration or atom_atm_dsolv1
!>  is used as initial guess, and it must therefore be reasonably
!>  accurate.
!>
!>  \author       Sverre Froyen, Jose Luis Martins
!>  \version      6.0.3
!>  \date         1980s, 12 September 2021
!>  \copyright    GNU Public License v2

  subroutine atom_atm_dsolv2(iter, iconv, ispp, ifcore,                  &
       nr, r, drdi, d2rodr,                                              &
       norb, ncore, no, lo, iso, zo, znuc,                               &
       cdv, cdc, vionic, vhxc_orb, ev, ek, ep, evi,                      &
       iowrite, mxdnr, mxdorb, mxdl)

! converted to fortran 90, March 1st 2018
! cleanup and new interface, July 2019. JLM
! exit bug. 22 June 2021. JLM
! so->iso, vionic, single orbital in orban, cdv, 12 September 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdorb                           !<  dimension of the number of orbitals
  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points
  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum l

  integer, intent(in)               ::  iter                             !<  iteration number
  integer, intent(in)               ::  iconv                            !<  convergence flag (if iconv = 1, calculates Hartree energy)
  character(len=1), intent(in)      ::  ispp                             !<  spin polarization  ' ', 's', 'r'
  integer, intent(in)               ::  ifcore                           !<  0 no partial core correction, 1 partial xc, 2 partial

  integer, intent(in)               ::  nr                               !<  dimension of the number of radial points

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i
  real(REAL64), intent(in)          ::  d2rodr(mxdnr)                    !<  (d^2 r(i) / d i^2) / (d r / d i)

  integer, intent(in)               ::  ncore                            !<  number of orbitals treated as core

  integer, intent(in)               ::  norb                             !<  number of orbitals

  integer, intent(in)               ::  no(mxdorb)                       !<  principal quantum number n
  integer, intent(in)               ::  lo(mxdorb)                       !<  angular quantum number l
  integer, intent(in)               ::  iso(mxdorb)                      !<  2*spin or 2*(j-l)
  real(REAL64), intent(in)          ::  zo(mxdorb)                       !<  orbital occupation

  real(REAL64), intent(in)          ::  znuc                             !<  nuclear charge

  real(REAL64), intent(in)          ::  vionic(mxdnr,0:mxdl,-1:1)        !<  r*ionic potential in Rydberg; -1:  j=l-1/2 or s= -1/2.  0:  average or spinless.  1:  j=l+1/2 or s=1/2

  real(REAL64), intent(in)          ::  vhxc_orb(mxdnr,mxdorb)           !<  screening potential in Rydberg for each orbital

  real(REAL64), intent(in)          ::  evi(mxdorb)                      !<  fixed orbital energy

! output

  real(REAL64), intent(out)         ::  cdv(mxdnr,-1:1)                  !<  4*pi*r**2 * valence charge density

  real(REAL64), intent(out)         ::  cdc(mxdnr)                       !<  4*pi*r**2 * core charge density

  real(REAL64), intent(out)         ::  ek(mxdorb)                       !<  orbital kinetic energy
  real(REAL64), intent(out)         ::  ep(mxdorb)                       !<  orbital potential energy

! input and output

  real(REAL64), intent(inout)       ::  ev(mxdorb)                       !<  orbital energy

! allocatable work arrays

  real(REAL64), allocatable         ::  v(:)                             !  work array (potential)
  real(REAL64), allocatable         ::  ar(:), br(:)                     !  work array (wave-functions)

! local variables

  integer        ::  lp
  integer        ::  llp
  integer        ::  iflag

  real(REAL64)   ::  denr
  real(REAL64)   ::  vzero

  real(REAL64)   ::  revi
  integer        ::  nrevi

  logical        ::  lsp, lnosp

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  SMEV = 0.0001_REAL64
  real(REAL64), parameter    ::  TOL = 1.0E-12_REAL64

! counters

  integer     ::  i, j, k


  allocate(v(mxdnr))
  allocate(ar(mxdnr), br(mxdnr))

  revi = 10.0

! Initialize arrays for charge density


  do i = 1,nr
    cdv(i,-1) = ZERO
    cdv(i, 0) = ZERO
    cdv(i, 1) = ZERO
  enddo
  if (ifcore == 0) then
    do i = 1,nr
      cdc(i)= ZERO
    enddo
  endif

! Start the loop over orbitals.
! Note that spin zero is treated as down.

  do i = 1,norb

    if (no(i) > 0 .and. (zo(i) /= 0.0 .or. iconv /= 0)) then

      if (ev(i) >= 0.0) ev(i) = -SMEV

!     Set up the potential, set the wave functionc array to zero-ar.

      lp  = lo(i)+1
      llp = lo(i)*lp

      do j = 1,nr
        ar(j)=ZERO
      enddo

      do j = 2,nr
        v(j) = vionic(j,lo(i),iso(i))/r(j) + vhxc_orb(j,i)
      enddo
      vzero = vhxc_orb(1,i)

      if (ispp /= 'r') then
        do j = 2,nr
          v(j) = v(j) + llp/(r(j)*r(j))
        enddo
      endif

!     Call the integration routine.

      if (ispp /= 'r') then

        if(evi(i) == ZERO) then

          call atom_atm_difnrl(nr, r, drdi, d2rodr, v, ar, br,           &
              no(i), lo(i), ev(i), iflag, TOL,                           &
              iowrite, mxdnr)

        else

          call atom_atm_difnrl_one(nr, r, drdi, d2rodr, v, ar, br,       &
              lo(i), evi(i), revi, nrevi, iflag,                         &
              iowrite, mxdnr)
          ev(i) = evi(i)

        endif

      else

        if(evi(i) == ZERO) then

          call atom_atm_difrel(nr, r, drdi, v, ar, br,                   &
              no(i), lo(i), iso(i), znuc, vzero, ev(i), iflag, TOL,      &
              iowrite, mxdnr)

        else

          call atom_atm_difrel_one(nr, r, drdi, v, ar, br,               &
              lo(i), iso(i), znuc, vzero, ev(i), revi, nrevi, iflag,     &
              iowrite, mxdnr)

        endif

      endif

      if(iflag /= 0) then
        write(iowrite,*)
        write(iowrite,'("  dsolv2: wave-function not found, iflag = ",   &
                  &  i3,"  iter = ",i3,"  iorb = ",i3)') iflag, iter, i
        write(iowrite,*)

        if(iflag > 4) stop

      endif

!     Add to the charge density.

      if (ispp == 'r') then

        do j = 1,nr
          denr = zo(i) *(br(j) * br(j) + ar(j) * ar(j))
          cdv(j,iso(i)) = cdv(j,iso(i)) + denr
        enddo

      else

        do  j = 1,nr
          denr = zo(i) * ar(j) * ar(j)
          cdv(j,iso(i)) = cdv(j,iso(i)) + denr
          enddo

      endif

      if (ifcore == 0 .and. i <= ncore) then

        do j = 1,nr
          denr = zo(i) * ar(j) * ar(j)
          cdc(j) = cdc(j) + denr
        enddo

      endif

!     Compute various quantitities if last iteration.

      if (iconv == 1) then

        call atom_atm_orban(ispp, ar, br, nr, r, drdi,                   &
            no(i), lo(i), zo(i), iso(i),                                 &
            vionic(:,lo(i),iso(i)), vhxc_orb(:,i),  ev(i), ek(i), ep(i), &
            iowrite, mxdnr, mxdl)

      endif

    endif

  enddo

! End loop over orbitals.

! redistribute density

  lsp = .FALSE.
  lnosp = .FALSE.
  do k = 1,norb
    if(iso(k) == 0) then
      lnosp = .TRUE.
    else
      lsp = .TRUE.
    endif
  enddo

! paranoid check

  if(lsp .and. lnosp) then
    write(6,*)
    write(6,*) '   Stopped in dsolv2.  Inconsistent spin structure'
    write(6,*)

    STOP

  endif

  if(lsp) then
    do i = 1,nr
      cdv(i, 0) = cdv(i,-1) + cdv(i, 1)
    enddo
  else
    do i = 1,nr
      cdv(i,-1) = cdv(i, 0) / 2
      cdv(i, 1) = cdv(i, 0) / 2
    enddo
  endif

  deallocate(v)
  deallocate(ar, br)

  return

end subroutine atom_atm_dsolv2
