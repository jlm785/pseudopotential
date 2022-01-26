!>     kb_test_dsolv2 finds the wave function and energies
!>     for a kleinman-Bylander type of potential.
!>     The energy level and wave-function from the previous iteration is used
!>     as initial guess, and they must therefore be reasonable
!>     accurate.
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.0.3
!>  \date         1980s, 22 June 2021, 17 September 2021.
!>  \copyright    GNU Public License v2

subroutine atom_kb_test_dsolv2(iter, iconv,                              &
       nr, r, drdi, d2rodr,                                              &
       norb, no, lo, iso, zo,                                            &
       cdv, vhxc, ev, ek, ep,                                            &
       vlocal, vkbproj, inorm, rpsi,                                     &
       iowrite, mxdnr, mxdorb, mxdl)

! converted to fortran 90, March 1st 2018
! cleanup and new interface, July 2019. JLM
! exit bug. 22 June 2021. JLM
! mxdnr, mxdl, 17 September 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdorb                           !<  dimension of the number of orbitals
  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points
  integer, intent(in)               ::  mxdl                             !<  dimension angular momentum

  integer, intent(in)               ::  iter                             !<  iteration number
  integer, intent(in)               ::  iconv                            !<  convergence flag (if iconv = 1, calculates Hartree energy)

  integer, intent(in)               ::  nr                               !<  number of radial points

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i
  real(REAL64), intent(in)          ::  d2rodr(mxdnr)                    !<  (d^2 r(i) / d i^2) / (d r / d i)

  integer, intent(in)               ::  norb                             !<  number of orbitals

  integer, intent(in)               ::  no(mxdorb)                       !<  principal quantum number n
  integer, intent(in)               ::  lo(mxdorb)                       !<  angular quantum number l
  integer, intent(in)               ::  iso(mxdorb)                      !<  2*spin or 2*(j-l)
  real(REAL64), intent(in)          ::  zo(mxdorb)                       !<  orbital occupation

  real(REAL64), intent(in)          ::  vkbproj(mxdnr,0:mxdl,-1:1)       !<  kb-projector
  real(REAL64), intent(in)          ::  vlocal(mxdnr)                    !<  local pseudopotential
  integer, intent(in)               ::  inorm(0:mxdl,-1:1)               !<  sign of denominator of KB operator

  real(REAL64), intent(in)          ::  vhxc(mxdnr,mxdorb)               !<  screening potential in Rydberg for each orbital

! output

  real(REAL64), intent(out)         ::  cdv(mxdnr,-1:1)                  !<  4*pi*r**2 * valence charge density

  real(REAL64), intent(out)         ::  ek(mxdorb)                       !<  orbital kinetic energy
  real(REAL64), intent(out)         ::  ep(mxdorb)                       !<  orbital potential energy

! input and output

  real(REAL64), intent(inout)       ::  ev(mxdorb)                       !<  orbital energy

  real(REAL64), intent(inout)       ::  rpsi(mxdnr,0:mxdl,-1:1)          !<  wavefunctions (r(i),l,2j-2l).

! allocatable work arrays

  real(REAL64), allocatable         ::  v(:)                             !  work array (potential)
  real(REAL64), allocatable         ::  ar(:), br(:)                     !  work array (wave-functions)
  real(REAL64), allocatable         ::  work(:)                          !  work array (projector integration)

! local variables

  integer        ::  llp
  integer        ::  iflag

  real(REAL64)   ::  denr

  real(REAL64)   ::  revi

  real(REAL64)   ::  prowav


  logical        ::  lsp, lnosp

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  SMEV = 0.0001_REAL64
  real(REAL64), parameter    ::  TOL = 1.0E-12_REAL64

! counters

  integer     ::  i, j, k, l


  allocate(v(mxdnr))
  allocate(ar(mxdnr), br(mxdnr))
  allocate(work(mxdnr))

  revi = 10.0

! Initialize arrays for charge density

  do i = 1,nr
    cdv(i,-1) = ZERO
    cdv(i, 0) = ZERO
    cdv(i, 1) = ZERO
  enddo

! Start the loop over orbitals.
! Note that spin zero is treated as down.

  do i = 1,norb

    if (no(i) > 0 .and. (zo(i) /= 0.0 .or. iconv /= 0)) then

      if (ev(i) >= 0.0) ev(i) = -SMEV

!     Set up the potential, set the wave functionc array to zero-ar.

!     find the projector operator factor prowav

      l  = lo(i)
      llp = l*(l+1)

!     assume previous wave-function is similar to final

      prowav = ZERO

      do k = 1,nr
        work(k) = drdi(k)*r(k)*rpsi(k,l,iso(i))*vkbproj(k,l,iso(i))
      enddo

      do k = 1,nr-4,4
          prowav = prowav +  7*(work(k  )+work(k+4))                     &
                          + 32*(work(k+1)+work(k+3))                     &
                          + 12* work(k+2)
      enddo
      prowav = 2*inorm(l,iso(i))*prowav / (45*ONE)

      do j = 1,nr
        ar(j) = ZERO
        br(j) = ZERO
      enddo

      do j = 2,nr
        v(j) = vlocal(j) + llp/(r(j)*r(j)) + vhxc(j,i) +                 &
                   prowav*vkbproj(j,l,iso(i))*r(j)/rpsi(j,l,iso(i))
      enddo

!     Call the integration routine.

      call atom_atm_difnrl(nr, r, drdi, d2rodr, v, ar, br,               &
          no(i), lo(i), ev(i), iflag, TOL,                               &
          iowrite, mxdnr)

      if(iflag /= 0) then
        write(iowrite,*)
        write(iowrite,'("  dsolv2: wave-function not found, iflag = ",   &
                  &  i3,"  iter = ",i3,"  iorb = ",i3)') iflag, iter, i
        write(iowrite,*)

        if(iflag > 4) stop

      endif

!     update wave-function for next try

      do k = 1,nr
        rpsi(k,l,iso(i)) = (ar(k) + 2*rpsi(k,l,iso(i))) / 3
      enddo

!     Add to the charge density.
        do  j = 1,nr
          denr = zo(i) * ar(j) * ar(j)
          cdv(j,iso(i)) = cdv(j,iso(i)) + denr
        enddo

!     Compute various quantitities if last iteration.

      if (iconv == 1) then

        call atom_kb_test_orban(i, ar, br, nr, r, drdi,                  &
            norb, no, lo, zo, iso ,v, ev, ek, ep,                        &
            iowrite, mxdnr, mxdorb, mxdl)

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
    write(6,*) '   Stopped in kb_test_dsolv2.  Inconsistent spin structure'
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
  deallocate(work)

  return

end subroutine atom_kb_test_dsolv2
