!>  Finds the (non)-relativistic wave function and eigenvalue
!>  using finite differences and matrix diagonalization.
!>  An initial guess for the eigenvalues need not be supplied.
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.0.4
!>  \date         1980s, 22 June 2021. 18 September 2021. 12 October 2021.
!>  \copyright    GNU Public License v2

subroutine atom_kb_test_dsolv1(nr, a, b, r, drdi, norb,                  &
     no, lo, iso, zo, cdv, vhxc, ev,                                     &
     vlocal, vkbproj, inorm, rpsi,                                       &
     iowrite, mxdnr, mxdorb, mxdl)

! converted to fortran 90, March 1st 2018
! cleanup and new interface, July 2019. JLM
! mxdnr, mxdl, z(6,nr) bug. 18 September 2021. JLM
! atom_atm_trieig.  12 October 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdorb                           !<  dimension of the number of orbitals
  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points
  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum

  integer, intent(in)               ::  nr                               !<  number of radial points

  real(REAL64), intent(in)          ::  a                                !<  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64), intent(in)          ::  b                                !<  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  integer, intent(in)               ::  norb                             !<  number of orbitals

  integer, intent(in)               ::  no(mxdorb)                       !<  principal quantum number n
  integer, intent(in)               ::  lo(mxdorb)                       !<  angular quantum number l
  integer, intent(in)               ::  iso(mxdorb)                      !<  2*spin or 2*(j-l)
  real(REAL64), intent(in)          ::  zo(mxdorb)                       !<  orbital occupation

  real(REAL64), intent(in)          ::  vhxc(mxdnr,-1:1)                 !<  screening potential (Ry)

  real(REAL64), intent(in)          ::  vkbproj(mxdnr,0:mxdl,-1:1)       !<  kb-projector
  real(REAL64), intent(in)          ::  vlocal(mxdnr)                    !<  local pseudopotential
  integer, intent(in)               ::  inorm(0:mxdl,-1:1)               !<  sign of denominator of KB operator

! output

  real(REAL64), intent(out)         ::  ev(mxdorb)                       !<  orbital energy


  real(REAL64), intent(out)         ::  cdv(mxdnr,-1:1)                  !<  4*pi*r**2 * valence charge density

  real(REAL64), intent(out)         ::  rpsi(mxdnr,0:mxdl,-1:1)          !<  wavefunctions (r(i),l,2j-2l).


! allocatable work arrays

  real(REAL64), allocatable         ::  dk(:), d(:)
  real(REAL64), allocatable         ::  sd(:)

  integer, allocatable              ::  nmax(:,:)

  real(REAL64), allocatable         ::  z(:,:)
  real(REAL64), allocatable         ::  e(:)

! local variables

  integer                           ::  lmax                             !  maximum angular momentum l

  real(REAL64)   ::  c1, c2

  integer        ::  nmaxmax

  integer        ::  nrm
  integer        ::  llp

  integer        ::  ki

  real(REAL64)   ::  denr

  real(REAL64)   ::  prowav

  logical        ::  lsp, lnosp

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  PONE = 0.1_REAL64
  real(REAL64), parameter    ::  OPF = 1.5_REAL64
  real(REAL64), parameter    ::  SMALL = 1.0E-6_REAL64

! counters

  integer     ::  i, j, k, l


  allocate(nmax(0:mxdl,-1:1))

  allocate(dk(mxdnr), d(mxdnr))
  allocate(sd(mxdnr))

  do i = 1,nr
    cdv(i,-1) = ZERO
    cdv(i, 0) = ZERO
    cdv(i, 1) = ZERO
  enddo

! Find the max n given l and s (or j-l).

  do l = 0,mxdl
    do i = -1,1
      nmax(l,i) = 0
    enddo
  enddo

  lmax = 0
  nmaxmax = 0
  do k = 1,norb
    l = lo(k)
    if (no(k) > 0) then
      nmax(l, iso(k)) = max(no(k),nmax(l, iso(k)))
      nmaxmax = max(nmax(l, iso(k)),nmaxmax)
    endif
    if(l > lmax) lmax = l
  enddo

  allocate(z(mxdnr,nmaxmax))
  allocate(e(nmaxmax))


! Set up hamiltonian matrix for kinetic energy.
! Only the diagonal depends on the potential.

  c2 = -ONE/b**2
  c1 = -2*ONE*c2 + ONE/4
  dk(1)  = c1 / (r(2)+a)**2
  sd(1)  = ZERO
  do i = 3,nr
    dk(i-1)  = c1 / (r(i)+a)**2
    sd(i-1)  = c2 / ((r(i)+a)*(r(i-1)+a))
  enddo

! Start loop over spins

  nrm = nr - 1

  do i = -1,1

! Start loop over s p d... states.

    do l = 0,lmax

      if (nmax(l,i) /= 0) then


        llp = l*(l+1)

!       assumes wave-function is similar to reference case

        prowav = ZERO

        do k = 1,nr-4,4
          prowav = prowav +  7*(vkbproj(k  ,l,i)*drdi(k)*r(k)*r(k) +          &
                                vkbproj(k+4,l,i)*drdi(k+4)*r(k+4)*r(k+4))     &
                          + 32*(vkbproj(k+1,l,i)*drdi(k+1)*r(k+1)*r(k+1) +    &
                                vkbproj(k+3,l,i)*drdi(k+3)*r(k+3)*r(k+3))     &
                          + 12* vkbproj(k+2,l,i)*drdi(k+2)*r(k+2)*r(k+2)
        enddo
        prowav = 2*inorm(l,i)*prowav / (45*ONE)

        do k = 2,nr
            d(k-1) = dk(k-1) + vlocal(k) + llp/(r(k)*r(k)) + vhxc(k,i)     &
                         + prowav*vkbproj(k,l,i)
        enddo

!       Diagonalize the matrix.

        call atom_atm_trieig(nrm, d, sd, nmax(l,i), e, z, iowrite, mxdnr)

!       Save the energy levels and add to charge density.

        ki = 1

        do k = 1,norb

          if (no(k) > 0 .and. lo(k) == l .and. iso(k) == i) then

            ev(k) = e(ki)

            rpsi(1,l,i) = ZERO
            do j = 2,nr
              rpsi(j,l,i) = z(j-1,ki)
              denr = zo(k) * z(j-1,ki)*z(j-1,ki) / drdi(j)
              cdv(j,i) = cdv(j,i) + denr
            enddo

            ki = ki + 1

          endif

        enddo

      endif

    enddo
  enddo

! End loop over s p and d states.

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
    write(6,*) '   Stopped in kb_test_dsolv1.  Inconsistent spin structure'
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

  deallocate(nmax)
  deallocate(dk, d)
  deallocate(sd)
  deallocate(z, e)

  return

end subroutine atom_kb_test_dsolv1
