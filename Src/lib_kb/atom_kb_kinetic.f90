!>  Calculates the kinetic energies of orbitals and suggests cutoff
!>
!>  \author       J.L.Martins
!>  \version      6.0.5
!>  \date         28 September 2021.
!>  \copyright    GNU Public License v2

subroutine atom_kb_kinetic(npot, lo, zo, irel, nr, r, drdi,              &
     nqbas, delqbas, rpsi, ektot,                                        &
     mxdl, mxdnr, nqmax)

! adapted from atom_kb_pot_four 28 September 2021. JLM
! warning for kinetic. 22 October 2021. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum
  integer, intent(in)               ::  mxdnr                            !<  dimension of radial grid points.
  integer, intent(in)               ::  nqmax                            !<  dimension of fourier grid points for potential

  integer, intent(in)               ::  npot(-1:1)                       !<  number of orbitals to be calculated
  integer, intent(in)               ::  lo(mxdl+1,-1:1)                  !<  angular momentum of orbital
  real(REAL64), intent(in)          ::  zo(mxdl+1)                       !<  orbital occupation (negative if abesent)

  character(len=3), intent(in)      ::  irel                             !<  flag for relativistic (r) or spin-polarized (s) original calculation

  integer, intent(in)               ::  nr                               !<  number of points in the radial grid
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points   r(i) = a*(exp(b*i)-1), i=1,...,nr
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  integer, intent(in)               ::  nqbas                            !<  number of points in the fourier grid for basis
  real(REAL64), intent(in)          ::  delqbas                          !<  spacing of the fourier grid for basis

  real(REAL64), intent(in)          ::  rpsi(mxdnr,0:mxdl,-1:1)          !<  wavefunctions (r(i),l,2j-2l).

! output

  real(REAL64), intent(out)         ::  ektot(0:nqmax)                   !<  Cumulative kineti energy of wave-functions in Fourier space


! allocatable arrays

  real(REAL64), allocatable         ::  fin(:)
  real(REAL64), allocatable         ::  drpsidr(:)

  real(REAL64), allocatable         ::  fout(:), yp(:), ypp(:), w(:,:)
  real(REAL64), allocatable         ::  ypwf(:), fkin(:), ans(:)

  real(REAL64), allocatable         ::  ek(:,:)

! local variables

  integer                      ::  l, jmax

  real(REAL64)                 ::  xk(-3:3)
  real(REAL64)                 ::  fk(-3:3)
  real(REAL64)                 ::  y(0:1), dy(0:1)
  integer                      ::  kref

  integer                      ::  ierr

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  real(REAL64), parameter    ::  SMALL = 2.0E-4_REAL64

! counters

  integer                  :: i, j, k, n

  allocate(fin(mxdnr))
  allocate(drpsidr(mxdnr))
  allocate(fout(0:nqmax), yp(0:nqmax), ypp(0:nqmax))
  allocate(ypwf(0:nqmax))
  allocate(fkin(0:nqmax))
  allocate(ans(0:nqmax))
  allocate(w(0:nqmax,3))
  allocate(ek(0:mxdl,-1:1))

  jmax = 0
  if(irel == 'rel') jmax = 1

  do k = 0,nqbas
    ypwf(k) = k*delqbas
  enddo

  do k = 0,nqbas
    ektot(k) = ZERO
  enddo

! Do loop over wave-functions

  do j = -jmax,jmax

    do i = 1,npot(j)
      l = lo(i,j)

      do k = 1,nr
        kref = k
        if(kref < 4) kref = 4
        if(kref > nr - 3) kref = nr - 3
        do n = -3,3
          xk(n) = r(kref+n) - r(k)
          fk(n) = rpsi(kref+n,l,j)
        enddo
        call poly_interp(y, dy, xk, fk, 6, 1)
        drpsidr(k) = y(1)
      enddo

      do k = 2,nr
        fin(k) = drpsidr(k)*drpsidr(k) +                                 &
                 rpsi(k,l,j)*rpsi(k,l,j)*(l*(l+1)) / (r(k)*r(k))
      enddo
      if(l == 0) then
        fin(1) = drpsidr(1)*drpsidr(1)
      else
        fin(1) = ZERO
      endif

      call  atom_atm_integ(ek(l,j), fin, drdi, nr)

    enddo

  enddo



! fourier transforms of wave functions

  do j = -jmax,jmax

    do i = 1,npot(j)
      l = lo(i,j)

      do k = 1,nr
        fin(k) = rpsi(k,l,j)
      enddo

      call ft_radial(nr-1, nqbas, r, drdi, ypwf, fin, fout, l)

      do k = 0,nqbas
        fkin(k) = fout(k)*fout(k)*ypwf(k)**4
      enddo

      call splift (ypwf,fkin,yp,ypp,nqbas+1,w,ierr,0,-0.5d0,-0.5d0,0.0d0,0.0d0)

      call spliq(ypwf,fkin,yp,ypp,nqbas+1,0.0d0,ypwf,nqbas+1,ans,ierr)

      if(abs(ans(nqbas)-ek(l,j)) / (ONE + abs(ans(nqbas))) > SMALL) then
        write(6,*)
        write(6,*) '    WARNING  in atom_kb_kinetic'
        write(6,*) '    inconsistent kinetic energies l = ',l,' j= ',j
        write(6,*) '    E-kinetic = ', ek(l,j), ans(nqbas)
        write(6,*)
      endif

      if(j == 0) then
        do k = 0,nqbas
          ektot(k) = ektot(k) + zo(l+1)*ans(k)
        enddo
      elseif(j == 1) then
        do k = 0,nqbas
          ektot(k) = ektot(k) + (l+1)*zo(l+1)*ans(k) / ((2*l+1)*ONE)
        enddo
      else
        do k = 0,nqbas
          ektot(k) = ektot(k) + l*zo(l+1)*ans(k) / ((2*l+1)*ONE)
        enddo
      endif

    enddo

  enddo

  deallocate(fin)
  deallocate(drpsidr)
  deallocate(fout, yp, ypp)
  deallocate(ypwf)
  deallocate(fkin)
  deallocate(ans)
  deallocate(w)
  deallocate(ek)

  return

end subroutine atom_kb_kinetic
