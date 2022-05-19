!>  Transforms to Fourier space all the radial quantities
!>
!>  \author       Norm Troullier, J.L.Martins
!>  \version      6.0.8
!>  \date         early 90s, May 2012, 1 September 2021, 19 May 2022.
!>  \copyright    GNU Public License v2

subroutine atom_kb_pot_four(npot, lo, irel, nicore, nr, r, drdi, zion,   &
     nql, nqnl, delql, nqbas, delqbas, vlocal, inorm, vkbproj,           &
     cdc, cdv, vlocft, vql0, vkbprft, cdcft, cdvft, basft,               &
     norbas, lo_b, rpsi,                                                 &
     mxdl, mxdnr, nqmax, nqmaxbas, mxdbas)

! adapted from the old program jlm 11/6/2012
! mxdl, etc... 18 September 2021. JLM
! nrm -> nrm-1 bug. 20 October 2021. JLM
! deallocate, 19 May 2022. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum
  integer, intent(in)               ::  mxdnr                            !<  dimension of radial grid points.
  integer, intent(in)               ::  nqmax                            !<  dimension of fourier grid points for potential
  integer, intent(in)               ::  nqmaxbas                         !<  dimension of fourier grid points for basis
  integer, intent(in)               ::  mxdbas                           !<  dimension of basis functions

  integer, intent(in)               ::  npot(-1:1)                       !<  number of orbitals to be calculated
  integer, intent(in)               ::  lo(mxdl+1,-1:1)                  !<  angular momentum of orbital
  character(len=3), intent(in)      ::  irel                             !<  flag for relativistic (r) or spin-polarized (s) original calculation
  character(len=4), intent(in)      ::  nicore                           !<  flag for core correction

  integer, intent(in)               ::  nr                               !<  number of points in the radial grid
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points   r(i) = a*(exp(b*i)-1), i=1,...,nr
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  real(REAL64), intent(in)          ::  zion                             !<  ionic charge

  integer, intent(in)               ::  nql                              !<  number of points for the local Fourier grid
  integer, intent(in)               ::  nqnl                             !<  number of points for the non-local Fourier grid
  real(REAL64), intent(in)          ::  delql                            !<  spacing of points in Fourier grid
  integer, intent(in)               ::  nqbas                            !<  number of points in the fourier grid for basis
  real(REAL64), intent(in)          ::  delqbas                          !<  spacing of the fourier grid for basis

  real(REAL64), intent(in)          ::  vlocal(mxdnr)                    !<  r*local pseudopotential
  integer, intent(in)               ::  inorm(0:mxdl,-1:1)               !<  sign of denominator of KB operator
  real(REAL64), intent(in)          ::  vkbproj(mxdnr,0:mxdl,-1:1)       !<  kb-projector
  real(REAL64), intent(in)          ::  cdc(mxdnr)                       !<  4*pi*r**2 charge density of core
  real(REAL64), intent(in)          ::  cdv(mxdnr)                       !<  4*pi*r**2 charge density of valence.  kb_conv is not spin polarized...

  integer, intent(in)               ::  norbas                           !<  number of basis functions
  integer, intent(in)               ::  lo_b(mxdbas)                     !<  angular momentum of basis function
  real(REAL64), intent(in)          ::  rpsi(mxdnr,mxdbas)               !<  wavefunctions (r(i),l,2j-2l).

! output

  real(REAL64), intent(out)         ::  vlocft(0:nqmax)                  !<  Fourier transform of local pseudopotential
  real(REAL64), intent(out)         ::  vql0                             !<  zero frequency of the integral without the Coulomb part
  real(REAL64), intent(out)         ::  vkbprft(0:nqmax,0:mxdl,-1:1)     !<  Fourier transform of kb-projector
  real(REAL64), intent(out)         ::  cdcft (0:nqmax)                  !<  Fourier transform of 4*pi*r**2 charge density of core
  real(REAL64), intent(out)         ::  cdvft(0:nqmax)                   !<  Fourier transform of 4*pi*r**2 charge density of valence.  kb_conv is not spin polarized...
  real(REAL64), intent(out)         ::  basft(0:nqmaxbas,mxdbas)         !< Fourier (Bessel) transform of the radial wavefunction

! allocatable arrays

  real(REAL64), allocatable         ::  fin(:)
  real(REAL64), allocatable         ::  fout(:),yp(:)
  real(REAL64), allocatable         ::  ypwf(:)

! local variables

  real(REAL64)                 ::  crnorm
  integer                      ::  l, nrm, jmax

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  TINYR = 1.0E-12_REAL64
  real(REAL64), parameter    ::  PI4 = 16*atan(1.0_REAL64)

! counters

  integer                  :: i, j, k

  allocate(fin(mxdnr))
  allocate(fout(0:nqmax),yp(0:nqmax))
  allocate(ypwf(0:nqmaxbas))

  jmax = 0
  if(irel == 'rel') jmax = 1

  do k = 0,nql
    yp(k) = k*delql
  enddo

! Fourier transform local potential,

! ft begins at 0 to have explicit r(1) = 0

  vql0 = ZERO
  vlocft = ZERO
  call ft_local(nr-1, nql, r, drdi, yp, zion, vlocal, vql0, vlocft)

! Do loop over K&B projector potentials

  vkbprft = ZERO
  do j = -jmax,jmax

    do i = 1,npot(j)
      l = lo(i,j)
      if(inorm(l,j) /= 0) then
        do k = nr,1,-1
          if (abs(vkbproj(k,l,j)) < TINYR) then
            nrm = k
          else
            exit
          endif
        enddo

        do k = 1,nrm
          fin(k) = vkbproj(k,l,j)
        enddo

        call ft_precise(nrm-1, nqnl, r, drdi, yp, fin, fout, l)

        do k = 0,nqnl
          vkbprft(k,l,j) = fout(k)
        enddo

      endif
    enddo

  enddo

! Fourier transform core charge density.
! if it is =/= 0.

  cdcft = ZERO
  if (nicore /= 'nc  ') then

    call ft_charge(nr-1, nql, r, drdi, yp, cdc, cdcft)

  endif

! Fourier transform the valence charge density.

  cdvft = ZERO
  call ft_charge(nr-1, nql, r, drdi, yp, cdv, cdvft)

  cdvft(0) = zion


! fourier transforms of wave functions

  do j = 0,nqbas
    ypwf(j) = j*delqbas
  enddo

  basft = ZERO
  do i = 1,norbas
    l = lo_b(i)
    do k = 1,nr
      fin(k) = rpsi(k,i)
    enddo

    call ft_radial(nr-1, nqbas, r, drdi, ypwf, fin, fout, l)

    crnorm = sqrt((2*l+1)/PI4)
    do k = 0,nqbas
      basft(k,i) = fout(k)*crnorm
    enddo


  enddo

  deallocate(fin)
  deallocate(fout,yp)
  deallocate(ypwf)

  return

end subroutine atom_kb_pot_four
