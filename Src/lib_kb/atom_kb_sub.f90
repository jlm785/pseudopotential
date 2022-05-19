!>  Chooses the local potential in KB transformation
!>
!>  \author       Norm Troullier, J.L.Martins
!>  \version      6.0.8
!>  \date         November 90, May 2012, July 2021, 19 May 2022.
!>  \copyright    GNU Public License v2

subroutine atom_kb_sub(llocal, nql, delql, nqbas, delqbas,               &
         lmax_bas, n_bas, r_bas, nz_bas,                                 &
         iowrite, ioparsec, fileparsec, iokb, sfilekb, ioupf, sfileupf,  &
         iopsdkb, filepsdkb, ioplot, fileplot, mxdnr, mxdl)

! converts a semi-local pseudopotential to the Kleinman and Bylander
! form, PRL 48, 1425 (1982) and writes files in several formats with the
! KB operators in both real and Fourier space.
!
!
! The program uses          RYDBERG UNITS
!
!
! the first version was written by Norm Troullier while at the U. of MN                 *
! and was Copyright Norm Troullier and Jose Luis Martins
! and  Version 1.25 was Dated Nov. 12, 1990
! Minor modifications by jlm, 2/2/96, 13/6/97, 11/2001, 10/2002, 1/2008
!
! This version was written in May 2012, breaking the long program in subroutines
! and adding the spin-orbit.
! Transform from program to subroutine, July 2021. JLM
! evl, vsiesta. 28 September 2021. JLM
! kinetic, 4 October 2021. JLM
! ev for scattering basis functions. 19 October 2021. JLM
! printing. 21 October 2021. JLM
! psdtitle. 19 May 2022. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points

  integer, intent(in)               ::  mxdl                             !<  dimension for angular momentum

  integer, intent(in)               ::  iowrite                          !<  default output tape

  integer, intent(in)               ::  ioparsec                         !<  default tape for pseudopotential in parsec format
  character(len=*), intent(in)      ::  fileparsec                       !<  name of default tape for reading pseudopotential in parsec format

  integer, intent(in)               ::  iokb                             !<  default tape for pseudopotential in KB format
  character(len=*), intent(in)      ::  sfilekb                          !<  suffix for default tape for writing pseudopotential in KB format

  integer, intent(in)               ::  ioupf                            !<  default tape for pseudopotential in UPF format
  character(len=*), intent(in)      ::  sfileupf                         !<  suffix for default tape for writing pseudopotential in UPF format

  integer , intent(in)              ::  iopsdkb                          !<  default tape for KB pseudopotential in real space
  character(len=15), intent(in)     ::  filepsdkb                        !<  name of default tape for writeing KB pseudopotential in real space

  integer, intent(in)               ::  ioplot                           !<  tape for later plotting
  character(len=*), intent(in)      ::  fileplot                         !<  name of file with plot information


  integer, intent(in)               ::  llocal                           !<  angular momentum for local potential (negative: maximum of l-dependent)
  integer, intent(in)               ::  nql                              !<  number of points for the Fourier grid of local potential and densities
  real(REAL64), intent(in)          ::  delql                            !<  spacing of points in Fourier grid

  integer, intent(in)               ::  nqbas                            !<  number of points for the Fourier grid of basis
  real(REAL64), intent(in)          ::  delqbas                          !<  spacing of points in Fourier grid for basis set

  integer,intent(in)                ::  lmax_bas                         !<  maximum angular momentum in basis
  integer,intent(in)                ::  n_bas(0:mxdl)                    !<  basis functions for angular momentum l
  real(REAL64), intent(in)          ::  r_bas(3,0:mxdl)                  !<  cutoff for the basis (up to triple zeta and l = 4)
  integer, intent(in)               ::  nz_bas(3,0:mxdl)                 !<  number of non-trivial zeroes in basis function

! Size of arrays

  integer                           ::  lmax_pot                         !  maximum number of angular momentum components
  integer                           ::  nrmax                            !  dimension of radial grid points.
  integer                           ::  nqmax                            !  maximum number of fourier grid points for potential

! general variables

  character(len=10)                 ::  dated                            !  date of calculation
  character(len=5)                  ::  version                          !  program version (should be the same as in atomic program)

  character(len=30)                 ::  filekb                           !  name of default tape for writing pseudopotential in KB format
  character(len=30)                 ::  fileupf                           !  name of default tape for writing pseudopotential in upf format

! variables from the output file of psd_gen

  character(len=2)                  ::  nameat                           !  chemical symbol of the element
  character(len=2)                  ::  icorr                            !  correlation used in the calculation
                                                                         !      ca -> Ceperley-Alder LDA parametrized by Perdew-Zunger
                                                                         !      pw -> Ceperley-Alder LDA parametrized by Perdew-Wang 92
                                                                         !      pb -> The pbe GGA
  character(len=3)                  ::  irel                             !  flag for relativistic (r) or spin-polarized (s) original calculation
                                                                         !      nrl -> non-relativistic and non-spin-polarized
                                                                         !      rel -> relativistic (includes spin-orbid)
                                                                         !      isp -> spin-polarized non-relativistic (no spin-orbit)
  character(len=4)                  ::  nicore                           !  flag for core correction
                                                                         !      nc   -> no core correction
                                                                         !      pcec -> partial core correction
                                                                         !      fcec -> full core correction
  character(len=10)                 ::  irdate,irvers                    !  date and version of original calculation
  character(len=10)                 ::  irayps(4)                        !  type of pseudopotential
  character(len=10)                 ::  psdtitle(20)                     !  pseudopotential parameters
  integer                           ::  npot(-1:1)                       !  number of orbitals (s,p,d,...). -1:  j=l-1/2.  0:  average.  1:  j=l+1/2
                                                                         !      meaning depends on irel
  integer                           ::  nr                               !  number of points in the radial grid
  real(REAL64)                      ::  a,b                              !  constants used to generate the radial grid
  real(REAL64), allocatable         ::  r(:)                             !  radial grid points   r(i) = a*(exp(b*i)-1), i=1,...,nr
  real(REAL64)                      ::  zion                             !  ionic charge
  integer, allocatable              ::  lo(:,:)                          !  angular momentum of orbital
  real(REAL64), allocatable         ::  vionic(:,:,:)                    !  r*pseudopotential. -1:  j=l-1/2.  0:  average.  1:  j=l+1/2
  real(REAL64), allocatable         ::  cdc(:)                           !  4*pi*r**2 charge density of core
  real(REAL64), allocatable         ::  cdv(:)                           !  4*pi*r**2 charge density of valence.  kb_conv is not spin polarized...

! variables for KB-pseudopotential

  integer                           ::  nqnl                             !  number of points for the Fourier grid of non-local potential
  real(REAL64), allocatable         ::  drdi(:)                          !  d r(i) / d i
  real(REAL64), allocatable         ::  d2rodr(:)                        !  (d^2 r / d i^2) / (d r / d i)

  integer                           ::  ifcore                            !  derived from nicore
  real(REAL64)                      ::  totvel                           !  total number of valence electrons
  real(REAL64), allocatable         ::  vscreen(:)                       !  screening potential
  real(REAL64), allocatable         ::  vlocal(:)                        !  r*local pseudopotential
  integer, allocatable              ::  inorm(:,:)                       !  sign of denominator of KB operator
  real(REAL64), allocatable         ::  vkbproj(:,:,:)                   !  kb-projector(r(i),l,2j-2l)

! variables for wavefunctions

  real(REAL64), allocatable         ::  ev(:,:)                          !  eigenvalues (l,2j-2l).
  real(REAL64), allocatable         ::  rpsi(:,:,:)                      !  wavefunctions (r(i),l,2j-2l).

  real(REAL64), allocatable         ::  zo(:)                            !  orbital occupation (negative if abesent)
  real(REAL64), allocatable         ::  rc(:)                            !  core radius

  integer                           ::  norbas                           !  number of basis functions

  real(REAL64), allocatable         ::  rpsi_b(:,:)                      !  basis wavefunctions (r(i),n).
  real(REAL64), allocatable         ::  drpsidr_b(:)                     !  drpsi_b / dr
  real(REAL64), allocatable         ::  veff_b(:)                        !  effective potential for basis
  integer, allocatable              ::  lo_b(:)                          !  angular momentum of basis function
  real(REAL64), allocatable         ::  ev_b(:)                          !  variational energy of the orbital
  integer, allocatable              ::  nrc_b(:)                         !  maximum radius of basis function

  integer                           ::  iflag

! variables for Fourier transform of wavefunctions and potentials

  real(REAL64), allocatable         ::  basft(:,:)                        !  Fourier (Bessel) transform of the radial wavefunction
  real(REAL64), allocatable         ::  vlocft(:)                         !  Fourier transform of local pseudopotential
  real(REAL64)                      ::  vql0                              !  zero frequency of the integral without the Coulomb part
  real(REAL64), allocatable         ::  vkbprft(:,:,:)                    !  Fourier transform of kb-projector
  real(REAL64), allocatable         ::  cdcft(:)                          !  Fourier transform of 4*pi*r**2 charge density of core
  real(REAL64), allocatable         ::  cdvft(:)                          !  Fourier transform of 4*pi*r**2 charge density of valence.  kb_conv is not spin polarized...

  real(REAL64), allocatable         ::  ektot(:)                          !  Cumulative kineti energy of wave-functions in Fourier space


! other local variables

  real(REAL64)                      ::  rtry
  real(REAL64)                      ::  evl
  integer                           ::  nn

  real(REAL64)                      ::  revi
  real(REAL64)                      ::  vsiesta

  real(REAL64)                      ::  qb, qa

  character(len=5)                  ::  scorr

! constants

  real(REAL64), parameter           ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  real(REAL64), parameter           ::  SMALL = 1.0E-6_REAL64
  real(REAL64), parameter           ::  TOL = 1.0E-10_REAL64

! counters

  integer                           ::  i, j, n, l


  write(6,*)
  write(6,*) '  Constructing KB operator and basis functions'
  write(6,*)

  call zedate(dated)
  call atom_version(version)

  write(iowrite,*)
  write(iowrite,*)
  write(iowrite,'("   Kleinman-Bylander pseudopotential conversion")')
  write(iowrite,'("   and Bessel/Fourier transform program. Version  ",a5)') version
  write(iowrite,*)
  write(iowrite,'("   Run on ",a10)') dated
  write(iowrite,*)

  call atom_kb_psd_in_parsec_size(ioparsec, fileparsec, lmax_pot, nrmax)

  write(iowrite,'("   Size of arrays:   lmax = ",i5,"    mxdnr = ",i8)')  &
               lmax_pot, mxdnr

! check dimensions and allocate the arrays (paranoid check)

  if(nrmax < 1 .or. lmax_pot < 0) then
    write(6,*)
    write(6,'("   NEGATIVE DIMENSIONS,  mxdnr = ",i7,"  LMAX = ",i7)')   &
                 mxdnr, lmax_pot
    write(6,*)

    STOP

  endif

  if(lmax_pot > mxdl) then
    write(6,*)
    write(6,'("   lmax = ",i7," is greater than mxdl = ",i7)')           &
                 lmax_pot, mxdl
    write(6,*)

    STOP

  endif

  if(nrmax > mxdnr) then
    write(6,*)
    write(6,'("   nrmax = ",i7," is greater than mxdnr = ",i7)')           &
                 nrmax, mxdnr
    write(6,*)

    STOP

  endif

  allocate(r(mxdnr))
  allocate(lo(mxdl+1,-1:1))
  allocate(vionic(mxdnr,0:mxdl,-1:1))
  allocate(cdc(mxdnr),cdv(mxdnr))

  allocate(zo(0:mxdl))
  allocate(rc(0:mxdl))

  allocate(ev(0:mxdl,-1:1))
  allocate(rpsi(mxdnr,0:mxdl,-1:1))

  r = ZERO
  lo = 0
  vionic = ZERO
  cdc = ZERO
  cdv = ZERO

! reads the pseudopotential data file


  call atom_kb_psd_in_parsec(ioparsec, fileparsec,                       &
         nameat, icorr, irel, nicore, irdate, irvers, irayps, psdtitle,  &
         npot, nr, a, b, r, zion, lo, vionic, cdc, cdv,                  &
         zo, rc, rpsi(:,:,0),                                            &
         lmax_pot, mxdnr, mxdl)


  write(iowrite,*)
  write(iowrite,'("   The semi-local pseudopotential for ",a2,           &
       &  " was created on ",a10)') nameat, irdate
  write(iowrite,'("   with atomic program version ",a10)') irvers
  write(iowrite,'("   Pseudopotential type: ",4a10)') (irayps(i),i=1,4)
  write(iowrite,*)
  write(iowrite,'("   Pseudo parameters: ")')
  write(iowrite,'(3x,20a10)') (psdtitle(i),i=1,20)
  write(iowrite,*)


  nqmax = nql

  nqnl = nql

  write(iowrite,*)
  write(iowrite,'("   The Fourier grid will range from q = 0.0000 to",  &
       &  f10.4, "  with steps of ",f10.4)') nql*delql, delql

! some setup

  write(irayps(4),'("kb-loc=",i2,1X)') llocal

  if (nicore == 'nc  ') then
    ifcore = 0
  elseif(nicore == 'pcec' .or. nicore == 'fcec') then
    ifcore = 1
  else
    ifcore = 2
  endif

  allocate(drdi(mxdnr))
  allocate(d2rodr(mxdnr))

! paranoid check

  do i = 1,nr
    rtry = a*(exp(b*(i-1))-ONE)

    if(abs(rtry - r(i)) > SMALL) then
      write(iowrite,'("   INCOMPATIBLE  a, b r(i) ")')

      STOP

    endif

    drdi(i) = (r(i)+a)*b
    d2rodr(i) = b
  enddo

  allocate(vscreen(mxdnr))

  call atom_kb_screen(nr, r, drdi, cdv, cdc, icorr, ifcore, totvel, vscreen,   &
         iowrite, mxdnr)

  write(iowrite,*)
  write(iowrite,'("   Total valence electrons = ",f14.6)') totvel

  if(abs(totvel - zion) > SMALL) then
    write(iowrite,*)
    write(iowrite,'("   Ionic charge:",f14.6)') zion - totvel
  endif

! make ev zero for scattering basis functions

  do i = -1,1
  do l = 1,mxdl
    ev(l,i) = ZERO
  enddo
  enddo


  call atom_kb_wvfct(npot, lo, irel, nr, r, drdi, d2rodr,                &
         vionic, vscreen, ev, rpsi,                                      &
         iowrite, mxdl, mxdnr)

  write(iowrite,'(/,"   Orbital energies",//,5x,"n",5x,"l",13x,"e",/)')
  do i=1,npot(0)
    write(iowrite,'(i6,i6,4x,f14.6)') lo(i,0)+1,lo(i,0),ev(lo(i,0),0)
  enddo

  if(irel == 'rel') then
    write(iowrite,'(/,5x,"n",5x,"l",5x,"j",15x,"e",/)')
    do i=1,npot(1)
      write(iowrite,'(i6,i6,i6,"/2",4x,f14.6)')                          &
               lo(i,1)+1,lo(i,1),2*lo(i,1)+1,ev(lo(i,1),1)
    enddo
    do i=1,npot(-1)
      write(iowrite,'(i6,i6,i6,"/2",4x,f14.6)')                          &
               lo(i,-1)+1,lo(i,-1),2*lo(i,-1)-1,ev(lo(i,-1),-1)
    enddo

  endif

  allocate(ektot(0:nqmax))


  write(iowrite,*)
  write(iowrite,'("  Calculating Bessel transforms may take a few seconds ")')
  write(iowrite,*)
  write(iowrite,*) '  Suggested range of kinetic energy cutoff'
  write(iowrite,*) '  for plane-wave calculations:'
  write(iowrite,*)

  call atom_kb_kinetic(npot, lo, zo, irel, nr, r, drdi,                  &
     nqbas, delqbas, rpsi, ektot,                                        &
     mxdl, mxdnr, nqmax)

  qa = nqbas*delqbas
  qb =  nqbas*delqbas
  do j = nqbas,0,-1
    if(abs(ektot(nqbas)-ektot(j)) < 0.0004) qa = j*delqbas
    qb = j*delqbas
    if(abs(ektot(nqbas)-ektot(j)) > 0.005) exit
  enddo

  write(iowrite,*) '  lower  cutoff (bands)   ',nint(qb*qb/2),' Ha'
  write(iowrite,*) '  higher cutoff (pressure)',nint(qa*qa/2),' Ha'
  write(iowrite,*)

  allocate(vlocal(mxdnr))

  call atom_kb_choose_local(npot, lo, llocal, nr, r, vionic, vlocal,     &
      mxdl, mxdnr)


  allocate(inorm(0:mxdl,-1:1))
  allocate(vkbproj(mxdnr,0:mxdl,-1:1))
  vkbproj = ZERO


  call atom_kb_projector(npot, lo, nr, r, drdi, irel, llocal,            &
      vlocal, vionic, rpsi, inorm, vkbproj,                              &
      iowrite, mxdl, mxdnr)


!   Find the 1st and 2nd eigenvalues for all angular
!   momentum using just the local potential.  Print
!   out results, this is used to find if any ghost
!   states exist for the potential.  See (preprint)paper
!   of Gonze, Kackell, and Scheffler, Phy. Rev. B.


  call atom_kb_ghost_test(npot, lo, llocal, irel, nr, r, drdi, d2rodr,   &
     vscreen, vlocal, ev, inorm, irayps,                                 &
     iowrite, mxdl, mxdnr)

! writes a file with the results in real space

  call atom_kb_psd_out_real(iopsdkb, filepsdkb,                          &
      nameat, icorr, irel, nicore, irdate, irvers, irayps, psdtitle,     &
      npot, lo, nr, a, b, r, zion,                                       &
      vlocal, inorm, vkbproj,                                            &
      cdc, cdv,                                                          &
      mxdl, mxdnr)

! finds the confined wavefunctions

  norbas = 0
  do l = 0,lmax_bas
    norbas = norbas + n_bas(l)
  enddo

  allocate(rpsi_b(mxdnr,norbas))
  allocate(drpsidr_b(mxdnr))
  allocate(veff_b(mxdnr))
  allocate(lo_b(norbas))
  allocate(ev_b(norbas))
  allocate(nrc_b(norbas))

  n = 0
  do l = 0,lmax_bas
  do i = 1,n_bas(l)
    n = n+1
    lo_b(n) = l

    evl =  ev(l,0)
    nn = l+1 + nz_bas(i,l)
    revi = r_bas(i,l)

!     call atom_kb_basis(nn, l, nr, r, drdi, d2rodr, r_bas(i,l), nrc_b(n), &
!         vionic(:,:,0), vscreen, evl, iflag, rpsi_b(:,n),                 &
!         iowrite, mxdl, mxdnr)

!   SIESTA recipe (only used up to revi)

    do j = 2,nr
      veff_b(j) = (vionic(j,l,0)  + (l*(l+1))/r(j)) / r(j)  + vscreen(j)
      if(r(j) > 0.9*revi .and. r(j) < revi) then
        vsiesta = 40*( exp(-0.1*revi/(r(j)-0.9*revi)) / (revi- r(j) + SMALL) )
        veff_b(j) = veff_b(j) + vsiesta
      endif
!       if(r(j) > 0.8*revi .and. r(j) < revi) then
!         vsiesta = 2.0*(r(j) - 0.8*revi)**3
!         veff_b(j) = veff_b(j) + vsiesta
!       endif
    enddo

    call atom_atm_difnrl_zr(nr, r, drdi, d2rodr,                         &
        veff_b, rpsi_b(:,n), drpsidr_b,                                  &
        nn, l, evl, iflag, r_bas(i,l), nrc_b(n), TOL,                     &
        iowrite, mxdnr)

    if(iflag /= 0) then
      write(6,*)
      write(6,*) '  Stopped in atom_kb_sub:'
      write(6,*) '  Could not get basis functions, iflag = ',iflag
      write(6,*)

      STOP

    endif

  enddo
  enddo


! calculates the Fourier transforms


  write(iowrite,*)
  write(iowrite,'("  Calculating Bessel transforms may take a few seconds ")')
  write(iowrite,*)


  allocate(vlocft(0:nqmax))
  allocate(vkbprft(0:nqmax,0:mxdl,-1:1))
  allocate(cdcft(0:nqmax))
  allocate(cdvft(0:nqmax))
  allocate(basft(0:nqbas,norbas))

  call atom_kb_pot_four(npot, lo, irel, nicore, nr, r, drdi, zion,       &
         nql, nqnl, delql, nqbas, delqbas, vlocal, inorm, vkbproj,       &
         cdc, cdv, vlocft, vql0, vkbprft, cdcft, cdvft, basft,           &
         norbas, lo_b, rpsi_b,                                           &
         mxdl, mxdnr, nqmax, nqbas, norbas)

! pseudo file names

  if(icorr == 'pw' .or. icorr == 'ca' .or. icorr == 'PW' .or. icorr == 'CA') then
    scorr = '_LDA_'
  elseif(icorr == 'pb' .or. icorr == 'PB') then
    scorr = '_GGA_'
  else
    scorr = '_UNK_'
  endif

  if(nameat(1:1) == ' ') then
    filekb = nameat(2:2)//sfilekb
    fileupf = nameat(2:2)//scorr//sfileupf
  elseif(nameat(2:2) == ' ') then
    filekb = nameat(1:1)//sfilekb
    fileupf = nameat(1:1)//scorr//sfileupf
  else
    filekb = nameat//sfilekb
    fileupf = nameat//scorr//sfileupf
  endif

  call atom_kb_psd_out_four(iokb, filekb,                                &
         nameat, icorr, irel, nicore, irdate, irvers, irayps, psdtitle,  &
         nql, nqnl, delql, nqbas, delqbas, zion, vql0,                   &
         npot, lo, ev, inorm, vkbprft, vlocft, cdcft, cdvft,             &
         norbas, lo_b, basft,                                            &
         mxdl, nqmax, nqbas, norbas)

  call atom_kb_psd_out_upf(ioupf, fileupf,                               &
         nameat, icorr, irel, nicore, irdate, irvers, irayps, psdtitle,  &
         nr, r,  zion, vlocal, cdc, cdv,                                 &
         llocal, lmax_pot, vkbproj, inorm, lmax_pot, rpsi,               &
         mxdl, mxdnr)

  call atom_kb_plot(ioplot, fileplot, irel, nr, r, vlocal,               &
      rpsi(:,:,0), rpsi_b, nqbas, delqbas, ektot,                        &
      npot, lo, norbas, lo_b,                                            &
      inorm, vkbproj, nqnl, delql, vkbprft,                              &
      mxdnr, nqmax, mxdl, norbas)


  deallocate(r)
  deallocate(lo)
  deallocate(vionic)
  deallocate(cdc,cdv)

  deallocate(zo)
  allocate(rc)

  deallocate(ev)
  deallocate(rpsi)

  deallocate(drdi)
  deallocate(d2rodr)

  deallocate(vscreen)

  deallocate(ektot)

  deallocate(vlocal)

  deallocate(inorm)
  deallocate(vkbproj)

  deallocate(rpsi_b)
  deallocate(drpsidr_b)
  deallocate(veff_b)
  deallocate(lo_b))
  deallocate(ev_b)
  deallocate(nrc_b)

  deallocate(vlocft)
  deallocate(vkbprft)
  deallocate(cdcft)
  deallocate(cdvft)
  deallocate(basft))

  return

end subroutine atom_kb_sub

