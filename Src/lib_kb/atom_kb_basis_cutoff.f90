!>  Chooses the cut-off for the atomic basis functions
!>  based on the SIESTA radii (energy shift ESHIFT) and
!>  the rule in atom_kb_basis_rule
!>
!>  \author       J.L.Martins+claude
!>  \version      6.1.0
!>  \date         30 August 2021, 3 November 2026.
!>  \copyright    GNU Public License v2

subroutine atom_kb_basis_cutoff(n_bsets, tbasis, lmax_pot, nameat,       &
         lmax_bas, n_bas, r_bas, nz_bas, r_siesta, r_99,                 &
         iowrite, ioreal, filereal,                                      &
         mxdl, mxdset, mxdnr)

! nrmax -> mxdnr.  17 September 2021. JLM
! jhard. 21 december 2021. JLM
! psdtitle, ifcore. 19 May 2022. JLM
! fileparsec -> filereal. 26 January 2026. JLM
! highest available l potential instead of l = 0. 24 October 2026. JLM+claude
! second zeta is the relaxed energy derivative (r_bas = 0), up to 4 zeta. 24 October 2026. JLM+claude
! r_siesta for all bound l, eshift parameter, basis from atom_kb_basis_rule. 3 November 2026. claude+JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdl                             !<  array dimension for basis angular momentum
  integer, intent(in)               ::  mxdnr                            !<  dimension of radial grid points.
  integer, intent(in)               ::  mxdset                           !<  dimension for number of atomic basis sets

  integer, intent(in)               ::  iowrite                          !<  default output tape

  integer, intent(in)               ::  ioreal                           !<  default tape for pseudopotential in real space (parse)c format
  character(len=*), intent(in)      ::  filereal                         !<  name of default tape for reading pseudopotential in real space (parsec) format

  integer, intent(in)               ::  n_bsets                          !<  number of atomic basis sets
  character(len=3), intent(in)      ::  tbasis(mxdset)                   !<  type of basis

  integer,intent(in)                ::  lmax_pot                         !<  maximum angular momentum for potential

! output

  character(len=2),intent(out)      ::  nameat                           !<  chemical symbol of the element

  integer,intent(out)               ::  lmax_bas(mxdset)                 !<  maximum angular momentum in basis
  integer,intent(out)               ::  n_bas(0:mxdl,mxdset)             !<  basis functions for angular momentum l
  real(REAL64), intent(out)         ::  r_bas(4,0:mxdl,mxdset)           !<  cutoff for the basis (up to quadruple zeta), 0: energy derivative of previous
  integer, intent(out)              ::  nz_bas(4,0:mxdl,mxdset)          !<  number of non-trivial zeroes in basis function
  real(REAL64), intent(out)         ::  r_siesta(0:mxdl)                 !<  cutoff using the SIESTA recipe
  real(REAL64), intent(out)         ::  r_99(0:mxdl)                     !<  radius with 99% of charge

! variables from the output file of psd_gen

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
  real(REAL64), allocatable         ::  tauc_total(:)                    !  core kinetic energy density

! variables for KB-pseudopotential

  real(REAL64), allocatable         ::  drdi(:)                          !  d r(i) / d i
  real(REAL64), allocatable         ::  d2rodr(:)                        !  (d^2 r / d i^2) / (d r / d i)

  integer                           ::  ifcore                           !  derived from nicore
  real(REAL64)                      ::  totvel                           !  total number of valence electrons
  real(REAL64), allocatable         ::  vscreen(:)                       !  screening potential

  real(REAL64), allocatable         ::  veff(:)                          !  effective potential

! variables for wavefunctions

  real(REAL64), allocatable         ::  ev(:,:)                          !  eigenvalues (l,2j-2l).
  real(REAL64), allocatable         ::  rpsi(:,:,:)                      !  wavefunctions (r(i),l,2j-2l).

  real(REAL64), allocatable         ::  zo(:)                            !  orbital occupation (negative if abesent)
  real(REAL64), allocatable         ::  rc(:)                            !  core radius

  real(REAL64), allocatable         ::  rpsib(:)                         !  basis functions
  real(REAL64), allocatable         ::  drpsibdr(:)                      !  d rpsib / d r


! other local variables

  real(REAL64)                      ::  rtry

  real(REAL64)                      ::  revi                             !  wave-function calculated up to revi
  integer                           ::  iflag                            !  iflag = 0: success; iflag = 1: failed to converge, iflag > 3: major error
  integer                           ::  nrevi                            !  wave-functions was calculated up to r(nrevi)

  integer                           ::  l_occ                            !  highest l of a bound occupied orbital

  real(REAL64)                      ::  factor
  integer                           ::  ll
  character(len=3)                  ::  ch_temp, ch_temp2

! constants

  real(REAL64), parameter           ::  ZERO = 0.0_REAL64
  real(REAL64), parameter           ::  SMALL = 1.0E-6_REAL64

! Hard coded constants

  real(REAL64), parameter           ::  ESHIFT = 0.02_REAL64             !  energy shift (Ry) of the SIESTA recipe

! counters

  integer                           ::  i, j, l


  INTEGER   ::  LMAX

! paranoid check

  if(lmax_pot > mxdl) then
    write(6,*)
    write(6,*) '  Stopped in atom_kb_basis_cutoff:'
    write(6,*) '  lmax_pot = ',lmax_pot,' mxdl = ',mxdl

    STOP

  endif

! check consistency of tbasis

  if(n_bsets > 1) then
    do i = 1,n_bsets-1
    do j = 2,n_bsets
      ch_temp = tbasis(i)
      ch_temp2 = tbasis(j)
      call chrcap(ch_temp,3)
      call chrcap(ch_temp2,3)
      if(ch_temp == ch_temp2) then
        write(6,*) '   STOPPED in atom_kb_basis_cutoff'
        write(6,*) '   Duplicate kind of atomic basis set', i, j
        write(6,*) tbasis(i), tbasis(j)

        stop

      endif
    enddo
    enddo
  endif

! check known type of basis

  do i = 1,n_bsets
    ch_temp = tbasis(i)
    call chrcap(ch_temp,3)
    if(ch_temp(1:2) /= 'SZ' .and. ch_temp(1:2) /= 'DZ') then
        write(6,*) '   STOPPED in atom_kb_basis_cutoff'
        write(6,*) '   Unrecognized type of atomic basis set', ch_temp

        stop

    endif
  enddo

! initialize

  do l = 0,mxdl
    r_siesta(l) = ZERO
    r_99(l) = ZERO
    do i = 1,mxdset
      n_bas(l,i) = 0
      do j = 1,4
        r_bas(j,l,i) = ZERO
        nz_bas(j,l,i) = 0
      enddo
    enddo
  enddo

! allocations

  allocate(r(mxdnr))
  allocate(lo(mxdl+1,-1:1))
  allocate(vionic(mxdnr,0:mxdl,-1:1))
  allocate(cdc(mxdnr),cdv(mxdnr))
  allocate(tauc_total(mxdnr))

  allocate(ev(0:mxdl,-1:1))
  allocate(rpsi(mxdnr,0:mxdl,-1:1))

  allocate(zo(0:mxdl))
  allocate(rc(0:mxdl))

! initialize to zero

  r = ZERO
  lo = 0
  vionic = ZERO
  cdc = ZERO
  cdv = ZERO

! reads the pseudopotential data file

  call atom_kb_psd_in_real(ioreal, filereal,                             &
         nameat, icorr, irel, nicore, irdate, irvers, irayps, psdtitle,  &
         npot, nr, a, b, r, zion, lo, vionic, cdc, cdv, tauc_total,      &
         zo, rc, rpsi(:,:,0),                                            &
         lmax, mxdnr, mxdl)

  allocate(drdi(mxdnr))
  allocate(d2rodr(mxdnr))

! paranoid check and auxiliary arrays

  do j = 1,nr
    rtry = a*(exp(b*(j-1))-1)

    if(abs(rtry - r(j)) > SMALL) then
      write(iowrite,'("   INCOMPATIBLE  a, b r(j) ")')

      STOP

    endif

    drdi(j) = (r(j)+a)*b
    d2rodr(j) = b
  enddo

  allocate(vscreen(mxdnr))

  ifcore = 0
  if(nicore == 'fcec'.or.nicore == 'pcec') ifcore = 1
  if(nicore == 'fche'.or.nicore == 'pche') ifcore = 2

  call atom_kb_screen(nr, r, drdi, cdv, cdc, icorr, ifcore, totvel, vscreen,   &
         iowrite, mxdnr)

  call atom_kb_wvfct(npot, lo, irel, nr, r, drdi, d2rodr,                &
         vionic, vscreen, ev, rpsi,                                      &
         iowrite, mxdl, mxdnr)

  write(iowrite,'(/,"   Orbital energies",//,5x,"n",5x,"l",13x,"e",/)')
  do i=1,npot(0)
    write(iowrite,'(i6,i6,4x,f14.6)') lo(i,0)+1,lo(i,0),ev(lo(i,0),0)
  enddo

! highest angular momentum of a bound occupied orbital

  allocate(veff(mxdnr))
  allocate(rpsib(mxdnr))
  allocate(drpsibdr(mxdnr))

  l_occ = 0
  do l = lmax_pot,0,-1
    if( ev(l,0) < -SMALL .and. zo(l) > SMALL) then
      l_occ = l

      exit

    endif
  enddo

! r_99 and r_siesta for all l of the potential (r_siesta only for bound orbitals)

  do l = 0,lmax_pot

    factor = ZERO
    ll = 4
    do j = 2,nr
      factor = factor + ll*rpsi(j,l,0)*rpsi(j,l,0)*drdi(j) / 3
       ll = 6 - ll
       r_99(l) = r(j)

       if(factor > 0.99) exit

    enddo

    if(ev(l,0) >= -SMALL) cycle

!   uses the highest available l in the case potential is not available

    ll = -1
    do i = 1,npot(0)
      if(lo(i,0) == l) ll = l
    enddo
    if(ll < 0) then
      do i = 1,npot(0)
        if(lo(i,0) > ll) ll = lo(i,0)
      enddo
    endif

    do j = 2,nr
      veff(j) = (vionic(j,ll,0)  + (l*(l+1)) /r(j)) / r(j)  + vscreen(j)
    enddo
    veff(1) = veff(2)

!   siesta recipe: last node of the solution at energy ev + ESHIFT

    revi = 1.5*r_99(l)
    call atom_atm_difnrl_one(nr, r, drdi, d2rodr, veff, rpsib, drpsibdr,   &
        l, ev(l,0)+ESHIFT, revi, nrevi, iflag,                             &
        iowrite, mxdnr)

    do j = 2,nrevi-1
      if(rpsib(j)*rpsib(j-1) < ZERO) then
        r_siesta(l) = (rpsib(j)*r(j-1)-rpsib(j-1)*r(j)) / (rpsib(j)-rpsib(j-1))
      endif
    enddo

  enddo

! angular momenta and cutoff radii of the atomic basis sets

  call atom_kb_basis_rule(n_bsets, tbasis, lmax_pot, l_occ,              &
      r_siesta, lmax_bas, n_bas, r_bas, nz_bas,                          &
      iowrite, mxdl, mxdset)

  deallocate(r)
  deallocate(lo)
  deallocate(vionic)
  deallocate(cdc,cdv)

  deallocate(zo)
  deallocate(rc)

  deallocate(ev)
  deallocate(rpsi)

  deallocate(vscreen)

  deallocate(veff)
  deallocate(rpsib)
  deallocate(drpsibdr)

  deallocate(drdi)
  deallocate(d2rodr)

  return

end subroutine atom_kb_basis_cutoff

