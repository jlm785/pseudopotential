!>  writes the KB pseudopotential in the Hamman-abinit psp8 format
!>
!>  \author       J.L.Martins+claude
!>  \version      6.1.0
!>  \date         29 January 2026, 25 October 2026.
!>  \copyright    GNU Public License v2

subroutine atom_kb_psd_out_psp8(iotape, fname,                           &
      nameat, icorr, irel, nicore, irvers, psdtitle,                     &
      nr, r,  zion, vlocal, cdc, cdv,                                    &
      llocal, lmax_pot, vkbproj, inorm, lmax_psi, rpsi,                  &
      mxdl, mxdnr)

! Adapted from atom_kb_psd_out_upf 29 January 2026. JLM
! spin-orbit, valence charge, normalized projectors. 25 October 2026. claude

! Used linout.f90 from oncvpsp-3.3.1. claude

! Format of the psp8 file as read by abinit (src/64_psp/m_psp8.F90),
! all reads are list directed:
!
!   title
!   zatom, zion, pspd
!   pspcod (=8), pspxc, lmax, lloc, mmax, r2well
!   rchrg, fchrg, qchrg
!   nproj(0:lmax)
!   extension_switch      0, 1 (+valence charge), 2 (spin-orbit), 3 (both)
!   nprojso(1:lmax)       only if extension_switch = 2 or 3
!   for l = 0,lmax:  if nproj(l) > 0:  l, ekb(1:nproj(l))
!                                      mmax lines: i, r, p_1(r),...
!                    else if l = lloc: l
!                                      mmax lines: i, r, vloc(r)
!   if lloc > lmax:  lloc
!                    mmax lines: i, r, vloc(r)
!   for l = 1,lmax with nprojso(l) > 0 (spin-orbit):
!                    l, ekbso(1:nprojso(l))
!                    mmax lines: i, r, pso_1(r),...
!   if fchrg > 0:    mmax lines: i, r, 4 pi rho_core and 4 derivatives
!   if extension_switch = 1 or 3:  mmax lines: i, r, 4 pi rho_valence
!
! The radial grid is linear starting at r = 0.  Energies in Hartree.
! The projectors p(r) = u(r) dV(r) are normalized to integral p**2 dr = 1,
! the non-local operator is sum_i |p_i> ekb_i <p_i| (+ L.S sum_i |pso_i> ekbso_i <pso_i|).
! lloc must not have projectors if lloc < 4.
!
! Spin-orbit: with p+, p- the projectors for j = l+1/2, l-1/2 (l > 0)
! and E+, E- their KB energies, V = V_SR + V_SO L.S with
!   V_SR:  p+ with (l+1)/(2l+1) E+,   p- with l/(2l+1) E-
!   V_SO:  p+ with 2/(2l+1) E+,       p- with -2/(2l+1) E-
! as in sr_so_r.f90 of ONCVPSP (D. R. Hamann).  All l may then have
! projectors, so a general local potential (lloc = 4) is used.

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! dimensions:

  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum components
  integer, intent(in)               ::  mxdnr                            !<  dimension of radial grid points.

! input:

  integer, intent(in)               ::  iotape                           !<  output tape number
  character(len=*), intent(in)      ::  fname                            !<  file name of psp8 potential

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbol of the element
  character(len=2), intent(in)      ::  icorr                            !<  correlation used in the calculation
  character(len=3), intent(in)      ::  irel                             !<  flag for relativistic (r) or spin-polarized (s) original calculation
  character(len=4), intent(in)      ::  nicore                           !<  flag for core correction
  character(len=10), intent(in)     ::  irvers                   !<  version of original calculation
  character(len=10), intent(in)     ::  psdtitle(20)                     !<  pseudopotential parameters

  integer, intent(in)               ::  nr                               !<  number of points in the radial grid
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points   r(i) = a*(exp(b*(i-1))-1), i=1,...,nr

  real(REAL64), intent(in)          ::  zion                             !<  ionic charge

  real(REAL64), intent(in)          ::  vlocal(mxdnr)                    !<  r*local pseudopotential
  real(REAL64), intent(in)          ::  cdc(mxdnr)                       !<  4*pi*r**2 charge density of valence.  kb_conv is not spin polarized...
  real(REAL64), intent(in)          ::  cdv(mxdnr)                       !<  4*pi*r**2 charge density of core

  integer, intent(in)               ::  llocal                           !<  angular momentum of local pseudopotential (<0 supremum)
  integer                           ::  lmax_pot                         !<  maximum angular momentum in potential
  real(REAL64), intent(in)          ::  vkbproj(mxdnr,0:mxdl,-1:1)       !<  kb-projector(r(i),l,2j-2l)
  integer, intent(in)               ::  inorm(0:mxdl,-1:1)               !<  sign of denominator of KB operator

  integer, intent(in)               ::  lmax_psi                         !<  maximum angular momentum in wave-functions
  real(REAL64), intent(in)          ::  rpsi(mxdnr,0:mxdl,-1:1)          !<  wavefunctions (r(i),l,2j-2l).

! local variables

  integer                 ::  iz
  integer                 ::  nrpsp8                                     !  number of points in the pseudo grid

  real(REAL64)            ::  delr

  logical                 ::  lso                                        !  indicates if spin-orbit is used (irel == 'rel')

  integer                 ::  isw                                        !  extension_switch, valence charge and spin-orbit
  integer                 ::  ixc                                        !  indicates the type of exchange correlation
  integer                 ::  lloc                                       !  new coding of angular momentum of local pseudopotential

  character(len=10)       ::  ctemp
  integer                 ::  ierr

  real(REAL64)            ::  fchrg

  character(len=8)        ::  date
  character(len=10)       ::  time
  integer                 ::  idate

  integer                 ::  ifail

  real(REAL64)            ::  pnorm                                      !  integral of projector squared
  real(REAL64)            ::  ekb                                        !  KB energy (Hartree) of normalized projector
  real(REAL64)            ::  fsr, fso                                   !  scalar and spin-orbit weights of j = l +/- 1/2

! allocatable local variables

  real(REAL64), allocatable        ::  rpsp8(:)                          !  radial grid points   rupf(i) = a*(i-1), i=1,...,nrpsp8
  real(REAL64), allocatable        ::  vlocpsp8(:)                       !  local pseudopotential (not r*vlocal) for psp8
  real(REAL64), allocatable        ::  vnlpsp8(:,:,:)                    !  kb-projectors (multiplied by r) for psp8
  real(REAL64), allocatable        ::  chi(:,:,:)                        !  wavefunctions (r(i),l,2j-2l) for psp8 (r*psi) (not used)
  real(REAL64), allocatable        ::  cdcpsp8(:)                        !  charge density of core.
  real(REAL64), allocatable        ::  cdvpsp8(:)                        !  4*pi*r**2 charge density of valence.

  integer, allocatable             ::  nproj(:)                          !  number of scalar projectors for each angular momentum
  integer, allocatable             ::  nprojso(:)                        !  number of spin-orbit projectors for each angular momentum
  real(REAL64), allocatable        ::  ekbsr(:,:)                        !  KB energies of scalar projectors (Hartree)
  real(REAL64), allocatable        ::  ekbso(:,:)                        !  KB energies of spin-orbit projectors (Hartree)
  real(REAL64), allocatable        ::  projsr(:,:,:)                     !  normalized scalar projectors
  real(REAL64), allocatable        ::  projso(:,:,:)                     !  normalized spin-orbit projectors
  real(REAL64), allocatable        ::  r_core(:)                         !  core radius

  real(REAL64), allocatable        ::  dcdcdrpsp8(:,:)                   !  d^n cdcpsp8 / d rpsp8^n
  real(REAL64), allocatable        ::  err_max(:)                        !  estimate of maximum error of dcdcdrpsp8
  real(REAL64), allocatable        ::  cdcor2(:)                         !  core charge withiut the r**2 term

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  character(len=1), parameter  ::  IL(0:6) = (/'s','p','d','f','g','h','i'/)
  real(REAL64), parameter   ::  PI   = 3.14159265358979312_REAL64
  real(REAL64), parameter   ::  TINY = 1.0E-50_REAL64                    !  smaller values are written as zero (3 digit exponents lose the D)

! counters

  integer       ::   i, j, k, l, n


  if(iotape == 0) return

! abinit also suffers from the year2000 bug;-)

  call date_and_time(date, time )
  read(date(3:8),'(i6)') idate

! Hard coded values (nrpsp8 must be odd for Simpson integration)

  nrpsp8 = 2001
  delr = 0.005_REAL64

  fchrg = ZERO
  if(nicore == 'nc  ') then
    fchrg = ZERO
  elseif(nicore == 'pcec') then
    fchrg = UM
  else

    write(6,*) '  Unrecognized partial core correction, not writing abinit pseudo'

    return

  endif

  if(icorr == 'pb') then
    ixc = 11
  elseif (icorr == 'ca') then
    ixc = 2
  elseif (icorr == 'pw') then
    ixc = 7
  else

    write(6,*) '  Unrecognized correlation, not writing abinit pseudo'

    return

  endif

! the valence charge is always written

  lso = .FALSE.
  isw = 1
  if(irel == 'rel') then
    lso = .TRUE.
    isw = 3
  endif

  call atom_p_tbl_charge(nameat,iz)

! reinvents the mesh as linear

  allocate(rpsp8(nrpsp8))

  do j = 1,nrpsp8
    rpsp8(j) = delr*(j-1)
  enddo

! interpolates on the psp8 mesh the local potential and non-local pseudopotentials,
! the wavefunctions and charge densities (core and valence)

  allocate(vlocpsp8(nrpsp8))
  allocate(vnlpsp8(nrpsp8,0:mxdl,-1:1))
  allocate(chi(nrpsp8,0:mxdl,-1:1))
  allocate(cdcpsp8(nrpsp8))
  allocate(cdvpsp8(nrpsp8))

  allocate(r_core(0:mxdl))

! recover information from psdtitle

  do l = 0,lmax_pot
    r_core(l) = ZERO
    ctemp = psdtitle(l*2+2)
    if(ctemp(1:1) == ' ' .or. ctemp(1:1) == 'r') then
      read(ctemp(6:10),*,iostat=ierr) r_core(l)
    else
      read(ctemp(7:10),*,iostat=ierr) r_core(l)
    endif
  enddo

  call atom_kb_psd_out_interp(nr, r, nrpsp8, rpsp8, lso,                 &
      vlocal, lmax_pot, vkbproj, lmax_psi, rpsi, cdv, cdc,               &
      vlocpsp8, vnlpsp8, chi, cdvpsp8, cdcpsp8,                          &
      mxdl, mxdnr)

! local potential

  if(lso) then
    lloc = max(4,lmax_pot+1)
  elseif(llocal < 0) then
    lloc = lmax_pot+1
  else
    lloc = llocal
  endif

! normalized projectors and their KB energies in Hartree

  allocate(nproj(0:lmax_pot), nprojso(0:lmax_pot))
  allocate(ekbsr(2,0:lmax_pot), ekbso(2,0:lmax_pot))
  allocate(projsr(nrpsp8,2,0:lmax_pot), projso(nrpsp8,2,0:lmax_pot))

  nproj(:) = 0
  nprojso(:) = 0
  ekbsr(:,:) = ZERO
  ekbso(:,:) = ZERO
  projsr(:,:,:) = ZERO
  projso(:,:,:) = ZERO

  do l = 0,lmax_pot
    do k = 1,2

!     k = 1:  j = l+1/2 (or no spin-orbit),  k = 2:  j = l-1/2

      if(lso) then
        j = 3 - 2*k
        if(l == 0 .and. j == -1) cycle
      else
        if(k == 2 .or. l == lloc) cycle
        j = 0
      endif

      if(inorm(l,j) == 0) cycle

!     Simpson integration on the uniform grid

      pnorm = vnlpsp8(1,l,j)**2 + vnlpsp8(nrpsp8,l,j)**2
      do i = 2,nrpsp8-1,2
        pnorm = pnorm + 4*vnlpsp8(i,l,j)**2
      enddo
      do i = 3,nrpsp8-2,2
        pnorm = pnorm + 2*vnlpsp8(i,l,j)**2
      enddo
      pnorm = pnorm*delr/3

!     V_KB = inorm |vkbproj><vkbproj| in Rydberg

      ekb = (inorm(l,j)*pnorm)/2

      if(.not. lso .or. l == 0) then
        fsr = UM
        fso = ZERO
      elseif(j == 1) then
        fsr = (l+1)/(2*l+UM)
        fso = 2/(2*l+UM)
      else
        fsr = l/(2*l+UM)
        fso = -2/(2*l+UM)
      endif

      nproj(l) = nproj(l) + 1
      n = nproj(l)
      ekbsr(n,l) = fsr*ekb
      do i = 1,nrpsp8
        projsr(i,n,l) = vnlpsp8(i,l,j) / sqrt(pnorm)
      enddo

      if(fso /= ZERO) then
        nprojso(l) = nprojso(l) + 1
        n = nprojso(l)
        ekbso(n,l) = fso*ekb
        do i = 1,nrpsp8
          projso(i,n,l) = vnlpsp8(i,l,j) / sqrt(pnorm)
        enddo
      endif

    enddo
  enddo

  where(abs(projsr) < TINY) projsr = ZERO
  where(abs(projso) < TINY) projso = ZERO

! writes heading

  open (unit = iotape, file = trim(fname), form = 'formatted', status = 'unknown')

  write(iotape,'(1x,a2,4x,a10,4x,"rcore=  ",8f10.5)') nameat,            &
           irvers, (r_core(l),l=0,lmax_pot)
  write(iotape,'(2f12.4,4x,i6,"    zatom,zion,pspd")') iz*UM, zion, idate
  write(iotape,'(i6,i8,i4,3i6, "    pspcod,pspxc,lmax,lloc,mmax,r2w")')   &
        8, ixc, lmax_pot,lloc, nrpsp8, 0
  write(iotape,'(3f12.8,"    rchrg fchrg qchrg")') rpsp8(nrpsp8),fchrg, 0.0
  write(iotape,'(8i6)', advance='no') (nproj(l),l=0,lmax_pot)
  write(iotape,'("    nproj")')
  write(iotape,'(i6,"              extension_switch")') isw
  if(lso) then
    write(iotape,'(8i6)', advance='no') (nprojso(l),l=1,lmax_pot)
    write(iotape,'("    nprojso")')
  endif

! write the scalar projectors and the local potential
! potential in hartree not rydberg

  do l = 0,lmax_pot
    if(nproj(l) > 0) then
      write(iotape,'(i4,23x,1p,2d21.13)') l, (ekbsr(k,l),k=1,nproj(l))
      do i = 1,nrpsp8
        write(iotape,'(i6,1p,3d21.13)') i, rpsp8(i), (projsr(i,k,l),k=1,nproj(l))
      end do
    elseif(l == lloc) then
      write(iotape,'(i4)') l
      do i = 1,nrpsp8
        write(iotape,'(i6,1p,2d21.13)') i, rpsp8(i), vlocpsp8(i)/2
      end do
    endif
  enddo

  if(lloc > lmax_pot) then
    write(iotape,'(i4)') lloc
    do i = 1,nrpsp8
      write(iotape,'(i6,1p,2d21.13)') i, rpsp8(i), vlocpsp8(i)/2
    end do
  endif

! write the spin-orbit projectors

  if(lso) then
    do l = 1,lmax_pot
      if(nprojso(l) > 0) then
        write(iotape,'(i4,23x,1p,2d21.13)') l, (ekbso(k,l),k=1,nprojso(l))
        do i = 1,nrpsp8
          write(iotape,'(i6,1p,3d21.13)') i, rpsp8(i), (projso(i,k,l),k=1,nprojso(l))
        end do
      endif
    enddo
  endif

! partial core charge (4 pi rho_core) and its derivatives

  if(fchrg > ZERO) then

    allocate(dcdcdrpsp8(nrpsp8,0:4))
    allocate(err_max(0:4))
    allocate(cdcor2(1:nr))

    do i = 2,nr
      cdcor2(i) = cdc(i) / (r(i)*r(i))
    enddo
    cdcor2(1) = cdcor2(2)

    call grid_interp(6, 4, nr, r, cdcor2, nrpsp8, rpsp8, dcdcdrpsp8, err_max, ifail, nrpsp8)

    if(ifail /= 0) then

      write(6,*)
      write(6,*) '   WARNING '
      write(6,*)

      write(6,*) '   calculation of partial core charge derivatives failed'
      write(6,*) '   will set them to zero and hope for the best...'

      dcdcdrpsp8(:,:) = ZERO
      do i = 2,nrpsp8
        dcdcdrpsp8(i,0) = cdcpsp8(i) / (rpsp8(i)*rpsp8(i))
      enddo
      dcdcdrpsp8(1,0) = dcdcdrpsp8(2,0)

    endif

    where(abs(dcdcdrpsp8) < TINY) dcdcdrpsp8 = ZERO

    do i = 1,nrpsp8
      write(iotape,'(i6,1p,6d21.13)') i, rpsp8(i), (dcdcdrpsp8(i,j),j=0,4)
    end do

    deallocate(cdcor2)
    deallocate(dcdcdrpsp8)
    deallocate(err_max)

  endif

! valence charge (4 pi rho_valence)

  do i = 2,nrpsp8
    cdvpsp8(i) = cdvpsp8(i) / (rpsp8(i)*rpsp8(i))
  enddo
  cdvpsp8(1) = cdvpsp8(2)
  where(abs(cdvpsp8) < TINY) cdvpsp8 = ZERO

  do i = 1,nrpsp8
    write(iotape,'(i6,1p,2d21.13)') i, rpsp8(i), cdvpsp8(i)
  end do

  close(unit = iotape)

  deallocate(rpsp8)
  deallocate(nproj, nprojso)
  deallocate(ekbsr, ekbso)
  deallocate(projsr, projso)
  deallocate(r_core)

  deallocate(vlocpsp8)
  deallocate(vnlpsp8)
  deallocate(chi)
  deallocate(cdcpsp8)
  deallocate(cdvpsp8)

  return

end subroutine atom_kb_psd_out_psp8
