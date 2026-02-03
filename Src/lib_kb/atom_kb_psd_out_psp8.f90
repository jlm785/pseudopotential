!>  writes the KB pseudopotential in the Hamman-abinit psp8 format
!>
!>  \author       J.L.Martins
!>  \version      6.1.0
!>  \date         29 January 2026.
!>  \copyright    GNU Public License v2

subroutine atom_kb_psd_out_psp8(iotape, fname,                           &
      nameat, icorr, irel, nicore, irvers, psdtitle,                     &
      nr, r,  zion, vlocal, cdc, cdv,                                    &
      llocal, lmax_pot, vkbproj, inorm, lmax_psi, rpsi,                  &
      mxdl, mxdnr)

! Adapted from atom_kb_psd_out_upf 29 January 2026. JLM
! Used linout.f90 from oncvpsp-3.3.1

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

  integer                 ::  isw                                        !  indicates if spin-orbit is included
  integer                 ::  ixc                                        !  indicates the type of exchange correlation
  integer                 ::  lloc                                       !  new coding of angular momentum of local pseudopotential

  character(len=10)       ::  ctemp
  integer                 ::  ierr

  real(REAL64)            ::  fchrg

  character(len=8)        ::  date
  character(len=10)       ::  time
  integer                 ::  idate

  integer                 ::  ifail

! allocatable local variables

  real(REAL64), allocatable        ::  rpsp8(:)                          !  radial grid points   rupf(i) = a*(i-1), i=1,...,nrpsp8
  real(REAL64), allocatable        ::  vlocpsp8(:)                       !  local pseudopotential (not r*vlocal) for psp8
  real(REAL64), allocatable        ::  vnlpsp8(:,:,:)                    !  kb-projectors (multiplied by r) for psp8
  real(REAL64), allocatable        ::  chi(:,:,:)                        !  wavefunctions (r(i),l,2j-2l) for psp8 (r*psi) (not used)
  real(REAL64), allocatable        ::  cdcpsp8(:)                        !  charge density of core.
  real(REAL64), allocatable        ::  cdvpsp8(:)                        !  4*pi*r**2 charge density of valence.

  integer, allocatable             ::  nproj(:)                          !  numbr of projectors for each angulat momentum (0 or 1)
  real(REAL64), allocatable        ::  r_core(:)                         !  core radius

  real(REAL64), allocatable        ::  dcdcdrpsp8(:,:)                   !  d^n cdcpsp8 / d rpsp8^n
  real(REAL64), allocatable        ::  err_max(:)                        !  estimate of maximum error of dcdcdrpsp8
  real(REAL64), allocatable        ::  cdcor2(:)                         !  core charge withiut the r**2 term

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  character(len=1), parameter  ::  IL(0:6) = (/'s','p','d','f','g','h','i'/)
  real(REAL64), parameter   ::  PI   = 3.14159265358979312_REAL64

! counters

  integer       ::   i, j, l


  if(iotape == 0) return

! abinit also suffers from the year2000 bug;-)

  call date_and_time(date, time )
  read(date(3:8),'(i6)') idate

! Hard coded values

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

  lso = .FALSE.
  isw = 0
  if(irel == 'rel') then
    lso = .TRUE.
    isw = 2
  endif

  ISW = 0     !    TEMPORARY BEFORE IMPLEMENTING SPIN ORBIT

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

! gets the configuration as psp8 needs

  allocate(nproj(0:lmax_pot))

  do l = 0,lmax_pot
    nproj(l) = 0
    if(inorm(l,-1) /= 0 .or. inorm(l, 0) /= 0 .or. inorm(l, 1) /= 0) nproj(l) = 1
  enddo

  if(llocal < 0) then
    lloc = lmax_pot+1
  else
    lloc = llocal
    nproj(lloc) = 0                        !  This has to be rethought for spin-orbit
  endif

! writes heading

  open (unit = iotape, file = trim(fname), form = 'formatted', status = 'unknown')

  write(iotape,'(1x,a2,4x,a10,4x,"rcore=  ",8f10.5)') nameat,            &
           irvers, (r_core(l),l=0,lmax_pot)
  write(iotape,'(2f12.4,4x,i6,"    zatom,zion,pspd")') iz*UM, zion, idate
  write(iotape,'(i6,i8,i4,3i6, "    pspcod,pspxc,lmax,lloc,mmax,r2w")')   &
        8, ixc, lmax_pot,lloc, nrpsp8, 0
  write(iotape,'(3f12.8,"    rchrg fchrg qchrg")') rpsp8(nrpsp8),fchrg, 0.0
  write(iotape,'(8i6,"    nproj")') (nproj(l),l=0,lmax_pot)
  write(iotape,'(i6,"              extension_switch")') isw

! write the VKB projectors and the local potential
! potential in hartree not rydberg

  do l = 0,max(lmax_pot,lloc)
    if(l == lloc) then
      write(iotape,'(i4)') l
      do i = 1,nrpsp8
        write(iotape,'(i6,1p,2d21.13)') i, rpsp8(i), vlocpsp8(i)/2
      end do
    else
      write(iotape,'(i4,23x,1p,d21.13)') l,(UM*inorm(l,0))/2
      do i = 1,nrpsp8
        write(iotape,'(i6,1p,2d21.13)') i,rpsp8(i),vnlpsp8(i,l,0)
      end do
    endif
  enddo

! calculates the derivatives of the partial core charge density

  allocate(dcdcdrpsp8(nrpsp8,0:4))
  allocate(err_max(0:4))
  allocate(cdcor2(1:nr))

  do i = 2,nr
    cdcor2(i) = cdc(i) / (r(i)*r(i))
  enddo
  cdcor2(1) = cdcor2(2)

  call grid_interp(6, 4, nr, r, cdcor2, nrpsp8, rpsp8, dcdcdrpsp8, err_max, ifail, nrpsp8)

!   do j = 0,4
!   do i = 1,nrpsp8
!     dcdcdrpsp8(i,j) = 4*PI*dcdcdrpsp8(i,j)
!   enddo
!   enddo

  DO I = 1,nrpsp8
    WRITE(87,'(2E13.5)') rpsp8(i), 4*PI*cdcpsp8(I)
  ENDDO

  DO I = 1,nrpsp8
    WRITE(88,'(2E13.5)') rpsp8(i), dcdcdrpsp8(i,0)
  ENDDO

  DO I = 1,nr
    WRITE(89,'(2E13.5)') r(i), cdcor2(i)
  ENDDO


  if(ifail /= 0) then

    write(6,*)
    write(6,*) '   WARNING '
    write(6,*)

    write(6,*) '   calculation of partial core charge derivatives failed'
    write(6,*) '   will set them to zero and hope for the best...'

    dcdcdrpsp8(:,:) = ZERO
    dcdcdrpsp8(1:nrpsp8,0) = 4*PI*cdcpsp8(1:nrpsp8)

  endif

! writes partial core charge

  do i = 1,nrpsp8
    write(iotape,'(i6,1p,6d21.13)') i, rpsp8(i), (dcdcdrpsp8(i,j),j=0,4)
    WRITE(90,*) i
  end do

  close(unit = iotape)

  deallocate(rpsp8)
  deallocate(nproj)

  deallocate(vlocpsp8)
  deallocate(vnlpsp8)
  deallocate(chi)
  deallocate(cdcpsp8)
  deallocate(cdvpsp8)

  deallocate(cdcor2)
  deallocate(dcdcdrpsp8)
  deallocate(err_max)

  return

end subroutine atom_kb_psd_out_psp8
