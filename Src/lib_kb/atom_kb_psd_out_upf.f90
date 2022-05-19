!>  writes the KB pseudopotential in the Quantum Espresso UPF format
!>
!>  \author       J.L.Martins
!>  \version      6.0.8
!>  \date         October 2018. 25 January 2022. 19 May 2022.
!>  \copyright    GNU Public License v2

subroutine atom_kb_psd_out_upf(iotape, fname,                            &
      nameat, icorr, irel, nicore, irdate, irvers, irayps, psdtitle,     &
      nr, r,  zion, vlocal, cdc, cdv,                                    &
      llocal, lmax_pot, vkbproj, inorm, lmax_psi, rpsi,                  &
      mxdl, mxdnr)

! Adapted by JLMartins October 2018
! Converted to new style. 25 January 2022. JLM
! psdtitle, bug/feature reported by Raymond Atta-Fynn, 19 May 2022. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! dimensions:

  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum components
  integer, intent(in)               ::  mxdnr                            !<  dimension of radial grid points.

! input:

  integer, intent(in)               ::  iotape                           !<  output tape number
  character(len=*), intent(in)      ::  fname                            !<  file name of UPF potential

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbol of the element
  character(len=2), intent(in)      ::  icorr                            !<  correlation used in the calculation
  character(len=3), intent(in)      ::  irel                             !<  flag for relativistic (r) or spin-polarized (s) original calculation
  character(len=4), intent(in)      ::  nicore                           !<  flag for core correction
  character(len=10), intent(in)     ::  irdate, irvers                   !<  date and version of original calculation
  character(len=10), intent(in)     ::  irayps(4)                        !<  type of pseudopotential
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
  integer                 ::  nrupf
  real(REAL64)            ::  xmin
  real(REAL64)            ::  rmax, a, b

  character(len=2)        ::  nl
  real(REAL64)            ::  occup, rc

  integer                 ::  nproj

  character(len=10)       ::  ctemp

! allocatable local variables

  real(REAL64), allocatable        ::  rupf(:), rabupf(:)
  real(REAL64), allocatable        ::  vlocupf(:)
  real(REAL64), allocatable        ::  vnlupf(:,:)
  real(REAL64), allocatable        ::  dij(:,:)
  real(REAL64), allocatable        ::  chi(:,:)
  real(REAL64), allocatable        ::  cdcupf(:)
  real(REAL64), allocatable        ::  cdvupf(:)

! spline interpolation

  real(REAL64), allocatable        ::  yp(:),ypp(:),wlu(:,:)
  real(REAL64), allocatable        ::  ypi(:),yppi(:)
  integer                          ::  ierr

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter    ::  PI4 = 16*atan(1.0_REAL64)

! counters

  integer       ::   i, j, l, np


  if(iotape == 0) return

! reinvents the mesh as purely logarithmic

  call atom_p_tbl_charge(nameat,iz)

  xmin = 4*UM
  a = exp(-xmin)/(iz*UM)
!  a = exp(-8*UM)/(iz*UM)
  b = UM/(40.0*UM)
!  b = UM/(80.0*UM)
  rmax = 61*UM
!  rmax = 120*UM
  nrupf = 1 + nint(log(rmax/a + UM)/b)

  allocate(rupf(nrupf),rabupf(nrupf))

  do j = 1,nrupf
    rupf(j) = a*(exp(b*(j-1)))
    rabupf(j) = rupf(j)*b
  enddo

! interpolates on the upf mesh.

  allocate(vlocupf(nrupf))
  allocate(vnlupf(nrupf,0:mxdl))
  allocate(dij(mxdl+1,mxdl+1))
  allocate(chi(nrupf,0:mxdl))
  allocate(cdcupf(nrupf))
  allocate(cdvupf(nrupf))

  allocate(yp(nr),ypp(nr),wlu(nr,3))
  allocate(ypi(nrupf),yppi(nrupf))

! interpolates local potential

  call splift (r, vlocal, yp, ypp, nr, wlu, ierr, 0, ZERO,ZERO,ZERO,ZERO)

  if(ierr /= 1) then
    write(6,*) '   error in splift'
    stop
  endif

  call splint (r, vlocal, ypp, nr, rupf, vlocupf, ypi, yppi, nrupf, ierr)

  if(ierr /= 1) then
    write(6,*) '   error in splint'
    stop
  endif

  do i = 1,nrupf
    vlocupf(i) = vlocupf(i)/rupf(i)
  enddo

! interpolates KB projectors

  do l = 0,lmax_pot

    call splift (r, vkbproj(:,l,0), yp, ypp, nr, wlu, ierr, 1, ZERO,ZERO,ZERO,ZERO)

    if(ierr /= 1) then
      write(6,*) '   error in splift'
      stop
    endif

    call splint (r, vkbproj(:,l,0), ypp, nr, rupf, vnlupf(1,l), ypi, yppi, nrupf, ierr)

    if(ierr /= 1) then
      write(6,*) '   error in splint'
      stop
    endif

    do i = 1,nrupf
      vnlupf(i,l) = vnlupf(i,l)*rupf(i)
    enddo

  enddo

  if(llocal < 0) then
    nproj = lmax_pot+1
  else
    nproj = lmax_pot
  endif

  np = 0
  do l = 0,lmax_pot
    if(inorm(l,0) /= 0) then
      np = np + 1
    endif
  enddo

  if(np /= nproj) then
    write(6,*) '  stopped in kb_psd_out_upf'
    write(6,*) '  inconsistent number of projectors',np,nproj

    STOP

  endif

  do i = 1,nproj
    do j = 1,nproj
      dij(i,j) = ZERO
    enddo
  enddo
  np = 0
  do l = 0,lmax_pot
    if(inorm(l,0) /= 0) then
      np = np + 1
      dij(np,np) = UM * inorm(l,0)
    endif
  enddo

! interpolates wave-functions

  do l = 0,lmax_psi

    call splift (r, rpsi(:,l,0), yp, ypp, nr, wlu, ierr, 1, ZERO,ZERO,ZERO,ZERO)

    if(ierr /= 1) then
      write(6,*) '   error in splift'
      stop
    endif

    call splint (r, rpsi(:,l,0), ypp, nr, rupf, chi(:,l), ypi, yppi, nrupf, ierr)

    if(ierr /= 1) then
      write(6,*) '   error in splint'
      stop
    endif

  enddo

! interpolates charge densities

  call splift (r, cdv, yp, ypp, nr, wlu, ierr, 1, ZERO,ZERO,ZERO,ZERO)

  if(ierr /= 1) then
    write(6,*) '   error in splift'
    stop
  endif

  call splint (r, cdv, ypp, nr, rupf, cdvupf, ypi, yppi, nrupf, ierr)

  if(ierr /= 1) then
    write(6,*) '   error in splint'
    stop
  endif


  call splift (r, cdc, yp, ypp, nr, wlu, ierr, 1, ZERO,ZERO,ZERO,ZERO)

  if(ierr /= 1) then
    write(6,*) '   error in splift'
    stop
  endif

  call splint (r, cdc, ypp, nr, rupf, cdcupf, ypi, yppi, nrupf, ierr)

  if(ierr /= 1) then
    write(6,*) '   error in splint'
    stop
  endif

  do i = 1,nrupf
    cdcupf(i) = cdcupf(i) / (PI4*rupf(i)*rupf(i))
  enddo

! Writes the file

!  filename = adjustl(trim(nameat))//"_LDA_ncpp.UPF"

  open (unit = iotape, file = trim(fname), form = 'formatted', status = 'unknown')

  write(iotape,'(a21)') "<UPF version='2.0.1'>"

  write(iotape,'(a9)') "<PP_INFO>"

    write(iotape,'("    Generated with atom program version:  ",4x,a10)') irvers
    write(iotape,'("    Author: Froyen,Troullier,Martins et al.")')
    write(iotape,'("    Generation date:  ",a10)')  irdate
    write(iotape,'("    Pseudopotential type:   ",4a10)') (irayps(j),j=1,4)
    write(iotape,'("    Element: ",a2)') nameat
    write(iotape,'("    Functional: ",a2)') icorr
    write(iotape,'("    Suggested minimum cutoff for wavefunctions:",    &
         &    f10.4, " Ry")')  ZERO
    write(iotape,'("    Suggested minimum cutoff for charge density:",   &
         &    f10.4, " Ry")')  ZERO

    if(irel == 'nrl') then
      write(iotape,'("    The Pseudo was generated with a Non-",         &
         &    "Relativistic Calculation")')
    elseif(irel == 'isp') then
      write(iotape,'("    The Pseudo was generated with a spin-",        &
         &    "polarized Calculation")')
    else
      write(iotape,'("    The Pseudo was generated with a ",            &
         &    "Relativistic Calculation")')
    endif

!   the use of psdtitle is an hack.  Should be done right...

    if(llocal >= 0) then
      ctemp = psdtitle(2*llocal+2)
      rc = ZERO
      read(ctemp(6:10),*,iostat=ierr) rc
      write(iotape,'("    L of local component and cutoff radius",      &
           5x,i3,f11.5)') llocal, rc
    else
       write(iotape,'("    Local component is maximum of all",          &
            " potentials (smoothed)")')
    endif

    write(iotape,'("    Valence configuration:")')
    write(iotape,'("    nl",9x,"occ",12x,"Rcut")')
    do l = 0,lmax_pot
      ctemp = psdtitle(l*2+1)
      nl = '  '
      read(ctemp(1:2),'(a2)',iostat=ierr) nl
      occup = ZERO
      read(ctemp(4:9),*,iostat=ierr) occup
      ctemp = psdtitle(l*2+2)
      rc = ZERO
      read(ctemp(6:10),*,iostat=ierr) rc
      if(nl /= '  ') then
        write(iotape,'(4x,a2,5x,f10.4,5x,f10.4)') nl,occup,rc
      endif
    enddo
    write(iotape,'("    Generation configuration: not available")')

  write(iotape,'(a10)') "</PP_INFO>"


  write(iotape,'(a4)') "<!--"
  write(iotape,'(a29)') "END OF HUMAN READABLE SECTION"
  write(iotape,'(a4)') "-->"


  write(iotape,'(a10)') "<PP_HEADER"

    write(iotape,'("  generated=''Generated by the atomic Kleinman",    &
       "-Bylander code version ",a10,"''")') irvers
    write(iotape,'("  author=''Froyen,Troullier,Martins et al.''")')
    write(iotape,'("  date=''",a10,"''")')  irdate
    write(iotape,'("  comment=''using alpha code, beware''")')
    write(iotape,'("  element=''",a2,"''")') nameat
    write(iotape,'("  pseudo_type=''NC''")')
    if(irel == 'rel') then
      write(iotape,'("  relativistic=''full''")')
    else
      write(iotape,'("  relativistic=''nonrelativistic''")')
    endif
    write(iotape,'("  is_ultrasoft=''F''")')
    write(iotape,'("  is_paw=''F''")')
    write(iotape,'("  is_coulomb=''F''")')
    if(irel == 'rel') then
      write(iotape,'("  has_so=''T''")')
    else
      write(iotape,'("  has_so=''F''")')
    endif
    write(iotape,'("  has_wfc=''F''")')
    write(iotape,'("  has_gipaw=''F''")')
    write(iotape,'("  paw_as_gipaw=''F''")')
    if(nicore == 'nc  ') then
      write(iotape,'("  core_correction=''F''")')
    else
      write(iotape,'("  core_correction=''T''")')
    endif
    if(icorr == 'pb') then
      write(iotape,'("  functional='' SLA  PW   PBX PBC''")')
    elseif (icorr == 'ca') then
      write(iotape,'("  functional='' SLA  PZ   NOGX NOGC''")')
    elseif (icorr == 'pw') then
      write(iotape,'("  functional='' SLA  PW   NOGX NOGC''")')
    else
      write(6,*) '  Unrecognized correlation'
      stop
    endif
    write(iotape,'("  z_valence=''",e22.15,"''")') zion
    write(iotape,'("  total_psenergy=''",e22.15,"''")') ZERO
    write(iotape,'("  wfc_cutoff=''",e22.15,"''")') ZERO
    write(iotape,'("  rho_cutoff=''",e22.15,"''")') ZERO
    write(iotape,'("  l_max=''",i2,"''")') lmax_pot
    write(iotape,'("  l_max_rho=''0''")')
    write(iotape,'("  l_local=''",i2,"''")') llocal
    write(iotape,'("  mesh_size=''",i6,"''")') nrupf
    write(iotape,'("  number_of_wfc=''",i2,"''")') lmax_pot+1
    write(iotape,'("  number_of_proj=''",i2,"''")') nproj
  write(iotape,'(a2)') "/>"


  write(iotape,'(a8)') "<PP_MESH"
    write(iotape,'("  dx=''",e22.15,"''")') b
    write(iotape,'("  mesh=''",i6,"''")') nrupf
    write(iotape,'("  xmin=''",e22.15,"''")') xmin
    write(iotape,'("  rmax=''",e22.15,"''")') rupf(nrupf)
    write(iotape,'("  zmesh=''",e22.15,"''")') iz*UM
    write(iotape,'(a1)') ">"
    write(iotape,'(a6)') "<PP_R>"
      write(iotape,'(4(1x,e24.15,1x))') rupf
    write(iotape,'(a7)') "</PP_R>"
    write(iotape,'(a8)') "<PP_RAB>"
      write(iotape,'(4(1x,e24.15,1x))') rabupf
    write(iotape,'(a9)') "</PP_RAB>"
  write(iotape,'(a10)') "</PP_MESH>"


  if(nicore /= 'nc  ') then
    write(iotape,'(a22)') "<PP_NLCC columns='4'>"
      write(iotape,'(4(1x,e24.15,1x))') cdcupf
    write(iotape,'(a10)') "</PP_NLCC>"
  endif


  write(iotape,'(a24)') "<PP_LOCAL columns='4'>"
    write(iotape,'(4(1x,e24.15,1x))') vlocupf
  write(iotape,'(a11)') "</PP_LOCAL>"

  write(iotape,'(a13)') "<PP_NONLOCAL>"
    np = 0
    do l = 0,lmax_pot
      if(l /= llocal) then
        np = np+1
        write(iotape,'(a9,i1)') "<PP_BETA.",np
        write(iotape,'("  index=''",i1,"''")') np
        ctemp = psdtitle(l*2+1)
        nl = '  '
        read(ctemp(1:2),'(a2)',iostat=ierr) nl
        write(iotape,'("  label=''",a2,"''")') nl
        write(iotape,'("  angular_momentum=''",i1,"''")') l
        write(iotape,'("  cutoff_radius_index=''",i6,"''")') nrupf
        write(iotape,'("  cutoff_radius=''",e22.15,"''")') ZERO
        write(iotape,'("  ultrasoft_cutoff_radius=''",e22.15,"''")') ZERO
        write(iotape,'(a1)') ">"
        write(iotape,'(4(1x,e24.15,1x))') vnlupf(:,l)
        write(iotape,'(a10,i1,a1)') "</PP_BETA.",np,">"
      endif
    enddo
    if(np /= nproj) then
      write(6,*) 'inconsistency in number of projectors'
      stop
    endif
    write(iotape,'(a8)') "<PP_DIJ>"
      write(iotape,'(4(1x,e24.15,1x))') ((dij(i,j),i=1,nproj),j=1,nproj)
    write(iotape,'(a9)') "</PP_DIJ>"
  write(iotape,'(a14)') "</PP_NONLOCAL>"


  write(iotape,'(a10)') "<PP_PSWFC>"
    do l = 0,lmax_pot
      write(iotape,'(a9,i1)') "<PP_CHI.",l+1
      write(iotape,'("  index=''",i1,"''")') l+1
      ctemp = psdtitle(l*2+1)
      nl = '  '
      read(ctemp(1:2),'(a2)',iostat=ierr) nl
      write(iotape,'("  label=''",a2,"''")') nl
      write(iotape,'("  l=''",i1,"''")') l
      occup = ZERO
      read(ctemp(4:9),*,iostat=ierr) occup
      if(ierr == 0) then
        write(iotape,'("  occupation=''",e22.15,"''")') occup
      endif
      write(iotape,'("  pseudo_energy=''",e22.15,"''")') ZERO
      write(iotape,'("  cutoff_radius=''",e22.15,"''")') ZERO
      write(iotape,'("  ultrasoft_cutoff_radius=''",e22.15,"''")') ZERO
      write(iotape,'(a1)') ">"
        write(iotape,'(4(1x,e24.15,1x))') chi(:,l)
      write(iotape,'(a9,i1,a1)') "</PP_CHI.",l+1,">"
    enddo
  write(iotape,'(a11)') "</PP_PSWFC>"


  write(iotape,'(a12)') "<PP_RHOATOM>"
    write(iotape,'(4(1x,e24.15,1x))') cdvupf
  write(iotape,'(a13)') "</PP_RHOATOM>"


  write(iotape,'(a6)') "</UPF>"


  close(unit=iotape)

  deallocate(rupf,rabupf)

  deallocate(vlocupf)
  deallocate(vnlupf)
  deallocate(dij)
  deallocate(chi)
  deallocate(cdcupf)
  deallocate(cdvupf)

  deallocate(yp,ypp,wlu)
  deallocate(ypi,yppi)

  return

end subroutine atom_kb_psd_out_upf
