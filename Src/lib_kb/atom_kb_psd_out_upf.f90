!>  writes the KB pseudopotential in the Quantum Espresso UPF format
!>
!>  \author       J.L.Martins
!>  \version      6.0.9
!>  \date         October 2018. 17 November 2024.
!>  \copyright    GNU Public License v2

subroutine atom_kb_psd_out_upf(iotape, fname,                            &
      nameat, icorr, irel, nicore, irdate, irvers, irayps, psdtitle,     &
      nr, r,  zion, vlocal, cdc, cdv,                                    &
      llocal, lmax_pot, vkbproj, inorm, lmax_psi, rpsi,                  &
      mxdl, mxdnr)

! Adapted by JLMartins October 2018
! Converted to new style. 25 January 2022. JLM
! psdtitle, bug/feature reported by Raymond Atta-Fynn, 19 May 2022. JLM
! Valid now for pseudopotentials wit spin-orbit. 16 November 2024.

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
  integer                 ::  nrupf                                      !  number of points in the UPF grid
  real(REAL64)            ::  xmin
  real(REAL64)            ::  rmax, a, b

  real(REAL64)            ::  occup, rc

  integer                 ::  nproj, nchi

  character(len=10)       ::  ctemp

  logical                 ::  lso                                        !  indicates if spin-orbit is used (irel == 'rel')

  integer                 ::  ierr

! allocatable local variables

  real(REAL64), allocatable        ::  rupf(:)                           !  radial grid points   rupf(i) = a*(exp(b*(i-1))), i=1,...,nrupf
  real(REAL64), allocatable        ::  rabupf(:)                         !  d rupf(i) / d i
  real(REAL64), allocatable        ::  vlocupf(:)                        !  local pseudopotential (not r*vlocal) for UPF
  real(REAL64), allocatable        ::  vnlupf(:,:,:)                     !  kb-projectors (multiplied by r) for UPF
  real(REAL64), allocatable        ::  dij(:,:)                          !  d_ij = inorm(i)*delta_ij
  real(REAL64), allocatable        ::  chi(:,:,:)                        !  wavefunctions (r(i),l,2j-2l) for UPF (r*psi)
  real(REAL64), allocatable        ::  cdcupf(:)                         !  charge density of core.
  real(REAL64), allocatable        ::  cdvupf(:)                         !  4*pi*r**2 charge density of valence.

  integer, allocatable             ::  l_vnl(:)                          !  angular momentum for projector i
  integer, allocatable             ::  is_vnl(:)                         !  2(j-l) for projector i

  integer, allocatable             ::  l_chi(:)                          !  angular momentum for wave-function i
  integer, allocatable             ::  is_chi(:)                         !  2(j-l) for wave-function i

  integer, allocatable             ::  n_conf(:)                         !  principal quantum number
  logical, allocatable             ::  l_skip(:)                         !  skip as information is missing
  real(REAL64), allocatable        ::  o_conf(:,:)                       !  occupation of orbital
  real(REAL64), allocatable        ::  r_core(:)                         !  core radius

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  character(len=1), parameter  ::  IL(0:6) = (/'s','p','d','f','g','h','i'/)

! counters

  integer       ::   i, j, l, np, is


  if(iotape == 0) return

  lso = .FALSE.
  if(irel == 'rel') lso = .TRUE.

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

! interpolates on the upf mesh the local potential and non-local pseudopotentials,
! the wavefunctions and charge densities (core and valence)

  allocate(vlocupf(nrupf))
  allocate(vnlupf(nrupf,0:mxdl,-1:1))
  allocate(chi(nrupf,0:mxdl,-1:1))
  allocate(cdcupf(nrupf))
  allocate(cdvupf(nrupf))

  call atom_kb_psd_out_upf_interp(nr, r, nrupf, rupf, lso,               &
      vlocal, lmax_pot, vkbproj, lmax_psi, rpsi, cdv, cdc,               &
      vlocupf, vnlupf, chi, cdvupf, cdcupf,                              &
      mxdl, mxdnr)

! gets the configuration as UPF needs

  allocate(n_conf(0:lmax_pot))
  allocate(l_skip(0:lmax_pot))
  allocate(o_conf(0:lmax_pot,2))
  allocate(r_core(0:lmax_pot))

  allocate(l_vnl(2*lmax_pot+1))
  allocate(is_vnl(2*lmax_pot+1))

  allocate(l_chi(2*lmax_pot+1))
  allocate(is_chi(2*lmax_pot+1))

  call atom_kb_psd_out_upf_config(lmax_pot, psdtitle, lso,               &
      llocal, inorm,                                                     &
      n_conf, l_skip, o_conf, r_core,                                    &
      nproj, l_vnl, is_vnl, nchi, l_chi, is_chi,                         &
      mxdl)


! The coding of inorm

  allocate(dij(nproj,nproj))

  do i = 1,nproj
    do j = 1,nproj
      dij(i,j) = ZERO
    enddo
  enddo
  do np = 1,nproj
    dij(np,np) = UM * inorm(l_vnl(np),is_vnl(np))
  enddo


! Writes the file

!  filename = adjustl(trim(nameat))//"_LDA_ncpp.UPF"

! nchi and nproj are written with 'i1'

  if(nchi < 10 .and. nproj < 10) then

    open (unit = iotape, file = trim(fname), form = 'formatted', status = 'unknown')

    write(iotape,'(a21)') "<UPF version='2.0.1'>"
      write(iotape,*)



      write(iotape,'(a9)') "<PP_INFO>"
        write(iotape,'("    Generated with atom program version:  ",4x,a10)') irvers
        write(iotape,'("    https://github.com/jlm785/pseudopotential")')
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

! ! !   the use of psdtitle is an hack.  Should be done right...

        if(llocal >= 0) then
          ctemp = psdtitle(2*llocal+2)
          rc = ZERO
          read(ctemp(6:10),*,iostat=ierr) rc
          write(iotape,'("    L of local component and cutoff radius",      &
             &     5x,i3,f11.5)') llocal, rc
        else
           write(iotape,'("    Local component is maximum of all",          &
             &    " potentials (smoothed)")')
        endif

        write(iotape,'("    Valence configuration:")')
        write(iotape,'("    nl",9x,"occ",12x,"Rcut")')
        do l = 0,lmax_pot
          if(.not. l_skip(l)) then
            write(iotape,'(4x,i1,a1,5x,f10.4,5x,f10.4)') n_conf(l),IL(l),  &
                 o_conf(l,1)+o_conf(l,2),r_core(l)
          endif
        enddo
        write(iotape,'("    UPF generation configuration: not available")')

      write(iotape,'(a10)') "</PP_INFO>"
      write(iotape,*)



      write(iotape,'(a4)') "<!--"
      write(iotape,'(a29)') "END OF HUMAN READABLE SECTION"
      write(iotape,'(a4)') "-->"
      write(iotape,*)



      write(iotape,'(a10)') "<PP_HEADER"
        write(iotape,'("  generated=''Generated by atom version ",a10,"''")') irvers
        write(iotape,'("  author=''Froyen,Troullier,Martins et al.''")')
        write(iotape,'("  date=''",a10,"''")')  irdate
        write(iotape,'("  comment=''using beta code, beware''")')
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
        if(lso) then
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
        write(iotape,'("  number_of_wfc=''",i2,"''")') nchi
        write(iotape,'("  number_of_proj=''",i2,"''")') nproj
      write(iotape,'(a2)') "/>"
      write(iotape,*)



      write(iotape,'(a8)') "<PP_MESH"
        write(iotape,'("  dx=''",e22.15,"''")') b
        write(iotape,'("  mesh=''",i6,"''")') nrupf
        write(iotape,'("  xmin=''",e22.15,"''")') xmin
        write(iotape,'("  rmax=''",e22.15,"''")') rupf(nrupf)
        write(iotape,'("  zmesh=''",e22.15,"''")') iz*UM
        write(iotape,'(a1)') ">"
        write(iotape,'(a18)') "<PP_R columns='4'>"

          write(iotape,'(4(1x,e24.15,1x))') rupf

        write(iotape,'(a7)') "</PP_R>"
        write(iotape,'(a20)') "<PP_RAB columns='4'>"

          write(iotape,'(4(1x,e24.15,1x))') rabupf

        write(iotape,'(a9)') "</PP_RAB>"
      write(iotape,'(a10)') "</PP_MESH>"
      write(iotape,*)



      if(nicore /= 'nc  ') then
        write(iotape,'(a22)') "<PP_NLCC columns='4'>"
          write(iotape,'(4(1x,e24.15,1x))') cdcupf
        write(iotape,'(a10)') "</PP_NLCC>"
      endif
      write(iotape,*)


      write(iotape,'(a24)') "<PP_LOCAL columns='4'>"
        write(iotape,'(4(1x,e24.15,1x))') vlocupf
      write(iotape,'(a11)') "</PP_LOCAL>"
      write(iotape,*)



      if(nproj > 0) then
        write(iotape,'(a13)') "<PP_NONLOCAL>"

          do np = 1,nproj
            l = l_vnl(np)
            is = is_vnl(np)

            write(iotape,'(a9,i1)') "<PP_BETA.",np
              write(iotape,'("  index=''",i1,"''")') np
              write(iotape,'("  label=''",i1,a1,"''")') n_conf(l), IL(l)
              write(iotape,'("  angular_momentum=''",i1,"''")') l
              write(iotape,'("  cutoff_radius_index=''",i6,"''")') nrupf
              write(iotape,'("  cutoff_radius=''",e22.15,"''")') r_core(l)
              write(iotape,'("  ultrasoft_cutoff_radius=''",e22.15,"''")') ZERO
              write(iotape,'(a1)') ">"

              write(iotape,'(4(1x,e24.15,1x))') vnlupf(:,l,is)

            write(iotape,'(a10,i1,a1)') "</PP_BETA.",np,">"

          enddo

          write(iotape,'(a8)') "<PP_DIJ>"
            write(iotape,'(4(1x,e24.15,1x))') ((dij(i,j),i=1,nproj),j=1,nproj)
          write(iotape,'(a9)') "</PP_DIJ>"

        write(iotape,'(a14)') "</PP_NONLOCAL>"
        write(iotape,*)
      endif



      write(iotape,'(a10)') "<PP_PSWFC>"
        do np = 1,nchi
          l = l_chi(np)
          is = is_chi(np)
          write(iotape,'(a9,i1)') "<PP_CHI.",np

            write(iotape,'("  index=''",i1,"''")') np
            write(iotape,'("  label=''",i1,a1,"''")') n_conf(l), IL(l)
            write(iotape,'("  l=''",i1,"''")') l
            occup = o_conf(l,1)+o_conf(l,2)
            if(is == 1) then
              occup = occup*(l+1) / (UM*(2*l+1))
            elseif(is == -1) then
              occup = occup*l / (UM*(2*l+1))
            endif
            write(iotape,'("  occupation=''",e22.15,"''")') occup
            write(iotape,'("  pseudo_energy=''",e22.15,"''")') ZERO
            write(iotape,'("  cutoff_radius=''",e22.15,"''")') ZERO
            write(iotape,'("  ultrasoft_cutoff_radius=''",e22.15,"''")') ZERO
            write(iotape,'(a1)') ">"

              write(iotape,'(4(1x,e24.15,1x))') chi(:,l,is)

          write(iotape,'(a9,i1,a1)') "</PP_CHI.",np,">"
        enddo
      write(iotape,'(a11)') "</PP_PSWFC>"
      write(iotape,*)



      write(iotape,'(a12)') "<PP_RHOATOM>"
        write(iotape,'(4(1x,e24.15,1x))') cdvupf
      write(iotape,'(a13)') "</PP_RHOATOM>"
      write(iotape,*)



      if(lso) then
        write(iotape,'(a13)') "<PP_SPIN_ORB>"
          if(nproj > 0) then
            do np = 1,nproj
              write(iotape,'("  <PP_RELBETA.",i1,"   index=''",i1,         &
                 &  "''     lll=''",i1,"''     jjj=''",f3.1,"''/>")')      &
                    np,np,l_vnl(np),(UM*(2*l_vnl(np)+is_vnl(np)))/2
            enddo
          endif
          do np = 1,nchi
            write(iotape,'("  <PP_RELWFC.",i1,"    index=''",i1,           &
               &  "''     lchi=''",i1,"''    jchi=''",f3.1,"''/>")')       &
                  np,np,l_chi(np),(UM*(2*l_chi(np)+is_chi(np)))/2
          enddo
        write(iotape,'(a14)') "</PP_SPIN_ORB>"
        write(iotape,*)
      endif



    write(iotape,'(a6)') "</UPF>"


    close(unit=iotape)

  else
    write(6,*)
    write(6,*) '  file ',trim(fname),' not written'
    write(6,*) '  rewrite a couple of lines in atom_kb_psd_out_upf'
    write(6,*) '  nproj = ',nproj,'  nchi = ',nchi
    write(6,*)
  endif

  deallocate(rupf,rabupf)

  deallocate(vlocupf)
  deallocate(vnlupf)
  deallocate(dij)
  deallocate(chi)
  deallocate(cdcupf)
  deallocate(cdvupf)

  deallocate(l_vnl)
  deallocate(is_vnl)
  deallocate(l_chi)
  deallocate(is_chi)

  deallocate(n_conf)
  deallocate(l_skip)
  deallocate(o_conf)
  deallocate(r_core)

  return

end subroutine atom_kb_psd_out_upf
