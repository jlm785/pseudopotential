!>  main subroutine to generate a pseudopotential
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.0.6
!>  \date         April 1990, 30 June 2021, 3 November 2021.
!>  \copyright    GNU Public License v2

subroutine atom_psd_sub(rc, ifcore, cfac, rcfac, lrci,                   &
         iowrite, ioae, fileae,                                          &
         iopsd, filepsd, ioparsec, fileparsec, ioplot, fileplot)


! Written 12 April 2018
! Modified, July 2021. JLM
! so->iso, vionic, vhxc, cdpsd, vscr. 15, 23 September 2021. JLM
! printing. 21 October 2021. JLM
! indx, rpsi, 3 November 2021. JLM
! lrci, 4 November 2021. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter                ::  mxdlc = 4                        !  hard coded maximum valence angular momentum
  integer, parameter                ::  mxdnw = 2                        !  hard coded dimension of number of wave-functions same l

! input

  integer, intent(in)               ::  iowrite                          !<  default output tape

  integer, intent(in)               ::  ioae                             !<  default tape for writing
  character(len=*), intent(in)      ::  fileae                           !<  name of default tape for all-electron results

  integer, intent(in)               ::  iopsd                            !<  default tape for pseudopotential in old format
  character(len=*), intent(in)      ::  filepsd                          !<  name of default tape for reading pseudopotential in old format

  integer, intent(in)               ::  ioparsec                         !<  default tape for pseudopotential in parsec format
  character(len=*)                  ::  fileparsec                       !<  name of default tape for reading pseudopotential in parsec format

  integer, intent(in)               ::  ioplot                           !<  default tape for plot file
  character(len=*), intent(in)      ::  fileplot                         !<  name of default tape for plot file

  real(REAL64), intent(inout)       ::  rc(0:mxdlc)                      !<  core radius r_c(l).  Rounded up.
  integer , intent(in)              ::  ifcore                           !<  0: no partial core correction; 1 : xc partial core correction; 2 : not recommended.
  real(REAL64), intent(in)          ::  cfac                             !<  criteria for pseudo-core charge
  real(REAL64), intent(in)          ::  rcfac                            !<  pseudo-core radius

  logical, intent(in)               ::  lrci                             !<  uses rverse communication interface

! dimensions

  integer                           ::  mxdnr                            !  dimension of the number of radial points
  integer                           ::  mxdorb                           !  dimension of the number of orbitals
  integer                           ::  mxdl                             !  dimension of angular momentum


! variables

  integer                           ::  nr                               !  number of radial points

  integer                           ::  norb                             !  number of orbitals
  integer                           ::  lmax                             !  maximum angular momentum

  integer                           ::  itype                            !  unused
  character(len=2)                  ::  icorr                            !  correlation type
  character(len=1)                  ::  ispp                             !  spin polarization  ' ', 's', 'r'

  real(REAL64)                      ::  a                                !  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64)                      ::  b                                !  r(i) = a*(exp(b*(i-1))-1)

  character(len=2)                  ::  nameat                           !  chemical symbot of the atom

  integer                           ::  ncore                            !  number of orbitals treated as core

  real(REAL64)                      ::  znuc                             !  nuclear charge
  real(REAL64)                      ::  zel                              !  electron charge

   real(REAL64)                      ::  zion, zratio
  character(len = 5)                ::  vers

  integer              ::  npot(-1:1)

  character(len=3)     ::  irel
  character(len=4)     ::  nicore
  character(len=10)    ::  iray(6), ititle(7)

  integer              ::  isplot

  integer              ::  lcmax

! allocatable arrays

  real(REAL64), allocatable         ::  r(:)                             !  radial grid points
  real(REAL64), allocatable         ::  drdi(:)                          !  d r(i) / d i
  real(REAL64), allocatable         ::  d2rodr(:)                        !  (d^2 r(i) / d i^2) /  (d r(i) / d i)

  integer, allocatable              ::  no(:)                            !  principal quantum number n
  integer, allocatable              ::  lo(:)                            !  angular quantum number l
  integer, allocatable              ::  iso(:)                           !  2*spin or 2*(j-l)
  real(REAL64), allocatable         ::  zo(:)                            !  orbital occupation

  real(REAL64), allocatable         ::  cdc(:)                           !  core charge density (total)

  real(REAL64), allocatable         ::  vionic(:,:,:)                    !  r*ionic potential in Rydberg

  real(REAL64), allocatable         ::  vhxc(:,:)                        !  effective potential in Rydberg

  real(REAL64), allocatable         ::  vscr(:,:,:)                      !  screened-pseudo-potential in Rydberg (down or total)

  real(REAL64), allocatable         ::  vpsd(:,:,:)                      !  pseudo-potential in Rydberg

  real(REAL64), allocatable         ::  cdpsd(:,:)                       !  pseudo-charge density

  real(REAL64), allocatable         ::  ev(:)                            !  orbital energy
  real(REAL64), allocatable         ::  evi(:)                           !  orbital energy
  real(REAL64), allocatable         ::  ev_ae(:)                         !  orbital energy

  real(REAL64), allocatable         ::  rpsi_ps(:,:,:,:)                 !  r*pseudo-wave-function
  real(REAL64), allocatable         ::  rpsi_ae(:,:,:,:)                 !  r*all-electron-wave-function (major component if relativistic)
  real(REAL64), allocatable         ::  br_ae(:,:,:,:)                   !  d rpsi_ae / dr or minor component

  real(REAL64), allocatable         ::  v(:)

  integer, allocatable              ::  indv(:,:)                        !  main orbital for scattering channel l (legacy)
  integer, allocatable              ::  indx(:,:,:)                      !  orbitals for scattering channel l
  integer, allocatable              ::  nindx(:,:)                       !  number of orbitals for scattering channel l

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64

! counter

  integer     ::  is, i, j, k, l


  write(6,*)
  write(6,*) '  Constructing the pseudopotential'
  write(6,*)

! get array sizes

  call atom_psd_dat_in_size(ioae, fileae, mxdl, mxdnr, mxdorb)

! allocates arrays

  allocate(r(mxdnr))
  allocate(drdi(mxdnr))
  allocate(d2rodr(mxdnr))

  allocate(cdpsd(mxdnr,-1:1))

  allocate(no(mxdorb))
  allocate(lo(mxdorb))
  allocate(iso(mxdorb))
  allocate(zo(mxdorb))

  allocate(cdc(mxdnr))
  allocate(vionic(mxdnr,0:mxdl,-1:1))
  allocate(vhxc(mxdnr,-1:1))

  allocate(ev(mxdorb))
  allocate(evi(mxdorb))
  allocate(ev_ae(mxdorb))

  allocate(rpsi_ps(mxdnr,mxdnw,0:mxdlc,-1:1))
  allocate(rpsi_ae(mxdnr,mxdnw,0:mxdlc,-1:1))
  allocate(br_ae(mxdnr,mxdnw,0:mxdlc,-1:1))

  allocate(v(mxdnr))

  allocate(vscr(mxdnr,0:mxdlc,-1:1))

  allocate(vpsd(mxdnr,0:mxdlc,-1:1))

! reads atomic data

  call atom_psd_dat_in(ioae, fileae, itype, icorr, ispp, lmax,           &
      nr, a, b, r, drdi, d2rodr,                                         &
      nameat, norb, ncore, no, lo, iso, zo, znuc, zel,                   &
      cdc, vionic, vhxc, ev,                                             &
      mxdnr, mxdorb, mxdl)

  do j = 1,norb
    ev_ae(j) = ev(j)
    evi(j) = ZERO
  enddo


! fills indv, indx

  allocate(indv(0:mxdlc,-1:1))
  allocate(indx(mxdnw,0:mxdlc,-1:1))
  allocate(nindx(0:mxdlc,-1:1))

  call atom_psd_ind(indv, indx, nindx, lcmax, norb, ncore, lo, iso,  &
        mxdorb, mxdlc, mxdnw)

  if(ispp == 'r' .or. ispp == 's') then
    isplot = 1
  else
    isplot = 0
  endif

! reconstructs all-electron functions


  do is = -1,1
  do l = 0,lcmax
    if(nindx(l,is) > 0) then
      do j = 1,nindx(l,is)

        i = indx(j,l,is)
        call atom_psd_ae_orb(ispp, nr, r, drdi, d2rodr, v,               &
            rpsi_ae(:,j,l,is), br_ae(:,j,l,is),                          &
            no(i), l, is, znuc, ev(i), .FALSE.,                          &
            vionic(:,l,is), vhxc(:,is),                                  &
            iowrite, mxdnr)

      enddo

    endif
  enddo
  enddo


! generates screened pseudopotential

! print heading

  write(iowrite,*)
  write(iowrite,*)
  write(iowrite,'(1x,a2," pseudopotential generation using the ",        &
       &  "Improved Troullier and Martins method")') nameat
  write(iowrite,'(1x,60("-"))')
  write(iowrite,*)
  write(iowrite,'(" nl    s    eigenvalue",6x,"rc",10x,"cdrc",7x,        &
        &   "delta",7x,"gamma",7x,"alpha_j")')
  write(iowrite,*)

  do is = -1,1
  do l = 0,lcmax
    if(nindx(l,is) == 1) then

      i = indx(1,l,is)
      call atom_psd_tm2(ispp, nr, r, drdi, lrci,                         &
          no(i), l, is, rc(l),                                           &
          vionic(:,l,is), vhxc(:,is), ev(i),                             &
          rpsi_ae(:,1,l,is), br_ae(:,1,l,is),                            &
          rpsi_ps(:,1,l,is), vscr(:,l,is),                               &
          iowrite, mxdnr)

    elseif(nindx(l,is) == 2) then

      call atom_psd_mrpp(ispp, nr, r, drdi, d2rodr,                        &
          no(indx(1,l,is)), l, is, rc(l),                                  &
          vionic(:,l,is), vhxc(:,is), ev(indx(1,l,is)), ev(indx(2,l,is)),  &
          rpsi_ae(:,1:2,l,is), br_ae(:,1:2,l,is),                          &
          rpsi_ps(:,1:2,l,is), vscr(:,l,is),                               &
          iowrite, mxdnr)

    elseif(nindx(l,is) > 2) then

      write(6,*)
      write(6,'("   Stopped in atom_psd_sub:   for l,is = ",2i3)') l,is
      write(6,'("   the number of valence orbitals ",i5," > 2")')  nindx(l,is)
      write(6,*)

      STOP

    endif
  enddo
  enddo


! pseudo-charge density

  cdpsd = ZERO

  do is = -1,1
  do l = 0,lcmax
    if(nindx(l,is) > 0) then
      do k = 1,nindx(l,is)
        i = indx(k,l,is)
        do j = 1,nr
          cdpsd(j,is) = cdpsd(j,is) + zo(i)*rpsi_ps(j,k,l,is)*rpsi_ps(j,k,l,is)
        enddo
      enddo
    endif
  enddo
  enddo

! recasts cdpsd as expected by atom_atm_velect

  if(ispp == ' ') then
    do j = 1,nr
      cdpsd(j, 1) = cdpsd(j, 0) / 2
      cdpsd(j,-1) = cdpsd(j, 0) / 2
    enddo
  else
    do j = 1,nr
      cdpsd(j, 0) = cdpsd(j, 1) + cdpsd(j,-1)
    enddo
  endif


! unscreens pseudopotential

  call atom_psd_unscreen(ifcore, icorr, ispp, nr, r, drdi,               &
      norb, ncore, lo, iso, zo, znuc, zel,                               &
      cdpsd, cdc, vscr, vpsd,                                            &
      cfac, rcfac, zratio, zion,                                         &
      iowrite, mxdnr, mxdorb)


! writes the plotting data file

  call atom_psd_plot(ioplot, fileplot, nr, r, drdi,                      &
      nindx, lcmax, rpsi_ps, rpsi_ae, vpsd, zion,                        &
      mxdnr, mxdlc, mxdnw)


! performs a simple test

  call atom_psd_test_simple(ifcore, icorr, ispp, rc,                     &
      nr, r, drdi, d2rodr,                                               &
      nameat, norb, ncore, lo, iso, zo, znuc, zel,                       &
      cdpsd, cdc, vpsd,                                                  &
      ev, evi, indv, zratio, zion,                                       &
      iowrite, mxdnr, mxdorb)


! writes files

  call atom_version(vers)

  call atom_psd_titles('tm2', vers, indv, no, zo,                        &
      ispp, ifcore, rc, cfac, zratio, ncore, norb,                       &
      iray, ititle, nicore, irel, npot,                                  &
      mxdorb)

  call atom_psd_out_unfmt(iopsd, filepsd, nameat,                        &
      icorr, irel, nicore, iray, ititle,                                 &
      npot, nr, a, b, r, zion, indv, ifcore,                             &
      vpsd, cdc, cdpsd,                                                  &
      mxdnr)

  call atom_psd_out_parsec(ioparsec, fileparsec, nameat,                 &
      icorr, irel, nicore, iray, ititle,                                 &
      npot, nr, a, b, r, zion, indv, ifcore,                             &
      vpsd, cdc, cdpsd,                                                  &
      cfac, rcfac, ncore+1, norb, lo, rc, zo, rpsi_ps,                   &
      mxdnr, mxdorb, mxdlc, mxdnw)

! deallocates arrays

  deallocate(r)
  deallocate(drdi)
  deallocate(d2rodr)

  deallocate(cdpsd)

  deallocate(no)
  deallocate(lo)
  deallocate(iso)
  deallocate(zo)

  deallocate(cdc)
  deallocate(vionic)
  deallocate(vhxc)

  deallocate(ev)
  deallocate(evi)
  deallocate(ev_ae)

  deallocate(rpsi_ps)
  deallocate(rpsi_ae)
  deallocate(br_ae)

  deallocate(v)
  deallocate(vscr)

  deallocate(vpsd)

  deallocate(indv)
  deallocate(indx)
  deallocate(nindx)

  return

end subroutine atom_psd_sub
