!>  Performs a series of tests of the pseudopotentials
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.4
!>  \date         October 2021
!>  \copyright    GNU Public License v2

subroutine atom_atm_test(iowrite, ioread, filein,                        &
         iopsd, filepsd, iopsdkb, filepsdkb)

! translated to f90 from atm.f version 5.805
! configurations without self-consistent solution are skipped (ifail). 8 November 2026. claude+JLM

! copyright  J.L.Martins, INESC-MN.


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  ioread                           !<  default tape for reading
  character(len=*), intent(in)      ::  filein                           !<  name of default tape for reading

  integer, intent(in)               ::  iopsd                            !<  default tape for pseudopotential in parsec format
  character(len=*), intent(in)      ::  filepsd                          !<  name of default tape for reading pseudopotential in parsec format

  integer, intent(in)               ::  iopsdkb                          !<  default tape for KB pseudopotential in real space
  character(len=15), intent(in)     ::  filepsdkb                        !<  name of default tape for writeing KB pseudopotential in real space

! array dimensions

  integer                           ::  mxdorb                           !  dimension of the number of orbitals
  integer                           ::  mxdnr                            !  dimension of the number of radial points
  integer                           ::  mxdl                             !  dimension for angular momentum l

! local variables

  character(len=12)                 ::  fileae                           !  dummy filename

  character(len=2)                  ::  nameat                           !  chemical symbol of the atom

  character(len=2)                  ::  ctype                            !  type of calculation flag (converted to itype)
  character(len=10)                 ::  ititle(5)                        !  title of calculation

  character(len=3)                  ::  kerker                           !  type of pseudopotential flag

  character(len=2)                  ::  icorr                            !  correlation type
  character(len=1)                  ::  ispp                             !  spin polarization  ' ', 's', 'r'

  real(REAL64)                      ::  znuc                             !  nuclear charge

  real(REAL64)                      ::  rmax                             !  maximum mesh radius
  real(REAL64)                      ::  aa                               !  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64)                      ::  bb                               !  a = exp(-aa)/znuc, b = 1/bb

  integer                           ::  ncore                            !  number of orbitals treated as core
  integer                           ::  nval                             !  number of valence orbitals

  real(REAL64)                      ::  etotal                           !  total energy
  real(REAL64)                      ::  eae_ref                          !  all electron reference energy
  real(REAL64)                      ::  epsd_ref                         !  pseudopotential reference energy
  real(REAL64)                      ::  ekb_ref                          !  KB projector reference energy

  real(REAL64)                      ::  zui, zdi                         !  amount of ionization/excitation
  real(REAL64)                      ::  zux, zdx                         !  amount of ionization/excitation

  real(REAL64)                      ::  ztot                             !  valence charge
  real(REAL64)                      ::  zz                               !  excitation charge

  integer                           ::  nx                               !  number of excitations calculated

  integer, allocatable              ::  ni(:)                            !  valence principal quantum number
  integer, allocatable              ::  li(:)                            !  valence angular quantum number
  real(REAL64), allocatable         ::  zd(:), zu(:)                     !  occupation of valence down and up orbitals
  real(REAL64), allocatable         ::  evd(:)                           !  default eigenvalue

  real(REAL64), allocatable         ::  eion_ae(:)                       !  all-electron ionization
  real(REAL64), allocatable         ::  eion_psd(:)                      !  pseudopotential ionization
  real(REAL64), allocatable         ::  eion_kb(:)                       !  KB projector ionization

  real(REAL64), allocatable         ::  exct_ae(:)                       !  all-electron excitation
  real(REAL64), allocatable         ::  exct_psd(:)                      !  pseudopotential excitation
  real(REAL64), allocatable         ::  exct_kb(:)                       !  KB projector excitation

  logical, allocatable              ::  lempty(:), lfull(:)              !  indicates mostly empty or full orbital

  integer                           ::  ifail                            !  /= 0: no self-consistent solution
  character(len=6)                  ::  fail_ref                         !  calculation (ae, pseudo, KB) that failed for the reference
  character(len=6), allocatable     ::  fail_ion(:)                      !  calculation that failed for an ionization (blank if none)
  character(len=6), allocatable     ::  fail_exct(:)                     !  calculation that failed for an excitation (blank if none)
  real(REAL64), allocatable         ::  zexct(:), zexct2(:)              !  excitation charge

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  RYDBERG = 27.21138386_REAL64 / 2
  character(len=1), parameter  ::  IL(0:6) = (/'s','p','d','f','g','h','i'/)
! counter

  integer           ::  n, n2


  fileae = 'dummy.dat   '

! finds dimensions

  open(unit=ioread, file=trim(filein), status='OLD', form='FORMATTED')

  call atom_atm_input_size(ioread, iopsd, filepsd, mxdnr, mxdorb, mxdl)

  close(unit=ioread)

! gets configuration

  open(unit=ioread, file=trim(filein), status='OLD', form='FORMATTED')

  allocate(ni(mxdorb))
  allocate(li(mxdorb))
  allocate(zu(mxdorb),zd(mxdorb))
  allocate(evd(mxdorb))

  call atom_atm_read_input(nameat, ctype, ititle, kerker, icorr, ispp,   &
      znuc, rmax, aa, bb, ncore, nval, ni, li, zd, zu, evd,              &
      ioread, mxdorb)

  close(unit=ioread)

  allocate(eion_ae(nval))
  allocate(eion_psd(nval))
  allocate(eion_kb(nval))

  allocate(exct_ae(nval*nval))
  allocate(exct_psd(nval*nval))
  allocate(exct_kb(nval*nval))
  allocate(zexct(nval*nval))
  allocate(zexct2(nval*nval))

  allocate(lempty(nval))
  allocate(lfull(nval))

  allocate(fail_ion(nval))
  allocate(fail_exct(nval*nval))
  fail_ion(:) = '      '
  fail_exct(:) = '      '

! calculates the all-electron reference energy

  fail_ref = '      '

  ctype = 'ae'

  call atom_atm_scf(eae_ref, ifail,                                      &
      nameat, ctype, ititle, kerker, icorr, ispp,                        &
      znuc, rmax, aa, bb, ncore, nval, ni, li, zd, zu, evd,              &
      iowrite, 0, fileae, iopsd, filepsd,                                &
      mxdnr, mxdorb, mxdl)

  if(ifail /= 0) fail_ref = 'ae'

! calculates the pseudopotential reference energy

  if(fail_ref == ' ') then

    ctype = 'pt'

    call atom_atm_scf(epsd_ref, ifail,                                   &
        nameat, ctype, ititle, kerker, icorr, ispp,                      &
        znuc, rmax, aa, bb, ncore, nval, ni, li, zd, zu, evd,            &
        iowrite, 0, fileae, iopsd, filepsd,                              &
        mxdnr, mxdorb, mxdl)

    if(ifail /= 0) fail_ref = 'pseudo'

  endif

  if(fail_ref == ' ') then

    call atom_kb_test_scf(ekb_ref, ifail,                                &
           nameat, ititle, icorr, ispp,                                  &
           nval, ni, li, zd, zu, evd,                                    &
           iowrite, iopsdkb, filepsdkb,                                  &
           mxdnr, mxdorb, mxdl)

    if(ifail /= 0) fail_ref = 'KB'

  endif

! without the reference energies there is nothing to compare

  if(fail_ref /= ' ') then

    write(6,*)
    write(6,'("    Tests skipped, the ",a," calculation of the reference ",   &
       &     "configuration did not converge")') trim(fail_ref)
    write(6,*)
    if(iowrite /= 6) then
      write(iowrite,*)
      write(iowrite,'("    Tests skipped, the ",a," calculation of the ",   &
         &     "reference configuration did not converge")') trim(fail_ref)
      write(iowrite,*)
    endif

    deallocate(ni)
    deallocate(li)
    deallocate(zu,zd)
    deallocate(evd)

    deallocate(eion_ae)
    deallocate(eion_psd)
    deallocate(eion_kb)

    deallocate(exct_ae)
    deallocate(exct_psd)
    deallocate(exct_kb)

    deallocate(lempty)
    deallocate(lfull)
    deallocate(zexct)
    deallocate(zexct2)

    deallocate(fail_ion)
    deallocate(fail_exct)

    return

  endif

  ztot = ZERO
  do n = 1,nval
    ztot = ztot + zd(n) + zu(n)
  enddo

! loop over valence states for ionization energies

  do n = 1,nval

    zui = min(zu(n),ONE)
    zdi = min(zd(n),ONE-zui)

!   minimum of ionization to make it worth

    if(zui + zdi > ONE/10) then

      zd(n) = zd(n) - zdi
      zu(n) = zu(n) - zui

!     a calculation without self-consistent solution skips the configuration

      ctype = 'ae'

      call atom_atm_scf(etotal, ifail,                                   &
          nameat, ctype, ititle, kerker, icorr, ispp,                    &
          znuc, rmax, aa, bb, ncore, nval, ni, li, zd, zu, evd,          &
          iowrite, 0, fileae, iopsd, filepsd,                            &
          mxdnr, mxdorb, mxdl)

      eion_ae(n) = etotal - eae_ref
      if(ifail /= 0) fail_ion(n) = 'ae'

      if(fail_ion(n) == ' ') then

        ctype = 'pt'

        call atom_atm_scf(etotal, ifail,                                 &
            nameat, ctype, ititle, kerker, icorr, ispp,                  &
            znuc, rmax, aa, bb, ncore, nval, ni, li, zd, zu, evd,        &
            iowrite, 0, fileae, iopsd, filepsd,                          &
            mxdnr, mxdorb, mxdl)

        eion_psd(n) = etotal - epsd_ref
        if(ifail /= 0) fail_ion(n) = 'pseudo'

      endif

      if(fail_ion(n) == ' ') then

        call atom_kb_test_scf(etotal, ifail,                             &
            nameat, ititle, icorr, ispp,                                 &
            nval, ni, li, zd, zu, evd,                                   &
            iowrite, iopsdkb, filepsdkb,                                 &
            mxdnr, mxdorb, mxdl)

        eion_kb(n) = etotal - ekb_ref
        if(ifail /= 0) fail_ion(n) = 'KB'

      endif


      zd(n) = zd(n) + zdi
      zu(n) = zu(n) + zui

    endif

  enddo

! minimum of ionization to make it worth

  if(nval > 1) then

    do n = 1,nval
      if(zu(n) + zd(n) < ONE/10) then
        lempty(n) = .TRUE.
      else
        lempty(n) = .FALSE.
      endif
      if(zu(n) + zd(n) > 2*(2*li(n)+1)*ONE - ONE/10) then
        lfull(n) = .TRUE.
      else
        lfull(n) = .FALSE.
      endif
    enddo


!   double loop over valence states for excitation energies

    nx = 0
    do n = 1,nval
    do n2 = 1,nval

      if(n2 /= n .and. (.not. lempty(n)) .and. (.not. lfull(n2))) then

        nx = nx + 1

        zz = 2*(2*li(n2)+1)*ONE - zu(n2) - zd(n2)
        zz = min(zz,zu(n)+zd(n),ONE)

        zui = min(zu(n),zz)
        zdi = zz - zui

        zd(n) = zd(n) - zdi
        zu(n) = zu(n) - zui

        if(li(n2) > 2) then
          zdx = zdi / 2
          zux = zui / 2
        else
          zdx = zdi
          zux = zui
        endif

        zd(n2) = zd(n2) + zdx
        zu(n2) = zu(n2) + zux

!       a calculation without self-consistent solution skips the configuration

        ctype = 'ae'

        call atom_atm_scf(etotal, ifail,                             &
            nameat, ctype, ititle, kerker, icorr, ispp,              &
            znuc, rmax, aa, bb, ncore, nval, ni, li, zd, zu, evd,    &
            iowrite, 0, fileae, iopsd, filepsd,                      &
            mxdnr, mxdorb, mxdl)

        exct_ae((n-1)*nval+n2) = etotal - eae_ref
        if(ifail /= 0) fail_exct((n-1)*nval+n2) = 'ae'

        if(fail_exct((n-1)*nval+n2) == ' ') then

          ctype = 'pt'

          call atom_atm_scf(etotal, ifail,                           &
              nameat, ctype, ititle, kerker, icorr, ispp,            &
              znuc, rmax, aa, bb, ncore, nval, ni, li, zd, zu, evd,  &
              iowrite, 0, fileae, iopsd, filepsd,                    &
              mxdnr, mxdorb, mxdl)

          exct_psd((n-1)*nval+n2) = etotal - epsd_ref
          if(ifail /= 0) fail_exct((n-1)*nval+n2) = 'pseudo'

        endif

        if(fail_exct((n-1)*nval+n2) == ' ') then

          call atom_kb_test_scf(etotal, ifail,                       &
              nameat, ititle, icorr, ispp,                           &
              nval, ni, li, zd, zu, evd,                             &
              iowrite, iopsdkb, filepsdkb,                           &
              mxdnr, mxdorb, mxdl)

          exct_kb((n-1)*nval+n2) = etotal - ekb_ref
          if(ifail /= 0) fail_exct((n-1)*nval+n2) = 'KB'

        endif

        zexct((n-1)*nval+n2) = -zdi - zui
        zexct2((n-1)*nval+n2) = zdx + zux

        zd(n) = zd(n) + zdi
        zu(n) = zu(n) + zui

        zd(n2) = zd(n2) - zdx
        zu(n2) = zu(n2) - zux

      endif

    enddo
    enddo

  endif

  write(6,*)
  write(6,*) '   Comparison of ionization energies (eV)'
  write(6,*)
  write(6,*) '   nl      all-electron   pseudo     KB-proj'
  write(6,*)

  if(iowrite /= 6) then
    write(iowrite,*)
    write(iowrite,*) '   Comparison of ionization energies (eV)'
    write(iowrite,*)
    write(iowrite,*) '   nl      all-electron   pseudo     KB-proj'
    write(iowrite,*)
  endif

  do n = 1,nval

    zui = min(zu(n),ONE)
    zdi = min(zd(n),ONE-zui)

!   minimum of ionization to make it worth

    if(zui + zdi > ONE/10 .and. fail_ion(n) == ' ') then

      write(6,'(3x,i2,a1,3(5x,f8.3))') ni(n), IL(li(n)),                 &
              eion_ae(n)*RYDBERG, eion_psd(n)*RYDBERG,                   &
              eion_kb(n)*RYDBERG

      if(iowrite /= 6) then

        write(iowrite,'(3x,i2,a1,3(5x,f8.3))') ni(n),                    &
              IL(li(n)), eion_ae(n)*RYDBERG, eion_psd(n)*RYDBERG,        &
              eion_kb(n)*RYDBERG

      endif

    elseif(zui + zdi > ONE/10) then

      write(6,'(3x,i2,a1,5x,"skipped, no self-consistent ",a,           &
         &      " solution")') ni(n), IL(li(n)), trim(fail_ion(n))

      if(iowrite /= 6) then

        write(iowrite,'(3x,i2,a1,5x,"skipped, no self-consistent ",a,    &
           &      " solution")') ni(n), IL(li(n)), trim(fail_ion(n))

      endif

    endif

  enddo

  write(6,*)

  if(nx > 0 .and. nval > 1) then

    write(6,*)
    write(6,*) '   Comparison of excitation energies (eV)'
    write(6,*)
    write(6,*) '   nl      nl     all-electron   pseudo      KB-proj      charge'
    write(6,*)

    if(iowrite /= 6) then
      write(iowrite,*)
      write(iowrite,*) '   Comparison of excitation energies (eV)'
      write(iowrite,*)
      write(iowrite,*) '   nl      nl     all-electron   pseudo      KB-proj      charge'
      write(iowrite,*)
    endif


    do n = 1,nval
    do n2 = 1,nval

      if(n2 /= n .and. (.not. lempty(n)) .and. (.not. lfull(n2)) .and.  &
         fail_exct((n-1)*nval+n2) /= ' ') then

        write(6,'(3x,i2,a1,5x,i2,a1,5x,"skipped, no self-consistent ",a, &
           &      " solution")') ni(n), IL(li(n)), ni(n2), IL(li(n2)),    &
           trim(fail_exct((n-1)*nval+n2))

        if(iowrite /= 6) then

          write(iowrite,'(3x,i2,a1,5x,i2,a1,5x,                          &
             &     "skipped, no self-consistent ",a," solution")')        &
             ni(n), IL(li(n)), ni(n2), IL(li(n2)),                        &
             trim(fail_exct((n-1)*nval+n2))

        endif

      elseif(n2 /= n .and. (.not. lempty(n)) .and. (.not. lfull(n2))) then

        write(6,'(3x,i2,a1,5x,i2,a1,3(5x,f8.3),5x,2f6.2)')               &
            ni(n),                                                       &
            IL(li(n)), ni(n2), IL(li(n2)),                               &
            exct_ae((n-1)*nval+n2)*RYDBERG,                              &
            exct_psd((n-1)*nval+n2)*RYDBERG,                             &
            exct_kb((n-1)*nval+n2)*RYDBERG,                              &
            zexct((n-1)*nval+n2),                                        &
            zexct2((n-1)*nval+n2)

        if(iowrite /= 6) then

          write(iowrite,'(3x,i2,a1,5x,i2,a1,3(5x,f8.3),5x,2f6.2)')       &
            ni(n),                                                       &
            IL(li(n)), ni(n2), IL(li(n2)),                               &
            exct_ae((n-1)*nval+n2)*RYDBERG,                              &
            exct_psd((n-1)*nval+n2)*RYDBERG,                             &
            exct_kb((n-1)*nval+n2)*RYDBERG,                              &
            zexct((n-1)*nval+n2),                                        &
            zexct2((n-1)*nval+n2)

        endif

      endif

    enddo
    enddo

  endif

  write(6,*)

  deallocate(ni)
  deallocate(li)
  deallocate(zu,zd)
  deallocate(evd)

  deallocate(eion_ae)
  deallocate(eion_psd)
  deallocate(eion_kb)

  deallocate(exct_ae)
  deallocate(exct_psd)
  deallocate(exct_kb)

  deallocate(lempty)
  deallocate(lfull)
  deallocate(zexct)
  deallocate(zexct2)

  deallocate(fail_ion)
  deallocate(fail_exct)

  return

end subroutine atom_atm_test
