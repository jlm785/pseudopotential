!>  Identifies if the parameters used to generate a pseudopotential
!>  are those of one of the tested pseudopotentials of atom_p_tbl_tested.
!>
!>
!>  \author       claude
!>  \version      6.2.0
!>  \date         12 November 2026.
!>  \copyright    GNU Public License v2

subroutine atom_psd_tested_match(iowrite, nameat,                       &
      norb, ncore, no, lo, zo,                                          &
      rc, ifcore, cfac, rcfac, label, nvar, vlist, status, lrctab,      &
      mxdorb, mxdlc, mxdvar)

! Written 12 November 2026. JLM+claude

! Compared: valence configuration (n, l and occupation summed over
! spin or j), core radii of the valence l, partial core correction
! and, if it is used, cfac and rcfac.
! Not compared: exchange and correlation, relativistic or spin
! polarized calculation, number of core orbitals (fixed by the
! valence configuration).
!
! It also identifies if the core radii of the valence l are those
! of the old table (atom_p_tbl_psd_tm2), and if the pseudopotential
! is not a tested one returns the confidence on those radii.


  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)

  integer, parameter          ::  LC = 4                                 !  maximum angular momentum in the tables

! input

  integer, intent(in)               ::  mxdorb                           !<  dimension of the number of orbitals
  integer, intent(in)               ::  mxdlc                            !<  dimension of the angular momentum of core radii
  integer, intent(in)               ::  mxdvar                           !<  dimension of the list of labels

  integer, intent(in)               ::  iowrite                          !<  default output tape

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbol of the atom

  integer, intent(in)               ::  norb                             !<  number of orbitals
  integer, intent(in)               ::  ncore                            !<  number of core orbitals
  integer, intent(in)               ::  no(mxdorb)                       !<  principal quantum number n
  integer, intent(in)               ::  lo(mxdorb)                       !<  angular quantum number l
  real(REAL64), intent(in)          ::  zo(mxdorb)                       !<  orbital occupation

  real(REAL64), intent(in)          ::  rc(0:mxdlc)                      !<  core radius r_c(l) (before rounding to the grid)
  integer, intent(in)               ::  ifcore                           !<  0: no partial core correction; 1 : xc partial core correction; 2 : not recommended
  real(REAL64), intent(in)          ::  cfac                             !<  criteria for pseudo-core charge
  real(REAL64), intent(in)          ::  rcfac                            !<  pseudo-core radius

  real(REAL64)                      ::  cfactau                          !  SHOULD BE INTENT IN

! output

  character(len=4), intent(out)     ::  label                            !<  label of the tested pseudopotential, 'NONE' if not tested
  integer, intent(out)              ::  nvar                             !<  number of tested pseudopotentials of the element
  character(len=4), intent(out)     ::  vlist(mxdvar)                    !<  labels of the tested pseudopotentials
  character(len=30), intent(out)    ::  status                           !<  tests that were performed, if label = 'NONE' confidence on core radii of old table
  logical, intent(out)              ::  lrctab                           !<  the core radii are those of the old table

! local allocatable arrays

  integer, allocatable              ::  nv(:), lv(:)                     !  valence n, l
  real(REAL64), allocatable         ::  zv(:)                            !  valence occupation summed over spin or j

! local variables

  integer                           ::  nval                             !  number of different valence n, l

  character(len=4), allocatable     ::  vtmp(:)
  integer                           ::  nvtmp
  logical                           ::  lfound, lmatch

  integer                           ::  ncore_t, nval_t                  !  values in table
  integer                           ::  no_t(LC+1), lo_t(LC+1)
  real(REAL64)                      ::  zo_t(LC+1)
  real(REAL64)                      ::  rc_t(0:LC)
  character(len=2)                  ::  pcc_t
  real(REAL64)                      ::  cfac_t, rcfac_t
  integer                           ::  llocal_t
  character(len=30)                 ::  status_t

  integer                           ::  ifcore_t

  real(REAL64)                      ::  rc_tab(0:LC)                     !  core radii in the old table
  character(len=30)                 ::  status_tab                       !  confidence on the core radii of the old table
  integer                           ::  jm, l

  character(len=30)                 ::  stsq                             !  squeezed status

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  EPS = 0.001_REAL64                      !  tolerance for occupations and cfac
  real(REAL64), parameter    ::  EPSRC = 0.005_REAL64                    !  tolerance for core radii

! counters

  integer     ::  i, j, k


  label = 'NONE'
  status = '                              '

  call atom_p_tbl_tested(nameat, 'LIST', nvar, vlist, lfound,            &
      ncore_t, nval_t, no_t, lo_t, zo_t, rc_t, pcc_t, llocal_t,          &
      cfac, rcfac, cfactau, status_t,                                    &
      mxdvar)

! valence configuration summed over spin or j

  allocate(nv(mxdorb), lv(mxdorb), zv(mxdorb))

  nval = 0
  do i = ncore+1,norb
    jm = 0
    do j = 1,nval
      if(nv(j) == no(i) .and. lv(j) == lo(i)) jm = j
    enddo
    if(jm == 0) then
      nval = nval + 1
      nv(nval) = no(i)
      lv(nval) = lo(i)
      zv(nval) = zo(i)
    else
      zv(jm) = zv(jm) + zo(i)
    endif
  enddo

! are the core radii those of the old table?

  call atom_p_tbl_psd_tm2(nameat, rc_tab, status_tab)

  lrctab = .TRUE.
  do j = 1,nval
    l = lv(j)
    if(l > mxdlc .or. l > LC) then
      lrctab = .FALSE.
    else
      if(rc_tab(l) <= ZERO) lrctab = .FALSE.
      if(abs(rc(l) - rc_tab(l)) > EPSRC) lrctab = .FALSE.
    endif
  enddo

  if(lrctab) status = status_tab

  allocate(vtmp(mxdvar))

  do k = 1,nvar

    call atom_p_tbl_tested(nameat, vlist(k), nvtmp, vtmp, lfound,        &
        ncore_t, nval_t, no_t, lo_t, zo_t, rc_t, pcc_t, llocal_t,        &
        cfac_t, rcfac_t, cfactau, status_t,                              &
        mxdvar)

    lmatch = lfound .and. nval_t == nval

    if(lmatch) then
      do j = 1,nval_t
        jm = 0
        do i = 1,nval
          if(nv(i) == no_t(j) .and. lv(i) == lo_t(j)) jm = i
        enddo
        if(jm == 0) then
          lmatch = .FALSE.
        else
          if(abs(zv(jm) - zo_t(j)) > EPS) lmatch = .FALSE.
          l = lo_t(j)
          if(l > mxdlc) then
            lmatch = .FALSE.
          else
            if(abs(rc(l) - rc_t(l)) > EPSRC) lmatch = .FALSE.
          endif
        endif
      enddo
    endif

    if(lmatch) then
      ifcore_t = 0
      if(pcc_t == 'pe') ifcore_t = 1
      if(pcc_t == 'ph') ifcore_t = 2
      if(ifcore /= ifcore_t) lmatch = .FALSE.
      if(ifcore_t /= 0) then
        if(abs(cfac - cfac_t) > EPS) lmatch = .FALSE.
        if(rcfac > ZERO .or. rcfac_t > ZERO) then
          if(abs(rcfac - rcfac_t) > EPS) lmatch = .FALSE.
        endif
      endif
    endif

    if(lmatch) then
      label = vlist(k)
      status = status_t

      exit

    endif

  enddo

  deallocate(nv, lv, zv)
  deallocate(vtmp)

! prints the eventual match

  write(iowrite,*)

  if(label /= 'NONE') then

    write(iowrite,'("  The parameters are those of the tested ",         &
         &   "pseudopotential ",a,1x,a4)') trim(adjustl(nameat)), label
    write(iowrite,'("  tests: ",a)') trim(status)

  elseif(nvar > 0) then

    write(iowrite,'("  The parameters differ from the tested ",          &
         &   a," pseudopotentials:",9(2x,a4))') trim(adjustl(nameat)),   &
         (vlist(k),k=1,nvar)
    write(iowrite,'("  This pseudopotential was not tested")')

  else

    write(iowrite,'("  There are no tested pseudopotentials for ",a,     &
         &   ", this one was not tested")') trim(adjustl(nameat))

  endif

  if(label == 'NONE') then

    if(lrctab) then

      stsq = status
      call chr_squeeze(stsq)
      if(len_trim(stsq) == 0) stsq = 'no information'

      write(iowrite,'("  core radii: ",a,", test local and other choices")') trim(stsq)

    else

      write(iowrite,'("  core radii differ from the old table, test ",   &
           &   "all the choices")')

    endif

  endif

  write(iowrite,*)

  return

end subroutine atom_psd_tested_match
