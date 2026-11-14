!>  Pseudopotentials that were tested against all-electron calculations:
!>  configuration, core radii, partial core correction and
!>  local component of the Kleinman-Bylander form.
!>
!>  An element may have several tested pseudopotentials identified by
!>  a label ('HARD', 'SOFT',...), the first one is the default.
!>  If an element has none (nvar = 0) the old tables (atom_p_tbl_config,
!>  atom_p_tbl_psd_tm2, atom_p_tbl_psd_pcc, atom_p_tbl_kb_local)
!>  provide an educated guess.
!>
!>  The labels 'NONE' (old tables) and '    ' (default) are reserved,
!>  'YES', 'NO', 'LDA', 'GGA', 'PBE', anything with a 'Z' in the second
!>  position (basis set) and anything starting with 'XC' (functional)
!>  cannot be used as labels.
!>
!>  \author       Jose Luis Martins+claude
!>  \version      6.2.0
!>  \date         12 November 2026. 14 November 2026.
!>  \copyright    GNU Public License v2

subroutine atom_p_tbl_tested(name, vtype, nvar, vlist, lfound,          &
      ncore, nval, no, lo, zo, rc, pcc, llocal,                         &
      cfac, rcfac, cfactau, status,                                     &
      mxdvar)

! Written 12 November 2026. JLM+claude
! Si, Ge. 14 November 2026. JLM+claude

  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)

  integer, parameter          ::  LC = 4                                 !  maximum angular momentum in the tables
  integer, parameter          ::  MXDTAB = 9                             !  maximum number of tested pseudopotentials of an element

! input

  integer, intent(in)               ::  mxdvar                           !<  dimension of the list of labels

  character(len=2), intent(in)      ::  name                             !<  chemical symbol of the atom
  character(len=4), intent(in)      ::  vtype                            !<  label of the pseudopotential, 'LIST' only lists the labels

! output

  integer, intent(out)              ::  nvar                             !<  number of tested pseudopotentials (0: none)
  character(len=4), intent(out)     ::  vlist(mxdvar)                    !<  labels of the tested pseudopotentials, vlist(1) is the default
  logical, intent(out)              ::  lfound                           !<  vtype is one of the labels and the following output is defined

  integer, intent(out)              ::  ncore                            !<  number of core orbitals
  integer, intent(out)              ::  nval                             !<  number of valence orbitals
  integer, intent(out)              ::  no(LC+1), lo(LC+1)               !<  valence configuration
  real(REAL64), intent(out)         ::  zo(LC+1)                         !<  occupation
  real(REAL64), intent(out)         ::  rc(0:LC)                         !<  core radii for angular momentum l
  character(len=2), intent(out)     ::  pcc                              !<  no pcc: pg; pcc: pe.
  real(REAL64), intent(out)         ::  cfac                             !<  criteria for pseudo-core charge (if pcc = pe)
  real(REAL64), intent(out)         ::  rcfac                            !<  pseudo-core radius (if pcc = pe, <= 0 from cfac)
  real(REAL64), intent(out)         ::  cfactau                          !<  criteria for partial core kinetic energy radius
  integer, intent(out)              ::  llocal                           !<  angular momentum of the local potential (negative: maximum of l-dependent)
  character(len=30), intent(out)    ::  status                           !<  tests that were performed

! local variables

  character(len=4)                  ::  vtab(MXDTAB)                     !  labels in the table

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64

! counters

  integer     ::  i


  nvar = 0
  do i = 1,MXDTAB
    vtab(i) = '    '
  enddo

  lfound = .FALSE.
  ncore = 0
  nval = 0
  do i = 1,LC+1
    no(i) = 0
    lo(i) = 0
    zo(i) = ZERO
  enddo
  do i = 0,LC
    rc(i) = ZERO
  enddo
  pcc = 'pg'
  cfac = ONE
  rcfac = ZERO
  cfactau = ONE
  llocal = -1
  status = '                              '


  if (name == 'C ' .or. name == ' C') then


    nvar = 2
    vtab(1) = 'HARD'
    vtab(2) = 'SOFT'

    if(vtype == 'HARD') then

      lfound = .TRUE.
      ncore = 1
      nval = 3
      no(1) = 2
      lo(1) = 0
      zo(1) = 2*ONE
      no(2) = 2
      lo(2) = 1
      zo(2) = 2*ONE
      no(3) = 3
      lo(3) = 2
      zo(3) = 0*ONE
      rc(0) = 1.10_REAL64
      rc(1) = 1.10_REAL64
      rc(2) = 1.10_REAL64
      pcc = 'pe'
      cfac = ONE
      rcfac = ZERO
      cfactau = ONE
      llocal = 2
      status = 'diamond LAPW, C2H2 all-electr.'

    elseif(vtype == 'SOFT') then

      lfound = .TRUE.
      ncore = 1
      nval = 2
      no(1) = 2
      lo(1) = 0
      zo(1) = 2*ONE
      no(2) = 2
      lo(2) = 1
      zo(2) = 2*ONE
      rc(0) = 1.50_REAL64
      rc(1) = 1.54_REAL64
      pcc = 'pg'
      cfac = ONE
      rcfac = ZERO
      cfactau = ONE
      llocal = 1
      status = 'diamond PRB 43 1993 (1991)    '

    endif


  elseif (name == 'Si') then


    nvar = 2
    vtab(1) = 'HARD'
    vtab(2) = 'SOFT'

    if(vtype == 'HARD') then

      lfound = .TRUE.
      ncore = 3
      nval = 3
      no(1) = 3
      lo(1) = 0
      zo(1) = 2*ONE
      no(2) = 3
      lo(2) = 1
      zo(2) = 2*ONE
      no(3) = 3
      lo(3) = 2
      zo(3) = 0*ONE
      rc(0) = 1.90_REAL64
      rc(1) = 1.90_REAL64
      rc(2) = 1.90_REAL64
      pcc = 'pe'
      cfac = ONE
      rcfac = ZERO
      cfactau = ONE/4
      llocal = -1
      status = 'LAPW comparison               '

    elseif(vtype == 'SOFT') then

      lfound = .TRUE.
      ncore = 3
      nval = 2
      no(1) = 3
      lo(1) = 0
      zo(1) = 2*ONE
      no(2) = 3
      lo(2) = 1
      zo(2) = 2*ONE
      rc(0) = 2.50_REAL64
      rc(1) = 2.50_REAL64
      pcc = 'pg'
      cfac = ONE
      rcfac = ZERO
      cfactau = ONE
      llocal = 1
      status = 'silicon  PRB 52 13281 (1995)  '

    endif


  elseif (name == 'Ge') then


    nvar = 2
    vtab(1) = 'HARD'
    vtab(2) = 'SOFT'

    if(vtype == 'HARD') then

      lfound = .TRUE.
      ncore = 6
      nval = 4
      no(1) = 4
      lo(1) = 0
      zo(1) = 2*ONE
      no(2) = 4
      lo(2) = 1
      zo(2) = 2*ONE
      no(3) = 4
      lo(3) = 2
      zo(3) = 0*ONE
      no(4) = 4
      lo(4) = 3
      zo(4) = 0*ONE
      rc(0) = 2.20_REAL64
      rc(1) = 2.20_REAL64
      rc(2) = 2.20_REAL64
      rc(3) = 2.20_REAL64
      pcc = 'pe'
      cfac = ONE
      rcfac = ZERO
      cfactau = 8*ONE
      llocal = -1
      status = 'LAPW comparison               '

    elseif(vtype == 'SOFT') then

      lfound = .TRUE.
      ncore = 6
      nval = 2
      no(1) = 4
      lo(1) = 0
      zo(1) = 2*ONE
      no(2) = 4
      lo(2) = 1
      zo(2) = 2*ONE
      rc(0) = 2.68_REAL64
      rc(1) = 2.68_REAL64
      pcc = 'pg'
      cfac = ONE
      rcfac = ZERO
      cfactau = 8*ONE
      llocal = 1
      status = 'germanium  PRB 52 13283 (1995)'

    endif

  endif

! paranoid check, nvar set here

  if(nvar > mxdvar) then
    write(6,*)
    write(6,'("  Stopped in atom_p_tbl_tested, nvar = ",i3,          &
          &   " > mxdvar = ",i3)') nvar, mxdvar
    write(6,*)

    STOP

  endif

  do i = 1,nvar
    vlist(i) = vtab(i)
  enddo

  return

end subroutine atom_p_tbl_tested
