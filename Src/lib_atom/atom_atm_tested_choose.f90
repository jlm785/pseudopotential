!>  Lets the user choose (interactive) among the tested pseudopotentials
!>  of an element or the old tables (educated guess).
!>  The default (empty line or end of input) is the first tested pseudopotential.
!>
!>  \author       claude
!>  \version      6.2.0
!>  \date         12 November 2026.
!>  \copyright    GNU Public License v2

subroutine atom_atm_tested_choose(nameat, vtype)

! Written 12 November 2026. JLM+claude

  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)

  integer, parameter          ::  LC = 4                                 !  maximum angular momentum in the tables
  integer, parameter          ::  MXDVAR = 9                             !  maximum number of tested pseudopotentials of an element

! input

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbol of the atom

! output

  character(len=4), intent(out)     ::  vtype                            !<  label of the chosen pseudopotential, 'NONE' for old tables

! local variables

  integer                           ::  nvar                             !  number of tested pseudopotentials
  character(len=4)                  ::  vlist(MXDVAR)                    !  labels of tested pseudopotentials
  character(len=4)                  ::  vtmp(MXDVAR)
  integer                           ::  nvtmp
  logical                           ::  lfound

  integer                           ::  ncore, nval
  integer                           ::  no(LC+1), lo(LC+1)
  real(REAL64)                      ::  zo(LC+1)
  real(REAL64)                      ::  rc(0:LC)
  character(len=2)                  ::  pcc
  real(REAL64)                      ::  cfac, rcfac
  real(REAL64)                      ::  cfactau
  integer                           ::  llocal
  character(len=30)                 ::  status

  character(len=80)                 ::  line
  integer                           ::  ichoice, ios

! counters

  integer     ::  k


  call atom_p_tbl_tested(nameat, 'LIST', nvar, vlist, lfound,            &
      ncore, nval, no, lo, zo, rc, pcc, llocal,                          &
      cfac, rcfac, cfactau, status,                                      &
      MXDVAR)

  if(nvar == 0) then
    vtype = 'NONE'

    return

  endif

  write(6,*)
  write(6,'("  Tested pseudopotentials for ",a,":")') trim(adjustl(nameat))
  write(6,*)
  do k = 1,nvar
    call atom_p_tbl_tested(nameat, vlist(k), nvtmp, vtmp, lfound,        &
        ncore, nval, no, lo, zo, rc, pcc, llocal,                        &
        cfac, rcfac, cfactau, status,                                    &
        MXDVAR)
    write(6,'(4x,i2,3x,a4,3x,a30)') k, vlist(k), status
  enddo
  call atom_p_tbl_psd_tm2(nameat, rc, status)
  call chr_squeeze(status)
  if(len_trim(status) == 0) status = 'no information'
  write(6,'(4x,i2,3x,"NONE",3x,"old tables, core radii: ",a)')           &
       nvar+1, trim(status)
  write(6,*)
  write(6,'("  Enter your choice (default 1)")')
  write(6,*)

  ichoice = 1
  read(5,'(a)',iostat=ios) line
  if(ios == 0) then
    read(line,*,iostat=ios) ichoice
    if(ios /= 0) ichoice = 1
  endif

  if(ichoice < 1 .or. ichoice > nvar+1) then
    write(6,'("  Choice ",i5," out of range, using 1")') ichoice
    ichoice = 1
  endif

  if(ichoice == nvar+1) then
    vtype = 'NONE'
  else
    vtype = vlist(ichoice)
  endif

  return

end subroutine atom_atm_tested_choose
