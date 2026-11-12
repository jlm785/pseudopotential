!>  Default local component of the Kleinman-Bylander pseudopotential.
!>  If the label written by atom_psd_sub in the title of the real space
!>  pseudopotential file (psd.pot) identifies a tested pseudopotential
!>  its local component is used, otherwise the one from the old table.
!>
!>  \author       claude
!>  \version      6.2.0
!>  \date         12 November 2026.
!>  \copyright    GNU Public License v2

subroutine atom_kb_llocal_default(nameat, llocal,                       &
      iowrite, ioreal, filereal)

! Written 12 November 2026. JLM+claude

  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)

  integer, parameter          ::  LC = 4                                 !  maximum angular momentum in the tables
  integer, parameter          ::  MXDVAR = 9                             !  maximum number of tested pseudopotentials of an element

! input

  integer, intent(in)               ::  iowrite                          !<  default output tape

  integer, intent(in)               ::  ioreal                           !<  tape for pseudopotential in real space (psd.pot) format
  character(len=*), intent(in)      ::  filereal                         !<  name of the file with pseudopotential in real space

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbol of the atom

! output

  integer, intent(out)              ::  llocal                           !<  angular momentum for local potential (negative: maximum of l-dependent)

! local variables

  character(len=10)                 ::  psdtitle(20)                     !  second title of psd.pot
  character(len=4)                  ::  label                            !  label of the tested pseudopotential
  integer                           ::  ios
  integer                           ::  jhard                            !  flag for accuracy/speed compromise

  integer                           ::  nvar                             !  number of tested pseudopotentials
  character(len=4)                  ::  vlist(MXDVAR)                    !  labels of tested pseudopotentials
  logical                           ::  lfound

  integer                           ::  ncore, nval
  integer                           ::  no(LC+1), lo(LC+1)
  real(REAL64)                      ::  zo(LC+1)
  real(REAL64)                      ::  rc(0:LC)
  character(len=2)                  ::  pcc
  real(REAL64)                      ::  cfac, rcfac
  real(REAL64)                      ::  cfactau
  character(len=30)                 ::  status

! counters

  integer     ::  i


  do i = 1,20
    psdtitle(i) = '          '
  enddo

  open(unit = ioreal, file = trim(filereal), status='OLD', form='FORMATTED')

  read(ioreal,*)
  read(ioreal,*)
  read(ioreal,'(1x,20a10)', iostat = ios) psdtitle

  close(unit = ioreal)

  label = 'NONE'
  if(ios == 0 .and. psdtitle(11)(2:5) == 'var=') label = psdtitle(11)(6:9)

  lfound = .FALSE.
  if(label /= 'NONE') then
    call atom_p_tbl_tested(nameat, label, nvar, vlist, lfound,           &
        ncore, nval, no, lo, zo, rc, pcc, llocal,                        &
        cfac, rcfac, cfactau, status,                                    &
        MXDVAR)
  endif

  write(iowrite,*)
  if(lfound) then
    write(iowrite,'("  Local potential l = ",i2," of the tested ",       &
         &   "pseudopotential ",a,1x,a4)') llocal,                       &
             trim(adjustl(nameat)), label
  else
    if(label /= 'NONE') then
      write(iowrite,'("  ",a4," is not a tested pseudopotential for ",   &
           &   a," (table changed?)")') label, trim(adjustl(nameat))
    endif
    jhard = 0
    call atom_p_tbl_kb_local(nameat, llocal, jhard)
    write(iowrite,'("  Local potential l = ",i2," from the old table, ", &
         &   "not tested (negative: smoothed maximum)")') llocal
  endif

  return

end subroutine atom_kb_llocal_default
