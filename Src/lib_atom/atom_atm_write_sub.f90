!>  Writes a default input configuration file from the chemical symbol
!>  and the exchange and correlation functional, for the atomic ground
!>  state configuration (non-spin-polarized up to Ca, relativistic above).
!>  User may edit that file to have the other cases.
!>  The configuration, core radii and partial core correction come from
!>  a tested pseudopotential (atom_p_tbl_tested) if there is one,
!>  otherwise from the old tables (educated guess).
!>
!>  \author       Jose Luis Martins
!>  \version      6.2.0
!>  \date         22 June 2021, 27 March 2026, 14 November 2026.
!>  \copyright    GNU Public License v2

subroutine atom_atm_write_sub(ioread, filein, nameat, vtype, icorr)

! nval2, jhard. 21 December 2021. JLM
! default mesh. 21 April 2022. JLM
! correct format of core radii and cfactau. 27 March 2026. JLM
! MXDVAL, core radii by angular momentum, old convention (cfac in
! columns 41-50) unless there is a g orbital. 27 October 2026. JLM+claude
! rcval up to 9 fields. 27 October 2026. claude
! tested pseudopotentials, vtype. 12 November 2026. JLM+claude
! cfactau from the tested pseudopotential. 12 November 2026. claude
! confidence on core radii only printed by atom_psd_tested_match. 12 November 2026. claude
! icorr is an argument. 14 November 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter          :: LC = 4                                  !  maximum angular momentum in the tables
  integer, parameter          :: MXDVAL = LC + 1                         !  maximum number of valence orbitals in the tables
  integer, parameter          :: MXDVAR = 9                              !  maximum number of tested pseudopotentials of an element

! input

  integer, intent(in)               ::  ioread                           !<  default tape for reading
  character(len=*), intent(in)      ::  filein                           !<  name of default tape for reading

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbol of the atom
  character(len=4), intent(in)      ::  vtype                            !<  label of tested pseudopotential, '    ' default, 'NONE' old tables
  character(len=2), intent(in)      ::  icorr                            !<  exchange and correlation functional (ca, pw, pb)

! local variables

  character(len=10)                 ::  ititle(5)                        !  title of calculation
  character(len=9)                  ::  bdate                            !  date
  character(len=1)                  ::  ispp                             !  spin polarization  ' ', 's', 'r'

  integer                           ::  iz                               !  atomic number

  real(REAL64)                      ::  rmax, aa, bb                     !  logarithmic mesh parameters

  character(len=2)                  ::  pcc

  integer                           ::  ncore                            !  canonical number of core orbitals
  integer                           ::  nval                             !  canonical number of interesting valence orbitals
  integer                           ::  jhard                            !  flag for accuracy/speed compromise
  integer                           ::  no(MXDVAL), lo(MXDVAL)           !  configuration
  real(REAL64)                      ::  zo(MXDVAL)                       !  occupation
  real(REAL64)                      ::  rc(0:LC)                         !  core radii from the table
  real(REAL64)                      ::  rcval(0:8)                       !  core radii of the valence angular momenta (up to 9 fields)
  integer                           ::  lmax                             !  maximum valence angular momentum
  integer                           ::  nrc                              !  number of core radii fields (4 or 5)
  character(len=30)                 ::  status                           !  quality of pseudopotential in table

  character(len=128)                ::  cline                            !  better way to code the last line

  integer                           ::  nvar                             !  number of tested pseudopotentials
  character(len=4)                  ::  vlist(MXDVAR)                    !  labels of tested pseudopotentials
  character(len=4)                  ::  vuse                             !  label of the pseudopotential that is used
  logical                           ::  lfound                           !  vtype is one of the labels
  real(REAL64)                      ::  cfac, rcfac                      !  pseudo-core parameters
  real(REAL64)                      ::  cfactau                          !  partial core kinetic energy parameter
  integer                           ::  llocal                           !  local component (not used here)
  character(len=30)                 ::  status_tst                       !  tests performed on the tested pseudopotential

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64

! counters

  integer     ::  n, l


! tested pseudopotentials

  call atom_p_tbl_tested(nameat, 'LIST', nvar, vlist, lfound,            &
      ncore, nval, no, lo, zo, rc, pcc, llocal,                          &
      cfac, rcfac, cfactau, status_tst,                                  &
      MXDVAR)

  if(vtype == '    ') then
    if(nvar > 0) then
      vuse = vlist(1)
    else
      vuse = 'NONE'
    endif
  elseif(vtype == 'NONE') then
    vuse = 'NONE'
  else
    vuse = 'NONE'
    do n = 1,nvar
      if(vtype == vlist(n)) vuse = vtype
    enddo
    if(vuse == 'NONE') then
      if(nvar == 0) then
        write(6,*)
        write(6,'("  There are no tested pseudopotentials for ",a,       &
             &   ", ",a4," ignored")') trim(adjustl(nameat)), vtype
      else
        write(6,*)
        write(6,'("  Stopped in atom_atm_write_sub: ",a4," is not a ",   &
             &   "tested pseudopotential for ",a)') vtype,               &
             trim(adjustl(nameat))
        write(6,'("  Available:",9(2x,a4))') (vlist(n),n=1,nvar), 'NONE'
        write(6,*)

        STOP

      endif
    endif
  endif

  if(vuse /= 'NONE') then

    call atom_p_tbl_tested(nameat, vuse, nvar, vlist, lfound,            &
        ncore, nval, no, lo, zo, rc, pcc, llocal,                        &
        cfac, rcfac, cfactau, status_tst,                                &
        MXDVAR)

    write(6,*)
    write(6,'("  Using the tested pseudopotential ",a,1x,a4,             &
         &   " (",a,")")') trim(adjustl(nameat)), vuse, trim(status_tst)
    if(nvar > 1) write(6,'("  Tested pseudopotentials for ",a,":",       &
         &   9(2x,a4))') trim(adjustl(nameat)), (vlist(n),n=1,nvar)

  else

    jhard = 0
    call atom_p_tbl_psd_pcc(nameat, pcc, jhard)

    jhard = 0
    call atom_p_tbl_config(nameat, ncore, nval, no, lo, zo, jhard)

    call atom_p_tbl_psd_tm2(nameat, rc, status)

    cfac = ONE
    rcfac = ZERO
    cfactau = ONE

    write(6,*)
    write(6,'("  Pseudopotential parameters for ",a," from the old ",    &
         &   "tables (educated guess)")') trim(adjustl(nameat))
    if(nvar > 0) write(6,'("  Tested pseudopotentials for ",a,":",       &
         &   9(2x,a4))') trim(adjustl(nameat)), (vlist(n),n=1,nvar)

  endif

  call zedate(bdate)

  ititle(1) = '     '//nameat//'   '
  ititle(2) = '   file ge'
  ititle(3) = 'nerated by'
  ititle(4) = ' atm_write'
  ititle(5) = ' '//bdate

  call atom_p_tbl_charge(nameat, iz)

  ispp = ' '
  if(iz > 20) ispp = 'r'

  call atom_atm_default_mesh(rmax, aa, bb)

  open(unit=ioread, file=trim(filein), status='NEW', form='FORMATTED')

  write(ioread,'(3x,a2,5a10)') pcc,ititle
  write(ioread,'(8x,a3)') 'tm2'
  write(ioread,'(" n=",a2," c=",a2,a1)') nameat,icorr,ispp

  write(ioread,'(6f10.1)') ONE*iz, ZERO, ZERO, rmax, aa, bb

  write(ioread,'(1x,i4,1x,i4)') ncore,nval

  do n = 1,nval
    write(ioread,'(1x,i4,1x,i4,1x,f8.2,2x,f8.2)') no(n), lo(n), zo(n), ZERO
  enddo

  do n = 1,128
    cline(n:n) = ' '
  enddo

! core radii by angular momentum (zero if l is not in the valence).
! Old convention, 4 fields and cfac in columns 41-50, unless there
! is a g (or higher) orbital.

  do l = 0,8
    rcval(l) = ZERO
  enddo
  lmax = 0
  do n = 1,nval
    if(lo(n) <= LC) rcval(lo(n)) = rc(lo(n))
    if(lo(n) > lmax) lmax = lo(n)
  enddo

  nrc = 4
  if(lmax > 3) nrc = lmax+1

  if(nrc > 9) then
    write(6,*)
    write(6,*) '   Stopped in atom_atm_write_sub, lmax = ',lmax

    stop

  endif

  write(cline,'(9f10.5)') (rcval(l),l=0,lmax)
  write(cline(10*nrc+1:10*nrc+30),'(3f10.5)') cfac, rcfac, cfactau
  write(ioread,'(a128)') cline

! for back-compatibility you need to finish with an empty line

  write(ioread,'(50x)')

  close(unit=ioread)

  return

end subroutine atom_atm_write_sub
