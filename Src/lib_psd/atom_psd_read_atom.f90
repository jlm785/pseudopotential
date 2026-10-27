!>  Reads the pseudopotential core radii
!>
!>  The line after the valence orbitals has the core radii by angular
!>  momentum, rc(s), rc(p), rc(d), rc(f) in 4f10.5 (old convention),
!>  followed by cfac, rcfac and cfactau.  If there is a valence g
!>  orbital rc(g) is a fifth field and cfac etc. start in column 51.
!>  Files written by version 6.1.0 (columns 41-50 blank) are also read.
!>
!>  \author       Jose Luis Martins
!>  \version      6.1.0
!>  \date         27 August 2021, 25 March 2026, 27 October 2026.
!>  \copyright    GNU Public License v2

subroutine atom_psd_read_atom(ioread, filein, nameat, rc,                &
      cfac, rcfac, cfactau, lres)

! written 27 August 2021 from previous code. JLM
! cfactau, 25 March 2026. JLM
! old convention (cfac in columns 41-50) unless there is a g orbital,
! cfac, rcfac, cfactau optional. 27 October 2026. JLM+claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter                ::  mxdlc = 4                        !  maximum valence angular momentum

! input

  integer, intent(in)               ::  ioread                           !  default tape for reading
  character(len=*), intent(in)      ::  filein                           !  name of default tape for reading

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbot of the atom

! output

  real(REAL64), intent(out)         ::  rc(0:mxdlc)                      !<  core radius r_c(l)
  real(REAL64), intent(out)         ::  cfac                             !<  criteria for pseudo-core charge
  real(REAL64), intent(out)         ::  rcfac                            !<  pseudo-core radius
  real(REAL64), intent(out)         ::  cfactau                          !<  core correction factor for kinetic energy density

  logical, intent(out)              ::  lres                             !<  lecture seems successful

! local variables

  character(len=2)     ::  ctype                                         !  type of calculation flag
  integer              ::  nc, nv
  integer              ::  ni, li                                        !  valence orbital
  integer              ::  nrc                                           !  number of core radii fields (4 or 5)
  integer              ::  icfac                                         !  column where cfac starts
  real(REAL64)         ::  c1, c2, c3

  integer              ::  ioerror
  character(len=2)     ::  name                                          !  element on the atom.dat file

  character(len=128)   ::  cline

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  real(REAL64), parameter    ::  EPS = 1.0E-12_REAL64

! counters

  integer        ::  i, l


! default values

  do l = 0,mxdlc
    rc(l) = ZERO
  enddo
  cfac = ONE
  rcfac = ZERO
  cfactau = ONE

! read rc(s),rc(p),rc(d),rc(f),rc(g),cfac,rcfac

! cfac is used for the pseudocore - the pseudocore stops where
! the core charge density equals cfac times the renormalized
! valence charge density (renormalized to make the atom neutral).
! If cfac is input as negative, the full core charge is used,
! if cfac is input as zero, it is set equal to one.
! rcfac is used for the pseudocore cut off radius.  If set
! to less then or equal to zero cfac is used.  cfac must be
! set to greater then zero.

  open(unit=ioread, file=trim(filein), status='OLD', form='FORMATTED')

  read(ioread,'(3x,a2)') ctype
  if(ctype /= 'ae') then
    read(ioread,*)
  endif
  read(ioread,'(3x,a2)') name

  if(nameat /= name) then
    if(nameat(1:1) == ' ' .and. name(2:2) == ' ' .and. nameat(2:2) /= name(1:1)) then
      if(nameat(2:2) == ' ' .and. name(1:1) == ' ' .and. nameat(1:1) /= name(2:2)) then

        write(6,*) '  Stopped in atom_psd_read_atom:'
        write(6,*) '  Inconsistent elements in files: ',name,'  ',nameat

        STOP

      endif
   endif
  endif


  read(ioread,*)
  read(ioread,*) nc,nv

! a g orbital in the valence needs a fifth field

  nrc = 4
  do i = 1,nv
    read(ioread,'(2i5)') ni, li
    if(li > 3) nrc = li+1
  enddo

! uses cline to be able to keep back-compatibility

  do i = 1,128
    cline(i:i) = ' '
  enddo

  read(ioread,'(a128)')  cline

  close(unit=ioread)

! core radii by angular momentum

  read(cline,'(9f10.5)',iostat = ioerror) (rc(l),l=0,nrc-1)
  if(ioerror /= 0) then
    do l = 0,mxdlc
      rc(l) = ZERO
    enddo
  endif

! cfac, rcfac and cfactau, missing values keep the defaults.
! Version 6.1.0 wrote them from column 51 with columns 41-50 blank.

  icfac = 10*nrc+1

  read(cline(icfac:128),*,iostat = ioerror) c1, c2, c3
  if(ioerror == 0) then
    cfac = c1
    rcfac = c2
    cfactau = c3
  else
    read(cline(icfac:128),*,iostat = ioerror) c1, c2
    if(ioerror == 0) then
      cfac = c1
      rcfac = c2
    else
      read(cline(icfac:icfac+9),'(f10.5)',iostat = ioerror) c1
      if(ioerror == 0) cfac = c1
    endif
  endif

  lres = .TRUE.
  if(abs(rc(0)) < EPS) lres = .FALSE.

  if(cfac < EPS) cfac = ONE
  if(cfactau < EPS) cfactau = ONE

  return

end subroutine atom_psd_read_atom

