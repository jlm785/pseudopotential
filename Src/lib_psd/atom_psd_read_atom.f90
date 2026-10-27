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
! fixed 10 column fields as the 1982 code, negative cfac kept,
! rc(l) stored only up to MXDLC. 27 October 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter                ::  MXDLC = 4                        !  maximum valence angular momentum

! input

  integer, intent(in)               ::  ioread                           !  default tape for reading
  character(len=*), intent(in)      ::  filein                           !  name of default tape for reading

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbot of the atom

! output

  real(REAL64), intent(out)         ::  rc(0:MXDLC)                      !<  core radius r_c(l)
  real(REAL64), intent(out)         ::  cfac                             !<  criteria for pseudo-core charge
  real(REAL64), intent(out)         ::  rcfac                            !<  pseudo-core radius
  real(REAL64), intent(out)         ::  cfactau                          !<  core correction factor for kinetic energy density

  logical, intent(out)              ::  lres                             !<  lecture seems successful

! local variables

  character(len=2)     ::  ctype                                         !  type of calculation flag
  integer              ::  nc, nv
  integer              ::  ni, li                                        !  valence orbital
  integer              ::  nrc                                           !  number of core radii fields (4 or lmax+1)
  integer              ::  lmax                                          !  maximum angular momentum
  integer              ::  icfac                                         !  column where cfac starts
  real(REAL64)         ::  c1
  real(REAL64)         ::  rcin(0:8)                                     !  core radii as read (up to 9 fields)

  integer              ::  ioerror
  character(len=2)     ::  name                                          !  element on the atom.dat file

  character(len=128)   ::  cline

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  real(REAL64), parameter    ::  EPS = 1.0E-12_REAL64

! counters

  integer        ::  i, l


! default values

  do l = 0,MXDLC
    rc(l) = ZERO
  enddo
  cfac = ONE
  rcfac = ZERO
  cfactau = ONE

! read rc(s),rc(p),rc(d),rc(f),rc(g),cfac,rcfac or
! read rc(s),rc(p),rc(d), rc(f), cfac,rcfac or
! read rc(s),rc(p),rc(d), empty, cfac,rcfac or....

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
  lmax = 0
  do i = 1,nv
    read(ioread,'(2i5)') ni, li
    if(li > lmax) lmax = li
  enddo
  if(lmax > 3) nrc = lmax+1

  if(lmax > MXDLC .or. lmax > 8) then
    write(6,*)
    write(6,*) '    Stopped in atom_psd_read_atom, lmax = ',lmax,' > ',MXDLC, ' or 8'

    stop

  endif

! uses cline to be able to keep back-compatibility

  do i = 1,128
    cline(i:i) = ' '
  enddo

  read(ioread,'(a128)')  cline

  close(unit=ioread)

! core radii by angular momentum, stored up to MXDLC
! (atom_psd_ind stops if the valence has larger l)

  read(cline,'(9f10.5)',iostat = ioerror) (rcin(l),l=0,lmax)
  if(ioerror == 0) then
    do l = 0,lmax
      rc(l) = rcin(l)
    enddo
  endif

! cfac, rcfac and cfactau, each in its own 10 column field as in the
! 1982 code.  A blank field reads as zero and keeps the default,
! a field that cannot be read (e.g. a label) also keeps the default.

  icfac = 10*nrc+1

  read(cline(icfac:icfac+9),'(f10.5)',iostat = ioerror) c1
  if(ioerror == 0) cfac = c1
  read(cline(icfac+10:icfac+19),'(f10.5)',iostat = ioerror) c1
  if(ioerror == 0) rcfac = c1
  read(cline(icfac+20:icfac+29),'(f10.5)',iostat = ioerror) c1
  if(ioerror == 0) cfactau = c1

  lres = .TRUE.
  if(abs(rc(0)) < EPS) lres = .FALSE.

! negative cfac means full core, keep the option for cfactau

  if(abs(cfac) < EPS) cfac = ONE
  if(abs(cfactau) < EPS) cfactau = ONE

  return

end subroutine atom_psd_read_atom

