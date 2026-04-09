!>  Reads the pseudopotential core radii
!>
!>  \author       Jose Luis Martins
!>  \version      6.1.0
!>  \date         27 August 2021, 25 March 2026.
!>  \copyright    GNU Public License v2

subroutine atom_psd_read_atom(ioread, filein, nameat, rc,                &
      cfac, rcfac, cfactau, lres)

! written 27 August 2021 from previous code. JLM
! cfactau, 25 March 2026. JLM

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

  integer              ::  ioerror
  character(len=2)     ::  name                                          !  element on the atom.dat file

  character(len=128)   ::  cline
  real(REAL64)         ::  rctmp(0:mxdlc)

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64

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
  do i = 1,nv
    read(ioread,*)
  enddo

! uses cline to be able to keep back-compatibility

  do i = 1,128
    cline(i:i) = ' '
  enddo

  read(ioread,'(a128)')  cline

  close(unit=ioread)

! note that mxdlc is hard coded here, but it is the old format...

  read(cline,'(5f10.5)') rc(0)
  do l = 1,4
    read(cline,'(5f10.5)',iostat = ioerror) (rctmp(i),i=0,l-1),rc(l)
    if(ioerror /= 0) exit
  enddo

  read(cline(51:128),*,iostat = ioerror) cfac, rcfac

  lres = .TRUE.
  if(ioerror /= 0) lres = .FALSE.
  if(abs(rc(0)) < 0.001) lres = .FALSE.

  if(ioerror == 0) then
    read(cline(51:128),*,iostat = ioerror) rctmp(1), rctmp(2), cfactau
  endif

  if(cfac == ZERO) cfac = ONE
  if(cfactau == ZERO) cfactau = ONE

  return

end subroutine atom_psd_read_atom

