!>  Reads the pseudopotential core radii
!>
!>  \author       Jose Luis Martins
!>  \version      6.014
!>  \date         27 August 2021
!>  \copyright    GNU Public License v2

subroutine atom_psd_read_atom(ioread, filein, nameat, rc, cfac, rcfac, lres)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter                ::  lc = 4                           !  maximum valence angular momentum

! input

  integer, intent(in)               ::  ioread                           !  default tape for reading
  character(len=*), intent(in)      ::  filein                           !  name of default tape for reading

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbot of the atom

! output

  real(REAL64), intent(out)         ::  rc(0:lc)                         !<  core radius r_c(l)
  real(REAL64), intent(out)         ::  cfac                             !<  criteria for pseudo-core charge
  real(REAL64), intent(out)         ::  rcfac                            !<  pseudo-core radius

  logical, intent(out)              ::  lres                             !<  lecture seems successful

! local variables

  character(len=2)  ::  ctype                                            !  type of calculation flag
  integer           ::  nc, nv

  integer           ::  ioerror
  character(len=2)  ::  name                                             !  element on the atom.dat file

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64

! counters

  integer        ::  i, l


! default values

  do l = 0,4
    rc(l) = ZERO
  enddo
  cfac = ONE
  rcfac = ZERO

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
  read(ioread,'(7f10.5)', iostat = ioerror) (rc(l),l=0,4), cfac, rcfac

  close(unit=ioread)

  if (cfac == ZERO) cfac = ONE

  lres = .TRUE.
  if(ioerror /= 0) lres = .FALSE.
  if(abs(rc(1)) < 0.001) lres = .FALSE.

  return

end subroutine atom_psd_read_atom

