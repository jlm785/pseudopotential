!>  Reads the pseudopotential from a file for testing
!>  Note that vionic is the ionic potential times r.
!>
!>  \author       Sverre Froyen, Norm Troullier, Jose Luis Martins
!>  \version      6.0.3
!>  \date         22 June 2021, 11 September 2021.
!>  \copyright    GNU Public License v2

 subroutine atom_atm_vpseudo(ispp, icorr, ifcore,                        &
     nr, a, b, r, drdi, d2rodr, nameat,                                  &
     cdv, cdc, vionic,                                                   &
     iopsd, filepsd, iowrite, mxdnr, mxdl)

! njtj ***  major modifications  ***
!    If a potential does not exist, it is approximated
!    by an existing potential.
!    A nonspin or spin-polarized pseudo test, uses the
!    down(nonspin generation), weighted average(spin-
!    polarized), or averaged(relativistic) potentials.
!    A relativistic pseudo test, must use relativistic
!    generated potentials.  The Schroedinger equation is
!    used to integrate a relativistic pseudo test,
!    not the Dirac equation.
! njtj  ***  major modifications  ***

! converted to f90, February 2018
! cleanup and new interface, July 2019. JLM
! adapted from vionic.f
! remove Coulomb and rsh, zsh, 6 August 2021. JLM
! vionic, cdd, 10 September 2021. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, intent(in)               ::  iopsd                            !<  default tape for pseudopotential
  character(len=*), intent(in)      ::  filepsd                          !<  name of default tape for reading pseudopotential in parsec format

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points
  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum l

! input

  character(len=2), intent(in)      ::  icorr                            !<  correlation type

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbot of the atom

! output

  integer, intent(out)              ::  ifcore                           !<  0 no partial core correction, 1 partial xc, 2 partial hartree

  real(REAL64), intent(out)         ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(out)         ::  drdi(mxdnr)                      !<  d r(i) / d i
  real(REAL64), intent(out)         ::  d2rodr(mxdnr)                    !<  (d^2 r(i) / d i^2) / (d r / d i)

  real(REAL64), intent(out)         ::  cdc(mxdnr)                       !<  4*pi*r**2 * core charge density

  real(REAL64), intent(out)         ::  vionic(mxdnr,0:mxdl,-1:1)        !<  r*ionic potential in Rydberg; -1:  j=l-1/2 or s= -1/2.  0:  average or spinless.  1:  j=l+1/2 or s=1/2

! input and output

  integer, intent(inout)            ::  nr                               !<  dimension of the number of radial points

  character(len=1), intent(inout)   ::  ispp                             !<  spin polarization  ' ', 's', 'r'
  real(REAL64), intent(inout)       ::  a                                !<  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64), intent(inout)       ::  b                                !<  r(i) = a*(exp(b*(i-1))-1)

  real(REAL64), intent(inout)       ::  cdv(mxdnr,-1:1)                  !<  4*pi*r**2 * valence charge density

! local variables

  integer, allocatable       ::  np(:,:)
  integer                    ::  npot(-1:1)

  character(len=2)           ::  icorrt, namet
  character(len=3)           ::  irel
  character(len=4)           ::  nicore
  character(len=10)          ::  iray(6), ititle(7)

  integer                    ::  nrm
  real(REAL64)               ::  zion
  integer                    ::  loi                                     !  angular momentum of orbital i

  real(REAL64)               ::  vsum, vdiff
  real(REAL64)               ::  rtry

  integer                    ::  lmax                                    !  maximum angular momentum

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  real(REAL64), parameter    ::  SMALL = 1.0E-8_REAL64

! counters

  integer                    ::  i, j, l


! read pseudopotentials from iopsd

  allocate(np(0:mxdl,-1:1))

  open(unit=iopsd, file=trim(filepsd), form='unformatted', status='unknown')

  read(iopsd) namet, icorrt, irel ,nicore, (iray(i),i=1,6),              &
      (ititle(i),i=1,7), npot(1), npot(-1), nrm, a, b, zion

  ifcore = 0
  if(nicore == 'fcec'.or.nicore == 'pcec') ifcore = 1
  if(nicore == 'fche'.or.nicore == 'pche') ifcore = 2

! compatibility with very old code (where r(1) /= 0)

  nr = nrm+1
  if(nr > mxdnr) then
    write(6,*)
    write(6,'("   stopped in atm_vpseudo, nr = ",i10," mxdnr = ",i10)') nr,mxdnr
    write(6,*)

    STOP

  endif

  read(iopsd) (r(i),i = 2,nr)
  r(1) = ZERO

! paranoid check

  do i = 1,nr
    rtry = a*(exp(b*(i-1))-ONE)

    if(abs(rtry - r(i)) > SMALL) then
      write(6,'("   INCOMPATIBLE  a, b r(i) ")')

      STOP

    endif

    drdi(i) = (r(i)+a)*b
    d2rodr(i) = b
  enddo

! down potentials (or average relativistic potentials)

! njtj  ***  major start  ***
!   if a potential does not exist, it is replaced by the
!   next existing lower angular momentum potential or
!   the next existing higher if no lower exist.

  do l = 0,mxdl
    np(l, 1) = 0
  enddo

  lmax = 0
  do i = 1,npot( 1)
    read(iopsd) loi,(vionic(j,loi, 1),j = 2,nr)
    if(loi > mxdl) then
      write(6,*)
      write(6,'("   stopped in atm_vpseudo, l = ",i10," mxdl = ",i10)') loi, mxdl
      write(6,*)

      STOP

    endif
    vionic(1,loi, 1) = ZERO
    np(loi, 1) = 1
    if(lmax < loi) lmax = loi
  enddo

! if there is no potential for l=0, replaces with lowest available pseudopotential

  if (np(0, 1) == 0) then

    do l = 1,lmax
      if (np(l, 1) > 0) then
        do j = 1,nr
          vionic(j,0, 1) = vionic(j,l, 1)
        enddo

        exit

      endif

    enddo

  endif

! if there is no potential for l, replaces with l-1

  do l = 1,mxdl
    if (np(l, 1) == 0) then
      do j = 1,nr
        vionic(j,l, 1) = vionic(j,l-1, 1)
      enddo
    endif
  enddo

! up potentials (or spin orbit potentials)

  lmax = 0
  if (npot(-1) > 0) then

    do l = 0,lmax
      np(l,-1) = 0
    enddo

    do i = 1,npot(-1)
      read(iopsd) loi,(vionic(j,loi,-1),j = 2,nr)
      if(loi > mxdl) then
        write(6,*)
        write(6,'("   stopped in atm_vpseudo, l = ",i10," mxdl = ",i10)') loi, mxdl
        write(6,*)

        STOP

      endif
      vionic(1,loi,-1) = ZERO
      np(loi,-1) = 1
      if(lmax < loi) lmax = loi
    enddo

!   if there is no potential for l=0, replaces with lowest available pseudopotential

    if (np(0,-1) == 0) then
      do l = 1,lmax
        if (np(l,-1) > 0) then
          do j = 1,nr
            vionic(j,0,-1) = vionic(j,l,-1)
          enddo

          exit

        endif

      enddo
    endif

!   if there is no potential for l, replaces with l-1

    do l = 1,mxdl
      if (np(l,-1) == 0) then
        do j = 1,nr
          vionic(j,l,-1) = vionic(j,l-1,-1)
        enddo
      endif
    enddo

  endif

! njtj  ***  major end  ***


! core and valence charges

  read(iopsd) (cdc(i),i = 2,nr)
  cdc(1) = ZERO

  read(iopsd) (cdv(i, 0),i = 2,nr)
  cdv(1, 0) = ZERO

  close(unit=iopsd)

! njtj  ***   major start  ***
!   distribute charge as up and down charge
!   generate radial intergration grid
!   set up potentials equal to down potentials for
!   spin-polarized pseudo test of nonspin and relativistic
!   generated potentails.  Construct spin-orbit potentials
!   from relativistic sum and difference potentials and
!   change ispp='r' to ispp=' '.

  do i = 1,nr
   cdv(i, 1) = cdv(i, 0) / 2
   cdv(i,-1) = cdv(i, 1)
  enddo

! non-relativistic, non-spin-polarized

  if(ispp /= 'r' .and. ispp /= 's') then
    do l = 0,mxdl
      do j = 1,nr
        vionic(j,l, 0) = vionic(j,l, 1)
      enddo
    enddo
  endif

! calculation is spin-polarized but pseudopotential isn't

  if (ispp == 's' .and. irel /= 'isp') then
    do l = 0,mxdl
      do j = 1,nr
        vionic(j,l,-1) = vionic(j,l,0)
      enddo
    enddo
  endif

! pseudopotentials calculations are never relativistic

  if (ispp == 'r') then
    ispp = ' '
    if (irel /= 'rel') then

      write(iowrite,*)
      write(iowrite,'("  Pseudopotentail is not relativistic!!!!")')
      write(iowrite,'("  setting potentials as equal.  ",a3)') irel

      do l = 0,mxdl
        do j = 1,nr
          vionic(j,l,-1) = vionic(j,l, 1)
        enddo
      enddo

    else

      do j = 1,nr
        vionic(j,0,-1) = vionic(j,0, 1)
        vionic(j,0, 0) = vionic(j,0, 1)
      enddo
      do l = 1,mxdl
        do j = 1,nr
          vsum = vionic(j,l, 1)
          vdiff = vionic(j,l,-1)
          vionic(j,l, 1) = vsum - ((l+1)*vdiff)/2
          vionic(j,l, 0) = vsum
          vionic(j,l,-1) = vsum + (l*vdiff)/2
        enddo
      enddo

    endif
  endif

! njtj  ***  major end   ***


! printout

  write(iowrite,'(//,1x,a2,2x,a2,2x,a3,2x,a4,"  pseudopotential",        &
      & " read from tape",/,1x,2a10,5x,4a10,/,1x,7a10,//)')              &
       namet,icorrt,irel,nicore,(iray(i),i=1,6),(ititle(i),i=1,7)

  if (nameat /= namet) write(iowrite,'(" input element ",a2,             &
      & " not equal to element on tape ",a2,//)') nameat,namet

  if (icorr /= icorrt) write(iowrite,'(" input correlation ",a2,         &
      & " not equal to correlation from tape ",a2,//)') icorr,icorrt

!jlm
  write(iowrite,'(" radial grid parameters",//,                          &
        & " r(1) = .0 , r(2) =",e11.3," , ... , r(",i4,") =",f8.3,//)')  &
        r(2),nr,r(nr)
!jlm

  deallocate(np)

  return

end subroutine atom_atm_vpseudo
