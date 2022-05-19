!>  Reads a file (default pseudokb.dat) with the pseudopotential
!>  in Kleinman Bylander form in real space
!>
!>  \author       Norm Troullier, J.L.Martins
!>  \version      6.0.8
!>  \date         13 August 2021, 19 May 2022.
!>  \copyright    GNU Public License v2

subroutine atom_kb_test_in_real(iopsd, filepsd,                          &
      nameat, icorr, irel, ifcore, irdate, irvers,                       &
      npot, lpot, nr, a, b, r, drdi, d2rodr, zion,                       &
      vlocal, inorm, vkbproj, cdc, cdv,                                  &
      mxdl, mxdnr)

! adapted from the old program August 2021. JLM
! mxdl, mxdnr, 17 September 2021. JLM
! irayps, psdtitle, 19 May 2022. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum
  integer, intent(in)               ::  mxdnr                            !<  dimension of radial grid points.

  integer, intent(in)               ::  iopsd                            !<  io tape number
  character(len=*), intent(in)      ::  filepsd                          !<  file name of the output of the atomic program, parsec style

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbol of the element
  character(len=2), intent(in)      ::  icorr                            !<  correlation used in the calculation

! output

  character(len=3), intent(out)     ::  irel                             !<  flag for relativistic (r) or spin-polarized (s) original calculation
  integer, intent(out)              ::  ifcore                           !<  0 no partial core correction, 1 partial xc, 2 partial hartree
  character(len=10), intent(out)    ::  irdate, irvers                   !<  date and version of original calculation

  integer, intent(out)              ::  npot(-1:1)                       !<  number of orbitals to be calculated
  integer, intent(out)              ::  lpot(mxdl+1,-1:1)                !<  angular momentum of potential

  integer, intent(out)              ::  nr                               !<  number of points in the radial grid
  real(REAL64), intent(out)         ::  a,b                              !<  constants used to generate the radial grid
  real(REAL64), intent(out)         ::  r(mxdnr)                         !<  radial grid points   r(i) = a*(exp(b*(i-1))-1), i=1,...,nr
  real(REAL64), intent(out)         ::  drdi(mxdnr)                      !<  d r(i) / d i
  real(REAL64), intent(out)         ::  d2rodr(mxdnr)                    !<  (d^2 r(i) / d i^2) / (d r / d i)

  real(REAL64), intent(out)         ::  zion                             !<  ionic charge

  real(REAL64), intent(out)         ::  vlocal(mxdnr)                    !<  local pseudopotential
  integer, intent(out)              ::  inorm(0:mxdl,-1:1)               !<  sign of denominator of KB operator
  real(REAL64), intent(out)         ::  vkbproj(mxdnr,0:mxdl,-1:1)       !<  kb-projector
  real(REAL64), intent(out)         ::  cdc(mxdnr)                       !<  4*pi*r**2 charge density of core
  real(REAL64), intent(out)         ::  cdv(mxdnr,-1:1)                  !<  4*pi*r**2 charge density of valence.  kb_conv is not spin polarized...
!  real(REAL64), intent(out)         ::  psi(mxdnr,0:mxdl,-1:1)          !<  wavefunctions (r(i),l,2j-2l).

! local variables

  real(REAL64)          ::  anorm
  real(REAL64)          ::  rtry
  character(len=4)      ::  nicore                                       !  flag for core correction
  character(len=2)      ::  icorrt, namet
  integer               ::  nrm

  character(len=10)     ::  irayps(4)                                    !  type of pseudopotential
  character(len=10)     ::  psdtitle(20)                                 !  pseudopotential parameters

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter    ::  SMALL = 1.0E-8_REAL64

! counters

  integer                  :: i, j


  open(unit = iopsd, file = filepsd, status = 'old', form = 'unformatted')

  read(iopsd) namet, icorrt, irel, nicore, irdate, irvers,               &
     (irayps(i),i=1,4), (psdtitle(i),i=1,7),                             &
     npot(0), npot(-1), nrm, a, b, zion

  if(namet /= nameat) then
    write(6,*)
    write(6,*) '  Stopped in atom_kb_test_in_real:'
    write(6,*) '  Elements from configuration and pseudopotential file'
    write(6,*) '  are different: ',namet,nameat

    STOP

  endif

  if(icorr /= icorrt) then
    write(6,*)
    write(6,*) '  WARNING in atom_kb_test_in_real:'
    write(6,*) '  Correlation from configuration and pseudopotential file'
    write(6,*) '  are different: ', icorr, icorrt
    write(6,*)
  endif

  ifcore = 0
  if(nicore == 'fcec'.or.nicore == 'pcec') ifcore = 1
  if(nicore == 'fche'.or.nicore == 'pche') ifcore = 2

  nr = nrm + 1
  if(nr > mxdnr) then
    write(6,*)
    write(6,'("   stopped in atm_kb_test_in_real, nr = ",i10," mxdnr = ",i10)') nr,mxdnr
    write(6,*)

    STOP

  endif

  npot(1) = npot(0)

  read(iopsd) (r(i),i=2,nr)
  r(1) = ZERO

! paranoid check

  do i = 1,nr
    rtry = a*(exp(b*(i-1))-1)

    if(abs(rtry - r(i)) > SMALL) then
      write(6,'("   INCOMPATIBLE  a, b r(i) ")')

      STOP

    endif

    drdi(i) = (r(i)+a)*b
    d2rodr(i) = b
  enddo

! Read the potentials from the current pseudokb.dat file

  if(irel == 'rel') then

    do i = 1,npot(1)
      read(iopsd) lpot(i, 1), (vkbproj(j,lpot(i, 1), 1),j=2,nr)
      vkbproj(1,lpot(i, 1), 1) = ZERO
    enddo
    do i = 1,npot(-1)
      read(iopsd) lpot(i,-1), (vkbproj(j,lpot(i,-1),-1),j=2,nr)
      vkbproj(1,lpot(i,-1),-1) = ZERO
    enddo

  else

    do i = 1,npot(0)
      read(iopsd) lpot(i, 0), (vkbproj(j,lpot(i, 0), 0),j=2,nr)
      vkbproj(1,lpot(i, 0), 0) = ZERO
    enddo

  endif

! read the charge densities from the current pseudokb.dat file

  read(iopsd) (cdc(j),j=2,nr)

  read(iopsd) (cdv(j,0),j=2,nr)

  cdc(1) = ZERO
  cdv(1,0) = ZERO

! write for input2

  read(iopsd) (vlocal(i),i=2,nr)
  vlocal(1) = vlocal(2)

  if(irel == 'rel') then

    read(iopsd) npot(1)
    do i = 1,npot(1)
      read(iopsd) inorm(lpot(i,1),1), anorm
      inorm(lpot(i,1),-1) = inorm(lpot(i,1),1)
    enddo

  else

    read(iopsd) npot(0)
    do i = 1,npot(0)
      read(iopsd) inorm(lpot(i,0),0)   !, anorm
    enddo

  endif

  close(unit = iopsd)

  return

end subroutine atom_kb_test_in_real
