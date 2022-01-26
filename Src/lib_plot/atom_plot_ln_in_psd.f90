!>  Reads a file (default pseudo.dat or pseudokb.dat) with the pseudopotential
!>  (semi-local or KB non-local)
!>
!>  \author       Norm Troullier, J.L.Martins
!>  \version      6.0.3
!>  \date         22 August 2021. 25 September 2021.
!>  \copyright    GNU Public License v2

subroutine atom_plot_ln_in_psd(iopsd, filepsd, lkb,                      &
      nameat, icorr, irel, ifcore, iray, ititle,                         &
      norb, lo, iso, nr, r, drdi, d2rodr, znuc,                          &
      vionic, cdc, cdv, vlocal, vkbproj, inorm,                          &
      mxdorb, mxdl, mxdnr)

! adapted from the old program August 2021. JLM
! unaverage vionic. 25 September 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdorb                           !<  dimension of the number of orbitals
  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum
  integer, intent(in)               ::  mxdnr                            !<  dimension of radial grid points.

  integer, intent(in)               ::  iopsd                            !<  io tape number
  character(len=*), intent(in)      ::  filepsd                          !<  file name of the output of the atomic program, parsec style

  logical, intent(in)               ::  lkb                              !<  KB pseudopotential or semi-local

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbol of the element
  character(len=2), intent(in)      ::  icorr                            !<  correlation used in the calculation

! output

  character(len=3), intent(out)     ::  irel                             !<  flag for relativistic (r) or spin-polarized (s) original calculation
  integer, intent(out)              ::  ifcore                           !<  0 no partial core correction, 1 partial xc, 2 partial hartree
  character(len=10), intent(out)    ::  iray(6)                          !<  date version and type of pseudopotential
  character(len=10), intent(out)    ::  ititle(7)                        !<  pseudopotential parameters

  integer, intent(out)              ::  norb                             !<  number of orbitals to be calculated
  integer, intent(out)              ::  lo(mxdorb)                       !<  angular momentum of orbital
  integer, intent(out)              ::  iso(mxdorb)                      !<  2*spin of 2*(j-l)

  integer, intent(out)              ::  nr                               !<  number of points in the radial grid
  real(REAL64), intent(out)         ::  r(mxdnr)                         !<  radial grid points   r(i) = a*(exp(b*i)-1), i=1,...,nr
  real(REAL64), intent(out)         ::  drdi(mxdnr)                      !<  d r(i) / d i
  real(REAL64), intent(out)         ::  d2rodr(mxdnr)                    !<  (d^2 r(i) / d i^2) / (d r / d i)

  real(REAL64), intent(out)         ::  znuc                             !<  ionic charge

  real(REAL64), intent(out)         ::  vionic(mxdnr,0:mxdl,-1:1)        !<  r*ionic potential in Rydberg

  real(REAL64), intent(out)         ::  cdc(mxdnr)                       !<  4*pi*r**2 charge density of core
  real(REAL64), intent(out)         ::  cdv(mxdnr,-1:1)                  !<  4*pi*r**2 charge density of valence.

  real(REAL64), intent(out)         ::  vlocal(mxdnr)                    !<  local potentia in KB pseudopotentials
  real(REAL64), intent(out)         ::  vkbproj(mxdnr,0:mxdl,-1:1)       !<  kb-projector
  integer, intent(out)              ::  inorm(0:mxdl,-1:1)               !<  sign of denominator of KB operator


! local variables

  real(REAL64)          ::  anorm                                        !  anorm = 1 in new version.

  real(REAL64)          ::  vsum, vdiff
  real(REAL64)          ::  rtry
  character(len=4)      ::  nicore                                       !  flag for core correction
  character(len=2)      ::  icorrt, namet

  real(REAL64)          ::  a, b                                         !  constants used to generate the radial grid

  integer               ::  numu
  integer               ::  nrm
  integer               ::  norbs
  integer               ::  id

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter    ::  SMALL = 1.0E-8_REAL64

! counters

  integer                  :: i, j, l


  open(unit = iopsd, file = trim(filepsd), status = 'old', form = 'unformatted')

  read(iopsd) namet, icorrt, irel, nicore, (iray(i),i=1,6),              &
     (ititle(i),i=1,7), norb, numu, nrm, a, b, znuc

  nr = nrm + 1

  if(namet /= nameat) then
    write(6,*)
    write(6,*) '  Stopped in atom_plot_ln_in_psd:'
    write(6,*) '  Elements from configuration and pseudopotential file'
    write(6,*) '  are different: ',namet,nameat

    STOP

  endif

  if(icorr /= icorrt) then
    write(6,*)
    write(6,*) '  WARNING in atom_plot_ln_in_psd:'
    write(6,*) '  Correlation from configuration and pseudopotential file'
    write(6,*) '  are different: ', icorr, icorrt
    write(6,*)
  endif

  ifcore = 0
  if(nicore == 'fcec'.or.nicore == 'pcec') ifcore = 1
  if(nicore == 'fche'.or.nicore == 'pche') ifcore = 2

  if(nr > mxdnr) then
    write(6,*)
    write(6,*) '   Stopped in atom_plot_ln_in_psd:'
    write(6,'("   nr = ",i10," greater than mxdnr = ",i10)') nr,mxdnr
    write(6,*)

    STOP

  endif

  read(iopsd) (r(i),i=2,nr)
  r(1) = ZERO

! paranoid check

  do i = 1,nr
    rtry = a*(exp(b*(i-1))-1)

    if(abs(rtry - r(i)) > SMALL) then
      write(6,*)
      write(6,'("   Stopped in atom_plot_ln_in_psd")')
      write(6,'("   INCOMPATIBLE  a, b r(i) ")')

      STOP

    endif

    drdi(i) = (r(i)+a)*b
    d2rodr(i) = b
  enddo

  if(irel == 'nrl') then
    id = 0
  else
    id = 1
  endif

  if(numu /= 0 .and. irel == 'nrl') then
    write(6,*)
    write(6,*) '  Stopped in atom_plot_ln_in_psd'
    write(6,*) '  Scalar potential has spin components'
    write(6,*)

    STOP

  endif

  if(lkb) then

!   Read the KB projectors from the current pseudokb.dat file


    do i = 1,norb
      read(iopsd) lo(i), (vkbproj(j,lo(i),id),j=2,nr)
      vkbproj(1,lo(i),id) = ZERO
      iso(i) = id
    enddo

    if(numu > 0) then
      do i = norb+1,norb+numu
        read(iopsd) lo(i), (vkbproj(j,lo(i),-1),j=2,nr)
        vkbproj(1,lo(i),-1) = ZERO
        iso(i) = -1
      enddo
    endif

  else

!   Read the semi-local potentials from the current pseudo.dat file

    do i = 1,norb
      read(iopsd) lo(i), (vionic(j,lo(i),id),j=2,nr)
      vionic(1,lo(i),id) = ZERO
      iso(i) = id
    enddo

    if(numu > 0) then
      do i = norb+1,norb+numu
        read(iopsd) lo(i), (vionic(j,lo(i),-1),j=2,nr)
        vionic(1,lo(i),-1) = ZERO
        iso(i) = -1
      enddo
    endif

  endif

  read(iopsd) (cdc(j),j=2,nr)

  read(iopsd) (cdv(j, 0),j=2,nr)

  cdc(1) = ZERO
  cdv(1, 0) = ZERO
  do j = 1,nr
    cdv(j,-1) = cdv(j, 0) / 2
    cdv(j, 1) = cdv(j, 0) / 2
  enddo

  if(lkb) then
    read(iopsd) (vlocal(j),j=2,nr)
    vlocal(1) = vlocal(2)
    read(iopsd) norbs

!   In recent implementations anorm = 1.   Assumes inorm is the same for all "spins"

    do i=1,norbs
      read(iopsd) inorm(lo(i), 0),anorm
      inorm(lo(i), 1) = inorm(lo(i), 0)
      inorm(lo(i),-1) = inorm(lo(i), 0)
    enddo
  endif

  close(unit = iopsd)

  norb = norb + numu

! recovers pseudopotential

  if(irel == 'rel') then

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

  return

end subroutine atom_plot_ln_in_psd
