!>  Writes a file (default pseudokb.dat) with the pseudopotential
!>  in Kleinman Bylander form in real space
!>
!>  \author       Norm Troullier, J.L.Martins
!>  \version      6.0.8
!>  \date         90s, 13 August 2021, 19 May 2022.
!>  \copyright    GNU Public License v2

subroutine atom_kb_psd_out_real(iopsd, filepsd,                          &
      nameat, icorr, irel, nicore, irdate, irvers, irayps, psdtitle,     &
      npot, lo, nr, a, b, r, zion,                                       &
      vlocal, inorm, vkbproj, cdc, cdv,                                  &
      mxdl, mxdnr)

! adapted from the old program August 2021. JLM
! mxdnr, mxdl, 18 September 2021. JLM
! psdtitle, 19 May 2022. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum
  integer, intent(in)               ::  mxdnr                            !<  dimension of radial grid points.

  integer, intent(in)               ::  iopsd                            !<  io tape number
  character(len=*), intent(in)      ::  filepsd                          !<  file name of the output of the atomic program, parsec style

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbol of the element
  character(len=2), intent(in)      ::  icorr                            !<  correlation used in the calculation
  character(len=3), intent(in)      ::  irel                             !<  flag for relativistic (r) or spin-polarized (s) original calculation
  character(len=4), intent(in)      ::  nicore                           !<  flag for core correction
  character(len=10), intent(in)     ::  irdate, irvers                   !<  date and version of original calculation
  character(len=10), intent(in)     ::  irayps(4)                        !<  type of pseudopotential
  character(len=10), intent(in)     ::  psdtitle(20)                     !<  pseudopotential parameters

  integer, intent(in)               ::  npot(-1:1)                       !<  number of orbitals to be calculated
  integer, intent(in)               ::  lo(mxdl+1,-1:1)                  !<  angular momentum of orbital

  integer, intent(in)               ::  nr                               !<  number of points in the radial grid
  real(REAL64), intent(in)          ::  a,b                              !<  constants used to generate the radial grid
  real(REAL64), intent(in)          ::  r(mxdnr)                       !<  radial grid points   r(i) = a*(exp(b*i)-1), i=1,...,nr

  real(REAL64), intent(in)          ::  zion                             !<  ionic charge

  real(REAL64), intent(in)          ::  vlocal(mxdnr)                  !<  r*local pseudopotential
  integer, intent(in)               ::  inorm(0:mxdl,-1:1)               !<  sign of denominator of KB operator
  real(REAL64), intent(in)          ::  vkbproj(mxdnr,0:mxdl,-1:1)     !<  kb-projector
  real(REAL64), intent(in)          ::  cdc(mxdnr)                     !<  4*pi*r**2 charge density of core
  real(REAL64), intent(in)          ::  cdv(mxdnr)                     !<  4*pi*r**2 charge density of valence.  kb_conv is not spin polarized...

! local variables

  integer       ::  nrm

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64

! counters

  integer                  :: i, j


  if(iopsd == 0) RETURN

  nrm = nr-1

  open(unit = iopsd, file = filepsd, status = 'unknown', form = 'unformatted')

! hack that keeps compatibility with old versions

  write(iopsd) nameat, icorr, irel, nicore, irdate, irvers,              &
     (irayps(i),i=1,4), (psdtitle(i),i=1,7),                             &
     npot(0), npot(-1), nrm, a, b, zion, (psdtitle(i),i=8,20)

  write(iopsd) (r(i),i=2,nr)

! Write the potentials to the current pseudo.dat file (unit=1).

  if(irel == 'rel') then

    do i = 1,npot(1)
      write(iopsd) lo(i,1), (vkbproj(j,lo(i,1), 1),j=2,nr)
    enddo
    do i = 1,npot(-1)
      write(iopsd) lo(i,-1), (vkbproj(j,lo(i,-1),-1),j=2,nr)
    enddo

  else

    do i = 1,npot(0)
      write(iopsd) lo(i,0), (vkbproj(j,lo(i,0), 0),j=2,nr)
    enddo

  endif

! Write the charge densities to the current pseudo.dat file (unit=1).

!   if (ifcore == 0) then
!     write(iopsd) (ZERO,j=2,nr)
!   else
    write(iopsd) (cdc(j),j=2,nr)
!  endif

  write(iopsd) (cdv(j),j=2,nr)

! write for input2

  write(iopsd) (vlocal(i)/r(i),i=2,nr)

  if(irel == 'rel') then

    write(iopsd) npot(1)
    do i = 1,npot(1)
      write(iopsd) inorm(lo(i,1),1), UM
    enddo

  else

    write(iopsd) npot(0)
    do i = 1,npot(0)
      write(iopsd) inorm(lo(i,0),0), UM
    enddo

  endif

  close(unit = iopsd)

  return

end subroutine atom_kb_psd_out_real
