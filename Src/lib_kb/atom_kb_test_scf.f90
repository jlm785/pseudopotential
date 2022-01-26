!>  kb test calculation for a given atomic configuration
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.0.6
!>  \date         22 June 2021, 13 August 2021. 28 October 2021.
!>  \copyright    GNU Public License v2

subroutine atom_kb_test_scf(etotal,                                      &
         nameat, ititle, icorr, ispp_in,                                 &
         nval_in, ni, li, zd, zu, evd,                                   &
         iowrite, iopsd, filepsd,                                        &
         mxdnr, mxdorb, mxdl)

! translated to f90 from atm.f version 5.805
! mxdl. 18 September 2021. JLM
! number of iterations. 21 October 2021. JLM
! vhxc_orb in dsolv1. 28 October 2021. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input


  integer, intent(in)               ::  mxdorb                           !<  dimension of orbitals
  integer, intent(in)               ::  mxdnr                            !<  dimension of radial points
  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  iopsd                            !<  default tape for pseudopotential in parsec format
  character(len=*), intent(in)      ::  filepsd                          !<  name of default tape for reading pseudopotential in parsec format

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbol of the atom

  character(len=10), intent(in)     ::  ititle(5)                        !<  title of calculation

  character(len=2), intent(in)      ::  icorr                            !<  correlation type

  integer, intent(in)               ::  ni(mxdorb)                       !<  valence principal quantum number
  integer, intent(in)               ::  li(mxdorb)                       !<  valence angular quantum number
  real(REAL64), intent(in)          ::  zd(mxdorb), zu(mxdorb)           !<  occupation of valence down and up orbitals
  real(REAL64), intent(in)          ::  evd(mxdorb)                      !<  default eigenvalue

! modified locally, do not propagate the configuration!

  character(len=1), intent(in)      ::  ispp_in                          !<  spin polarization  ' ', 's', 'r'

  integer, intent(in)               ::  nval_in                          !<  number of valence orbitals

! output

  real(REAL64), intent(out)         ::  etotal                           !<  total energy

! used array dimensions

  integer                           ::  norb                             !  number of orbitals
  integer                           ::  nr                               !  number of radial points in mesh

! allocatable arrays

  real(REAL64), allocatable         ::  r(:)                             !  radial grid points
  real(REAL64), allocatable         ::  drdi(:)                          !  d r(i) / d i
  real(REAL64), allocatable         ::  d2rodr(:)                        !  (d^2 r(i) / d i^2) / (d r / d i)

  integer, allocatable              ::  no(:)                            !  principal quantum number n
  integer, allocatable              ::  lo(:)                            !  angular quantum number l
  integer, allocatable              ::  iso(:)                           !  2*spin or 2*(j-l)
  real(REAL64), allocatable         ::  zo(:)                            !  orbital occupation
  real(REAL64), allocatable         ::  evi(:)                           !  orbital default energy

  real(REAL64), allocatable         ::  cdc(:)                           !  core charge density (total)
  real(REAL64), allocatable         ::  cdv(:,:)                         !  valence charge density (total)

  real(REAL64), allocatable         ::  vhxc_orb(:,:)                     !  screening potential in input (Ry)

  real(REAL64), allocatable         ::  rpsi(:,:,:)                      !  wavefunctions (r(i),l,2j-2l).

  real(REAL64), allocatable         ::  vhxc_out(:,:)                    !  screening potential in Rydberg
  real(REAL64), allocatable         ::  vhxc_in(:,:)                     !  screening potential in Rydberg

  real(REAL64), allocatable         ::  ev(:)                            !  orbital energy
  real(REAL64), allocatable         ::  ek(:)                            !  orbital kinetic energy
  real(REAL64), allocatable         ::  ep(:)                            !  orbital potential energy

  integer, allocatable              ::  lpot(:,:)                        !  angular momentum l ot potential
  real(REAL64), allocatable         ::  vlocal(:)                        !  local pseudopotential
  integer , allocatable             ::  inorm(:,:)                       !  sign of denominator of KB operator
  real(REAL64), allocatable         ::  vkbproj(:,:,:)                   !  kb-projector

! allocatable arrays for dmixp

  real(REAL64), allocatable         ::  vmem(:,:,:)

! subroutine variables

  character(len=1)                  ::  ispp                             !  spin polarization  ' ', 's', 'r'
  character(len=3)                  ::  irel                             !  flag for relativistic (r) or spin-polarized (s) original calculation
  integer                           ::  ifcore                           !<  0 no partial core correction, 1 partial xc, 2 partial hartree
  character(len=10)                 ::  irdate, irvers                   !  date and version of original calculation
  character(len=10)                 ::  irayps(4)                        !  type of pseudopotential
  character(len=10)                 ::  ititle_kb(7)                     !  pseudopotential parameters

  real(REAL64)                      ::  znuc                             !  nuclear charge
  real(REAL64)                      ::  znuc_in                          !  nuclear charge

  integer                           ::  npot(-1:1)                       !  number of orbitals calculated
  real(REAL64)                      ::  zion                             !  ionic charge

  integer                           ::  nval                             !  number of valence orbitals

  integer                           ::  itype                            !  type of calculation (-1 signals end of calculation)

  real(REAL64)                      ::  a                                !  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64)                      ::  b                                !  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64)                      ::  zel                              !  electron charge

  real(REAL64)                      ::  etot(10)                         !  components of total energy

  integer                           ::  iter                             !  iteration number
  integer                           ::  iconv                            !  convergence flag (if iconv = 1, calculates Hartree energy)
  integer                           ::  maxit                            !  maximum number of iterations

! other variables

  real(REAL64)                      ::  t1, t2                           !  timings
  integer                           ::  itsm                             !  The screening potential is intially mixed with a percentage of old and new for the first itsm iterations
  integer                           ::  iiter

  real(REAL64)                      ::  xmixo                            !  potential mixing parameter
  real(REAL64)                      ::  dv, dvmax                        !  potential difference

  integer             ::  iz

! counters

  integer     ::  i, icon2, j

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  TOL = 1.0E-8_REAL64


  call zesec(t1)

! maximum number of iterations (bigger for KB)

  maxit = 300

! keep the configuration pristine

  call atom_p_tbl_charge(nameat, iz)

  znuc_in = iz*ONE

  ispp = ispp_in
  znuc = znuc_in
  nval = nval_in


  allocate(r(mxdnr))
  allocate(drdi(mxdnr))
  allocate(d2rodr(mxdnr))

  allocate(no(mxdorb))
  allocate(lo(mxdorb))
  allocate(iso(mxdorb))
  allocate(zo(mxdorb))
  allocate(evi(mxdorb))

  allocate(cdc(mxdnr))
  allocate(cdv(mxdnr,-1:1))

  allocate(vhxc_out(mxdnr,-1:1))
  allocate(vhxc_in(mxdnr,-1:1))

  allocate(vhxc_orb(mxdnr,mxdorb))

  allocate(ev(mxdorb))
  allocate(ek(mxdorb))
  allocate(ep(mxdorb))

  allocate(rpsi(mxdnr,0:mxdl,-1:1))

  allocate(vmem(mxdnr,4,-1:1))

  allocate(lpot(mxdl+1,-1:1))
  allocate(vlocal(mxdnr))
  allocate(inorm(0:mxdl,-1:1))
  allocate(vkbproj(mxdnr,0:mxdl,-1:1))

  call atom_kb_test_start(nameat, ititle, icorr, ispp,                   &
       nval, ni, li, zd, zu, evd,                                        &
       norb, no, lo, iso, zo, zel, evi,                                  &
       iowrite, mxdorb, mxdl)

! Set up the initial charge density.
! cdd and cdu  =  (4 pi r**2 rho(r))/2

! The charge density setup (aazn) is scaled with an empirical formula that reduces the
! number of iterations needed for the screening potential convergence.  njtj

  itsm = nint(znuc/9 + 3)

  call atom_kb_test_in_real(iopsd, filepsd,                              &
      nameat, icorr, irel, ifcore, irdate, irvers, irayps, ititle_kb,    &
      npot, lpot, nr, a, b, r, drdi, d2rodr, zion,                       &
      vlocal, inorm, vkbproj, cdc, cdv,                                  &
      mxdl, mxdnr)

  do i = 1,nr
    cdv(i,-1) = cdv(i, 0) / 2
    cdv(i, 1) = cdv(i, 0) / 2
  enddo


! Set up the electronic potential.

  call atom_atm_velect(0, 0, icorr, ispp, ifcore,                        &
      nr, r, drdi, zel, cdv, cdc, vhxc_out, etot,                        &
      iowrite, mxdnr)

  do j = -1,1
    do i = 1,nr
      vhxc_in(i,j) = vhxc_out(i,j)
    enddo
  enddo

! Start the iteration loop for electronic convergence.

  iconv = 0
  icon2 = 0

! The screening potential mixing parameter is an empirical function of the nuclear charge.
! Larger atoms require a slower convergence then smaller atoms.  njtj

  xmixo = ONE / log(znuc+7*ONE)

  call atom_kb_test_dsolv1(nr, a, b, r, drdi, norb,                      &
      no, lo, iso, zo, cdv, vhxc_out, ev,                                &
      vlocal, vkbproj, inorm, rpsi,                                      &
      iowrite, mxdnr, mxdorb, mxdl)

  do iter = 1,maxit

!   last iteration just compute orbitals

    if (iter == maxit) iconv = 2

!   compute orbitals


    do j = 1,norb
      do i = 1,nr
        vhxc_orb(i,j) = vhxc_in(i,iso(j))
      enddo
    enddo

    call atom_kb_test_dsolv2(iter, iconv,                                &
        nr, r, drdi, d2rodr,                                             &
        norb, no, lo, iso, zo,                                           &
        cdv, vhxc_orb, ev, ek, ep,                                        &
        vlocal, vkbproj, inorm, rpsi,                                    &
        iowrite, mxdnr, mxdorb, mxdl)

!    write(6,*)  '  finished dsolk2'
!    write(6,'("  ev(i) = ",10g10.3)') (ev(i),i=1,norb)

!   set up output electronic potential from charge density

    call atom_atm_velect(iter, iconv, icorr, ispp, ifcore,               &
        nr, r, drdi, zel, cdv, cdc, vhxc_out, etot,                          &
        iowrite, mxdnr)

!   check for convergence

    if (iconv > 0) exit

    dvmax = ZERO
    do j = -1,1
    do i = 1,nr
      dv = (vhxc_in(i,j)-vhxc_out(i,j)) / (ONE+vhxc_in(i,j)+vhxc_out(i,j))
      if (abs(dv) > dvmax) dvmax = abs(dv)
    enddo
    enddo
    icon2 = icon2 + 1
    if (dvmax <= TOL) iconv=1

!   Mix the input and output electronic potentials.

!   The screening potential is initially mixed with a percentage of old and new
!   for itsm iterations.  This brings the potential to a stable region after which
!   an Anderson's extrapolation scheme is used.  njtj

    if (iter < itsm) then
      iiter = 2
    else
      iiter = iter - itsm + 3
    endif

    do j = -1,1
      call atom_atm_dmixp(vhxc_out(:,j), vhxc_in(:,j), xmixo, iiter,     &
          3, nr, vmem(:,:,j), mxdnr)
    enddo

  enddo

  if(iconv == 2) then

!   End of iteration of loop without convergence.

    write(iowrite,*)
    write(iowrite,'(" potential not converged - dvmax =",                &
         & e12.4," xmixo =",f5.3 )') dvmax,xmixo

    stop

  endif

  write(iowrite,*)
  write(iowrite,'("Total number of iterations needed for",               &
     & " electron screening potential is ",i4)') icon2
  write(iowrite,*)

!  Find the total energy.  jlm

  call atom_atm_etotal(itype, nameat, norb,                              &
      no, lo, iso, zo, etot, ev, ek, ep,                                 &
      iowrite, mxdorb)

  etotal = etot(10)

  write(iowrite,*)
  write(iowrite,*)  '  Potential energy not calculated correctly '
  write(iowrite,*)  '  in test with KB projectors... '
  write(iowrite,*)

! prints charge density if fifth variable is not zero...

!  call atom_atm_rhoan(nr, r, cdd, cdu, 0, mxdnr)

! Writes file for later pseudopotential generation.


  call zesec(t2)

  write(iowrite,*)
  write(iowrite,'("  computing time for this configuration: ",f12.3)') t2-t1
  write(iowrite,*)

  deallocate(r)
  deallocate(drdi)
  deallocate(d2rodr)

  deallocate(no)
  deallocate(lo)
  deallocate(iso)
  deallocate(zo)
  deallocate(evi)

  deallocate(cdc)
  deallocate(cdv)

  deallocate(vhxc_out)
  deallocate(vhxc_in)

  deallocate(vhxc_orb)

  deallocate(ev)
  deallocate(ek)
  deallocate(ep)

  deallocate(vmem)

  deallocate(lpot)
  deallocate(vlocal)
  deallocate(inorm)
  deallocate(vkbproj)

  return

end subroutine atom_kb_test_scf
