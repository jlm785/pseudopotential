!>  atomic self-consistent calculation for a given atomic configuration
!>
!>  \author       Sverre Froyen, Norm Troullier, Jose Luis Martins
!>  \version      6.0.8
!>  \date         22 June 2021, 5 August 2021. 18 May 2022.
!>  \copyright    GNU Public License v2

subroutine atom_atm_scf(etotal,                                          &
         nameat, ctype, ititle, kerker, icorr, ispp_in,                  &
         znuc_in, rmax_in, aa_in, bb_in, ncore_in, nval_in,              &
         ni, li, zd, zu, evd,                                            &
         iowrite, ioae, fileae, iopsd, filepsd,                          &
         mxdnr, mxdorb, mxdl)

! translated to f90 from atm.f version 5.805
! vionic, so->iso, vhxc, 12 September 2021. JLM
! vhxc in atom_atm_datout. 18 October 2021. JLM
! deallocate vhxc_in. 18 May 2022. JLM




  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input


  integer, intent(in)               ::  mxdorb                           !<  dimension of the number of orbitals
  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points
  integer, intent(in)               ::  mxdl                             !<  dimension for angular momentum l

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  ioae                             !<  default tape for writing
  character(len=*), intent(in)      ::  fileae                           !<  name of default tape for all-electron results

  integer, intent(in)               ::  iopsd                            !<  default tape for pseudopotential in parsec format
  character(len=*), intent(in)      ::  filepsd                          !<  name of default tape for reading pseudopotential in parsec format

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbol of the atom

  character(len=2), intent(in)      ::  ctype                            !<  type of calculation flag (converted to itype)
  character(len=10), intent(in)     ::  ititle(5)                        !<  title of calculation

  character(len=3), intent(in)      ::  kerker                           !<  type of pseudopotential flag
  character(len=2), intent(in)      ::  icorr                            !<  correlation type

  integer, intent(in)               ::  ni(mxdorb)                       !<  valence principal quantum number
  integer, intent(in)               ::  li(mxdorb)                       !<  valence angular quantum number
  real(REAL64), intent(in)          ::  zd(mxdorb), zu(mxdorb)           !<  occupation of valence down and up orbitals
  real(REAL64), intent(in)          ::  evd(mxdorb)                      !<  default eigenvalue

! output

  real(REAL64), intent(out)         ::  etotal                           !<  total energy

! modified locally, do not propagate the configuration!

  character(len=1), intent(in)      ::  ispp_in                          !<  spin polarization  ' ', 's', 'r'

  real(REAL64), intent(in)          ::  znuc_in                          !<  nuclear charge

  real(REAL64), intent(in)          ::  rmax_in                          !<  maximum mesh radius
  real(REAL64), intent(in)          ::  aa_in                            !<  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64), intent(in)          ::  bb_in                            !<  a = exp(-aa)/znuc, b = 1/bb

  integer, intent(in)               ::  ncore_in                         !<  number of orbitals treated as core
  integer, intent(in)               ::  nval_in                          !<  number of valence orbitals

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

  real(REAL64), allocatable         ::  cdc(:)                           !  4*pi*r**2 * core charge density (not polarized)

  real(REAL64), allocatable         ::  vionic(:,:,:)                    !  r*ionic (pseudo-)potential in Rydberg; -1:  j=l-1/2 or s= -1/2.  0:  average or spinless.  1:  j=l+1/2 or s=1/2

  real(REAL64), allocatable         ::  cdv(:,:)                         !  4*pi*r**2 * valence charge density

  real(REAL64), allocatable         ::  vhxc_orb(:,:)                    !  screening potential in input (Ry)

  real(REAL64), allocatable         ::  ev(:)                            !  orbital energy
  real(REAL64), allocatable         ::  ek(:)                            !  orbital kinetic energy
  real(REAL64), allocatable         ::  ep(:)                            !  orbital potential energy

  real(REAL64), allocatable         ::  vhxc_out(:,:)                    !  screening potential (Ry)
  real(REAL64), allocatable         ::  vhxc_in(:,:)                     !  screening potential (Ry)

! allocatable arrays for dmixp

  real(REAL64), allocatable         ::  vmem(:,:,:)

! subroutine variables

  character(len=1)                  ::  ispp                             !  spin polarization  ' ', 's', 'r'

  real(REAL64)                      ::  znuc                             !  nuclear charge

  real(REAL64)                      ::  rmax                             !  maximum mesh radius
  real(REAL64)                      ::  aa                               !  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64)                      ::  bb                               !  a = exp(-aa)/znuc, b = 1/bb

  integer                           ::  ncore                            !  number of orbitals treated as core
  integer                           ::  nval                             !  number of valence orbitals

  integer                           ::  itype                            !  type of calculation (-1 signals end of calculation)
  integer                           ::  ikerk                            !  type of pseudopotential

  real(REAL64)                      ::  a                                !  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64)                      ::  b                                !  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64)                      ::  zel                              !  electron charge

  integer                           ::  ifcore                           !  0 no partial core correction, 1 partial xc, 2 partial hartree

  real(REAL64)                      ::  etot(10)                         !  components of total energy

  integer                           ::  iter                             !  iteration number
  integer                           ::  iconv                            !  convergence flag (if iconv = 1, calculates Hartree energy)
  integer                           ::  maxit                            !  maximum number of iterations

! other variables

  real(REAL64)                      ::  t1, t2                           !  timings
  integer                           ::  itsm                             !  The screening potential is intially mixed with a percentage of old and new for the first itsm iterations
  integer                           ::  iiter
  real(REAL64)                      ::  aazn

  real(REAL64)                      ::  xmixo                            !  potential mixing parameter
  real(REAL64)                      ::  dv, dvmax                        !  potential difference

! counters

  integer     ::  i, j, l, icon2

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  TOL = 1.0E-8_REAL64


  call zesec(t1)

! maximum number of iterations

  maxit = 100

! keep the configuration pristine

  ispp = ispp_in
  znuc = znuc_in
  rmax = rmax_in
  aa = aa_in
  bb = bb_in
  ncore = ncore_in
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

  allocate(vionic(mxdnr,0:mxdl,-1:1))

  allocate(vhxc_orb(mxdnr,mxdorb))

  allocate(ev(mxdorb))
  allocate(ek(mxdorb))
  allocate(ep(mxdorb))

  allocate(vmem(mxdnr,4,-1:1))

  allocate(vhxc_out(mxdnr,-1:1))
  allocate(vhxc_in(mxdnr,-1:1))

  call atom_atm_start(nameat, ctype, ititle, kerker, icorr, ispp,        &
       znuc, rmax, aa, bb, ncore, nval, ni, li, zd, zu, evd,             &
       itype, ikerk,  nr, a, b, r, drdi, d2rodr,                         &
       norb, no, lo, iso, zo, zel, evi,                                  &
       iowrite, mxdnr, mxdorb, mxdl)

  if (itype < 0) then

    write(iowrite,*)
    write(iowrite,*)
    write(iowrite,'("  Unrecognized type of calculation ")')

    stop

  endif

! Set up the initial charge density.
! cdv(1) and cdv(-1)  =  (4 pi r**2 rho(r))/2

! The charge density setup (aazn) is scaled with an empirical formula that reduces the
! number of iterations needed for the screening potential convergence.  njtj

  itsm = nint(znuc/9 + 3)

  if (itype < 4) then

    aazn = sqrt(sqrt(znuc))/2 + ONE
    do i = 1,nr
      cdv(i, 1) = zel * aazn**3 * exp(-aazn*r(i)) * r(i)**2 / 4
      cdv(i,-1) = cdv(i, 1)
      cdv(i, 0) = cdv(i, 1) + cdv(i,-1)
      cdc(i) = ZERO
    enddo

! set up Coulomb potentials

    do i = -1,1
    do l = 0,mxdl
    do j = 1,nr
      vionic(j,l,i) = -2*znuc
    enddo
    enddo
    enddo

    ifcore = 0

  else

   call atom_atm_vpseudo(ispp, icorr, ifcore,                            &
        nr, a, b, r, drdi, d2rodr, nameat,                               &
        cdv, cdc, vionic,                                                &
        iopsd, filepsd, iowrite, mxdnr, mxdl)

  endif


! Set up the screening potential (up and down).

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

  call atom_atm_dsolv1(nr, a, b, r, drdi, norb,                          &
        no, lo, iso, zo, cdv, vionic, vhxc_out, ev,                      &
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

    call atom_atm_dsolv2(iter, iconv, ispp, ifcore,                      &
          nr, r, drdi, d2rodr,                                           &
          norb, ncore, no, lo, iso, zo, znuc,                            &
          cdv, cdc, vionic, vhxc_orb, ev, ek, ep, evi,                   &
          iowrite, mxdnr, mxdorb, mxdl)

!   set up output screening potential from charge density

    call atom_atm_velect(iter, iconv, icorr, ispp, ifcore,               &
        nr, r, drdi, zel, cdv, cdc, vhxc_out, etot,                      &
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
     & " electron screening potential is ",i2)') icon2
  write(iowrite,*)

!  Find the total energy.  jlm

  call atom_atm_etotal(itype, nameat, norb,                              &
      no, lo, iso, zo, etot, ev, ek, ep,                                 &
      iowrite, mxdorb)

  etotal = etot(10)

! prints charge density if fifth variable is not zero...

!  call atom_atm_rhoan(nr, r, cdv, 0, mxdnr)

! Writes file for later pseudopotential generation.

  if(ioae /= 0) then

    call atom_atm_datout(ioae, fileae, itype, icorr, ispp,               &
        nr, a, b, r, drdi,                                               &
        nameat, norb, ncore, no, lo, iso, zo, znuc, zel,                 &
        cdc, vionic, vhxc_in, ev,                                        &
        mxdnr, mxdorb, mxdl)

  endif


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
  deallocate(vionic)

  deallocate(vhxc_orb)

  deallocate(ev)
  deallocate(ek)
  deallocate(ep)

  deallocate(vmem)

  return

end subroutine atom_atm_scf
