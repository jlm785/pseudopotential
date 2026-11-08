!>  reads an output file that is compatible with the parsec and siesta codes
!>
!>  \author       J.L.Martins
!>  \version      6.1.0
!>  \date         25 May 2012, 8 November 2026.
!>  \copyright    GNU Public License v2

subroutine atom_kb_psd_in_real(iotape, fname,                            &
    nameat, icorr, irel, nicore, irdate, irvers, irayps, psdtitle,       &
    npot, nr, a, b, r, zion, lo, vionic, cdc, cdv, tauc_partial,         &
    zo, rc, rpsi_ps,                                                     &
    lmax, mxdnr, mxdl)

! Written JLMartins
! it is based on FDP modifications in parsec's fork of the code.
! mxdl, mxdnr, 17 Sepember 2021. JLM
! psdtitle, 19 May 2022. JLM
! tauc, 7 January 2026. JLM
! in_parsec -> in_real in name. 26 January 2026. JLM
! tauc, 5 february 2026. JLM
! paranoid check of r against a, b (as other readers). 25 October 2026. claude
! list-directed reads, 4 numbers per line, of old (g20.12) and new (es24.16e3) psd.pot. 8 November 2026. claude+JLM

  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)

! input:

  integer, intent(in)               ::  mxdnr                            !<  dimension of radial grid points.
  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum

  integer, intent(in)               ::  iotape                           !<  io tape number
  character(len=*), intent(in)      ::  fname                            !<  file name of the output of the atomic program, in real (parsec) style

! output:

  character(len=2), intent(out)     ::  nameat                           !<  chemical symbol of the element
  character(len=2), intent(out)     ::  icorr                            !<  correlation used in the calculation
                                                                         !<      ca -> Ceperley-Alder LDA parametrized by Perdew-Zunger
                                                                         !<      pw -> Ceperley-Alder LDA parametrized by Perdew-Wang 92
                                                                         !<      pb -> The pbe GGA
  character(len=3), intent(out)     ::  irel                             !<  flag for relativistic (r) or spin-polarized (s) original calculation
                                                                         !<      nrl -> non-relativistic and non-spin-polarized
                                                                         !<      rel -> relativistic (includes spin-orbid)
                                                                         !<      isp -> spin-polarized non-relativistic (no spin-orbit)
  character(len=4), intent(out)     ::  nicore                           !<  flag for core correction
                                                                         !<      nc   -> no core correction
                                                                         !<      pcec -> partial core correction
                                                                         !<     fcec -> full core correction

  character(len=10), intent(out)    ::  irdate,irvers                    !<  date and version of original calculation
  character(len=10), intent(out)    ::  irayps(4)                        !<  type of pseudopotential
  character(len=10), intent(out)    ::  psdtitle(20)                     !<  pseudopotential parameters

  integer, intent(out)              ::  npot(-1:1)                       !<  number of orbitals (s,p,d,...).  -1:  j=l-1/2.  0:  average.  1:  j=l+1/2
                                                                         !<  meaning it depends on irel
  integer, intent(out)              ::  nr                               !<  number of points in the radial grid
  real(REAL64), intent(out)         ::  a,b                              !<  constants used to generate the radial grid
  real(REAL64), intent(out)         ::  r(mxdnr)                         !<  radial grid points   r(i) = a*(exp(b*i)-1), i=0,...,nr
  real(REAL64), intent(out)         ::  zion                             !<  ionic charge
  integer, intent(out)              ::  lo(1:mxdl+1,-1:1)                !<  angular momentum of orbital for each spin

  real(REAL64), intent(out)         ::  vionic(mxdnr,0:mxdl,-1:1)        !<  r*pseudopotential.  -1:  j=l-1/2.  0:  average.  1:  j=l+1/2
  real(REAL64), intent(out)         ::  cdc(mxdnr)                       !<  4*pi*r**2 charge density of core
  real(REAL64), intent(out)         ::  cdv(mxdnr)                       !<  4*pi*r**2 charge density of valence.
  real(REAL64), intent(out)         ::  tauc_partial(mxdnr)              !<  core kinetic energy density

  real(REAL64), intent(out)         ::  zo(0:mxdl)                       !<  orbital occupation (negative if not present!)
  real(REAL64), intent(out)         ::  rc(0:mxdl)                       !<  orbital occupation
  real(REAL64), intent(out)         ::  rpsi_ps(mxdnr,0:mxdl)            !<  r*wave-function

  integer, intent(out)              ::  lmax                             !<  maximum angular momentum for potential

! local variables

  integer                       ::  ifcore
  real(REAL64)                  ::  cfac, rcfac, rv0, rv0pr
  real(REAL64)                  ::  rtry
  integer                       ::  l
  integer                       ::  ios
  integer                       ::  nrm
  character(len=256)            ::  line                             !  line with nl nls nr a b zion (cfac rcfac)



! parameters

  real(REAL64), parameter                 :: ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  real(REAL64), parameter                 :: SMALL = 1.0E-8_REAL64

! counters

  integer                        :: i, j



  open(unit = iotape, file = trim(fname), status = 'UNKNOWN',form = 'FORMATTED')

  read(iotape,'(1x,a2,1x,a2,1x,a3,1x,a4)')                               &
                      nameat, icorr, irel, nicore
  do i = 1,20
    psdtitle(i) = '          '
  enddo
  read(iotape,'(1x,2a10,4a10,/,1x,20a10)', iostat = ios) irvers, irdate, &
          irayps, psdtitle
! keep compatibility with old files
  if(ios /= 0) then
    backspace iotape
    read(iotape,'(1x,2a10,4a10,/,1x,20a10)') irvers, irdate, irayps,     &
         (psdtitle(i),i=1,7)
  endif

  read(iotape,'(a)') line
  read(line,'(1x,2i3,i5)') npot(1), npot(-1), nrm
  ifcore = 1
  read(line(13:),*,iostat=ios) a, b, zion, cfac, rcfac
  if(ios /= 0) then
    ifcore = 0
    read(line(13:),*) a, b, zion
  endif

  nr = nrm+1
  if(nr > mxdnr) then
    write(6,'("  CHANGE DIMENSION OF MXDNR FROM ",I6," TO ",I6)')        &
               mxdnr,nr

    stop

  endif

  r(1) = ZERO
  read(iotape,*)                               !'Radial grid follows'

  call read4(r(2:nr))

! paranoid check

  do i = 1,nr
    rtry = a*(exp(b*(i-1))-ONE)

    if(abs(rtry - r(i)) > SMALL) then
      write(6,*)
      write(6,'("   Stopped in atom_kb_psd_in_real")')
      write(6,'("   INCOMPATIBLE  a, b r(i) ")')

      STOP

    endif
  enddo

! convention from atomic program:   npot(-1) > npot(1)

  npot(0) = npot(1)

  lmax = 0
  do i = 1, npot(0)
    read (iotape,*)                           !'Pseudopotential follows (l on next line)'
    read(iotape,'(1x,i2)') l
    lo(i,1) = l
    lo(i,0) = l
    if(lmax < l) lmax = l
    if(l > mxdl .or. l < 0) then
      write(6,'("  CHANGE DIMENSION OF MXDL FROM ",I3," TO ",I3)') mxdl,l

      stop

    endif

    call read4(vionic(2:nr,l,0))

    rv0pr = (vionic(3,l,0) - vionic(2,l,0))/(r(3)-r(2))
    rv0 = vionic(2,l,0) - rv0pr*r(2)
    if(abs(rv0 - nint(rv0)) < 1.0e-6) rv0 = nint(rv0)
    vionic(1,l,0) = rv0
  enddo

  do i = 1, npot(-1)
    read(iotape,*)                          !   'Spin-orbit. of Pseudop. follows (l on next line)'
    read(iotape,'(1x,i2)') l
    lo(i,-1) = l
    if(lmax < l) lmax = l
    if(l > mxdl .or. l < 0) then
      write(6,'("  CHANGE DIMENSION OF MXDL FROM ",I3," TO ",I3)') mxdl,l

      stop

    endif

    call read4(vionic(2:nr,l,-1))

    rv0pr = (vionic(3,l,-1) - vionic(2,l,-1)) / (r(3)-r(2))
    rv0 = vionic(2,l,-1) - rv0pr*r(2)
    if(abs(rv0 - nint(rv0)) < 1.0e-6) rv0 = nint(rv0)
    vionic(1,l,-1) = rv0
  enddo

! undo storage in average and spin-orbit

  if(lmax > 0) then
    do l = 1,lmax
      do j = 1,nr
        vionic(j,l,1)  = vionic(j,l,0) + (l*vionic(j,l,-1))/2
        vionic(j,l,-1) = vionic(j,l,0) - ((l+1)*vionic(j,l,-1))/2
      enddo
    enddo
  endif
  do j = 1,nr
    vionic(j,0,1)  = vionic(j,0,0)
  enddo

  read(iotape,*) ! 'Core charge follows'

  call read4(cdc(2:nr))
  cdc(1) = ZERO

  read(iotape,*) !'Valence charge follows'

  call read4(cdv(2:nr))
  cdv(1) = ZERO

! iotape still has information about wave-functions

  do l = 0,lmax
    zo(l) = -ONE
  enddo

  do i = 1, npot(0)
    read (iotape,*)                           ! Pseudo-wave-function follows (l, zelect, rc)
    read(iotape,'(1x,i2,1x,g20.12,1x,g20.12)') l, zo(l), rc(l)
    call read4(rpsi_ps(2:nr,l))
    rpsi_ps(1,l) = ZERO
  enddo

! reads the core kinetic energy

  read(iotape,*)                             !Core kinetic energy density follows
  call read4(tauc_partial(2:nr))
  tauc_partial(1) = ZERO


  close(unit=iotape)

  return

  contains

! reads an array written 4 numbers per line, list-directed, so that the
! old g20.12 and the new es24.16e3 files are both read

  subroutine read4(f)

    real(REAL64), intent(out)       ::  f(:)

    integer   ::  k, m

    do k = 1,size(f),4
      read(iotape,*) (f(m), m = k,min(k+3,size(f)))
    enddo

  end subroutine read4

end subroutine atom_kb_psd_in_real
