!>  writes the pseudopotential in the PSML format (version 1.1),
!>  semilocal potentials, local potential, KB projectors,
!>  pseudo wave-functions, valence and core charges.
!>
!>  \author       claude+J.L.Martins
!>  \version      6.1.0
!>  \date         28 October 2026.
!>  \copyright    GNU Public License v2

subroutine atom_kb_psd_out_psml(iotape, fname,                           &
      nameat, icorr, irel, nicore, irdate, irvers, irayps, psdtitle,     &
      nr, r, drdi, zion, npot, lo, vionic, rc,                           &
      vlocal, llocal, vkbproj, inorm, ev, rpsi, zo, cdc, cdv,            &
      mxdl, mxdnr)

! written from atom_kb_psd_out_psp8. 28 October 2026. claude

! Conventions of PSML (A. Garcia et al., Comput. Phys. Commun. 227, 51 (2018)
! and the libPSML documentation), checked with the files of libPSML:
!
!   energies in Hartree, lengths in bohr, the grid starts at r = 0,
!   radial functions without extra factors of r:
!     potentials V(r),
!     projectors chi(r) with integral r^2 chi^2 dr = 1, operator sum |chi> ekb <chi|,
!     pseudo wave-functions R(r) with integral r^2 R^2 dr = 1,
!     charges q(r) = 4 pi r^2 rho / r^2 with integral r^2 q dr = Q.
!
! Relativistic case (irel = 'rel'), with V+, V- the j = l+1/2, l-1/2 components:
!   set "lj":                   V+, V- (and projectors, wave-functions)
!   set "scalar_relativistic":  ((l+1) V+ + l V-) / (2l+1)
!   set "spin_orbit":           2 (V+ - V-) / (2l+1)
! the same weights as the .psf files.  The projector sets use the
! j projectors with the KB energies multiplied by (l+1)/(2l+1), l/(2l+1)
! (scalar) and 2/(2l+1), -2/(2l+1) (spin-orbit), as in the psp8 file.
! SIESTA reads the scalar_relativistic and spin_orbit sets (or lj
! semilocal potentials), ABINIT the scalar_relativistic and spin_orbit
! projectors.

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! dimensions:

  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum components
  integer, intent(in)               ::  mxdnr                            !<  dimension of radial grid points.

! input:

  integer, intent(in)               ::  iotape                           !<  output tape number
  character(len=*), intent(in)      ::  fname                            !<  file name of psml potential

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbol of the element
  character(len=2), intent(in)      ::  icorr                            !<  correlation used in the calculation
  character(len=3), intent(in)      ::  irel                             !<  flag for relativistic (rel), non-relativistic (nrl) or spin-polarized (isp)
  character(len=4), intent(in)      ::  nicore                           !<  flag for core correction
  character(len=10), intent(in)     ::  irdate                           !<  date of generation of the semi-local pseudopotential
  character(len=10), intent(in)     ::  irvers                           !<  version of original calculation
  character(len=10), intent(in)     ::  irayps(4)                        !<  type of pseudopotential
  character(len=10), intent(in)     ::  psdtitle(20)                     !<  pseudopotential parameters

  integer, intent(in)               ::  nr                               !<  number of points in the radial grid
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points, r(1) = 0
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  real(REAL64), intent(in)          ::  zion                             !<  ionic charge

  integer, intent(in)               ::  npot(-1:1)                       !<  number of potentials (0: average)
  integer, intent(in)               ::  lo(mxdl+1,-1:1)                  !<  angular momentum of the potentials
  real(REAL64), intent(in)          ::  vionic(mxdnr,0:mxdl,-1:1)        !<  r*pseudopotential (Ry). -1: j=l-1/2. 0: average. 1: j=l+1/2
  real(REAL64), intent(in)          ::  rc(0:mxdl)                       !<  core radii

  real(REAL64), intent(in)          ::  vlocal(mxdnr)                    !<  r*local pseudopotential (Ry)
  integer, intent(in)               ::  llocal                           !<  angular momentum of local pseudopotential (<0 supremum)
  real(REAL64), intent(in)          ::  vkbproj(mxdnr,0:mxdl,-1:1)       !<  kb-projector(r(i),l,2j-2l)
  integer, intent(in)               ::  inorm(0:mxdl,-1:1)               !<  sign of denominator of KB operator

  real(REAL64), intent(in)          ::  ev(0:mxdl,-1:1)                  !<  eigenvalues (l,2j-2l) (Ry)
  real(REAL64), intent(in)          ::  rpsi(mxdnr,0:mxdl,-1:1)          !<  r*wavefunctions (r(i),l,2j-2l)
  real(REAL64), intent(in)          ::  zo(0:mxdl)                       !<  orbital occupation (negative if absent)

  real(REAL64), intent(in)          ::  cdc(mxdnr)                       !<  4*pi*r**2 charge density of core
  real(REAL64), intent(in)          ::  cdv(mxdnr)                       !<  4*pi*r**2 charge density of valence

! local variables

  integer                 ::  iz
  logical                 ::  lso                                        !  spin-orbit (irel == 'rel')
  logical                 ::  lcore                                      !  partial core correction
  integer                 ::  lmax                                       !  maximum angular momentum of the potentials
  integer                 ::  nshell
  real(REAL64)            ::  zval                                       !  total valence occupation
  real(REAL64)            ::  qval                                       !  integral of the valence charge
  real(REAL64)            ::  pnorm, ekb, fsr, fso, xj
  integer                 ::  n, nseq
  integer                 ::  ierr

  character(len=8)        ::  date
  character(len=10)       ::  time
  character(len=36)       ::  uuid
  character(len=20)       ::  cset
  character(len=40)       ::  ltype

  integer, allocatable    ::  nq(:)                                      !  principal quantum number of each l (from psdtitle)
  logical, allocatable    ::  lhas(:)                                    !  potential l exists
  real(REAL64), allocatable  ::  f(:)                                    !  work array

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  character(len=1), parameter  ::  IL(0:6) = (/'s','p','d','f','g','h','i'/)
  real(REAL64), parameter    ::  EKBMIN = 1.0E-10_REAL64                 !  smaller KB energies (Hartree) are not written

! counters

  integer       ::  i, k, l


  if(iotape == 0) return

  if(nicore == 'nc  ') then
    lcore = .FALSE.
  elseif(nicore == 'pcec' .or. nicore == 'fcec') then
    lcore = .TRUE.
  else
    write(6,*) '  Unrecognized partial core correction, not writing psml pseudo'

    return

  endif

  if(icorr /= 'ca' .and. icorr /= 'pw' .and. icorr /= 'pb') then
    write(6,*) '  Unrecognized correlation, not writing psml pseudo'

    return

  endif

  if(irel == 'rel') then
    lso = .TRUE.
  elseif(irel == 'nrl') then
    lso = .FALSE.
  else
    write(6,*) '  Spin-polarized pseudopotential, not writing psml pseudo'

    return

  endif

  call atom_p_tbl_charge(nameat,iz)

! angular momenta present and their principal quantum numbers

  allocate(nq(0:mxdl), lhas(0:mxdl))
  allocate(f(nr))

  lhas(:) = .FALSE.
  lmax = 0
  do i = 1,npot(0)
    lhas(lo(i,0)) = .TRUE.
    lmax = max(lmax,lo(i,0))
  enddo

  do l = 0,mxdl
    nq(l) = l+1
  enddo
  do k = 1,10
    if(psdtitle(2*k-1)(3:3) == '(') then
      do l = 0,min(mxdl,6)
        if(psdtitle(2*k-1)(2:2) == IL(l)) then
          read(psdtitle(2*k-1)(1:1),'(i1)',iostat=ierr) n
          if(ierr == 0) nq(l) = n
        endif
      enddo
    endif
  enddo

  call date_and_time(date, time)
  call atom_kb_psml_uuid(uuid)

  open(unit = iotape, file = trim(fname), form = 'formatted', status = 'unknown')

  write(iotape,'(a)') '<?xml version="1.0" encoding="UTF-8" ?>'
  write(iotape,'(a)') '<psml version="1.1" energy_unit="hartree" length_unit="bohr"'
  write(iotape,'(a)') ' uuid="'//uuid//'"'
  write(iotape,'(a)') ' xmlns="http://esl.cecam.org/PSML/ns/1.1">'

! provenance

  write(iotape,'(a)') '<provenance creator="ATOM (J.L. Martins) '//trim(adjustl(irvers))//'"'
  write(iotape,'(a)') ' date="'//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//'">'
  write(iotape,'(a)') '<annotation action="semilocal-pseudopotential-generation"'
  write(iotape,'(a)') ' semilocal-date="'//trim(adjustl(irdate))//'"'
  write(iotape,'(a)') ' type="'//trim(xmlesc(trim(adjustl(irayps(1)//irayps(2)//irayps(3)))))//'"'
  write(iotape,'(a)') ' parameters="'//trim(xmlesc(trim(join(psdtitle))))//'"'
  write(iotape,'(a)') ' action-cont="kleinman-bylander-projectors" />'
  write(iotape,'(a)') '</provenance>'

! pseudo-atom-spec

  write(iotape,'(a)') '<pseudo-atom-spec atomic-label="'//trim(adjustl(nameat))//'"'
  write(iotape,'(a,i0,a,a,a)') ' atomic-number="', iz, '" z-pseudo="', trim(rfmt(zion)), '"'
  write(iotape,'(a)') ' flavor="Troullier-Martins (ATOM)"'
  if(lso) then
    write(iotape,'(a)') ' relativity="dirac" spin-dft="no"'
  else
    write(iotape,'(a)') ' relativity="no" spin-dft="no"'
  endif
  if(lcore) then
    write(iotape,'(a)') ' core-corrections="yes">'
  else
    write(iotape,'(a)') ' core-corrections="no">'
  endif

  write(iotape,'(a)') '<exchange-correlation>'
  write(iotape,'(a)') '<annotation atom-xc-code="'//icorr//'" />'
  write(iotape,'(a)') '<libxc-info number-of-functionals="2">'
  if(icorr == 'ca') then
    write(iotape,'(a)') '<functional name="Slater exchange" type="exchange" id="1" />'
    write(iotape,'(a)') '<functional name="Perdew &amp; Zunger" type="correlation" id="9" />'
  elseif(icorr == 'pw') then
    write(iotape,'(a)') '<functional name="Slater exchange" type="exchange" id="1" />'
    write(iotape,'(a)') '<functional name="Perdew &amp; Wang" type="correlation" id="12" />'
  else
    write(iotape,'(a)') '<functional name="Perdew, Burke &amp; Ernzerhof" type="exchange" id="101" />'
    write(iotape,'(a)') '<functional name="Perdew, Burke &amp; Ernzerhof" type="correlation" id="130" />'
  endif
  write(iotape,'(a)') '</libxc-info>'
  write(iotape,'(a)') '</exchange-correlation>'

  zval = ZERO
  nshell = 0
  do l = 0,lmax
    if(lhas(l) .and. zo(l) > ZERO) then
      zval = zval + zo(l)
      nshell = nshell + 1
    endif
  enddo
  write(iotape,'(a,a,a)') '<valence-configuration total-valence-charge="', trim(rfmt(zval)), '">'
  do l = 0,lmax
    if(lhas(l) .and. zo(l) > ZERO) then
      write(iotape,'(a,i0,a,a,a,a,a)') '<shell n="', nq(l), '" l="', IL(l),         &
          '" occupation="', trim(rfmt(zo(l))), '" />'
    endif
  enddo
  write(iotape,'(a)') '</valence-configuration>'
  write(iotape,'(a)') '</pseudo-atom-spec>'

! grid

  write(iotape,'(a,i0,a)') '<grid npts="', nr, '">'
  write(iotape,'(a)') '<annotation type="log" recipe="r(i) = a*(exp(b*(i-1))-1)" />'
  write(iotape,'(a)') '<grid-data>'
  write(iotape,'(1p,4e21.13)') (r(i),i=1,nr)
  write(iotape,'(a)') '</grid-data>'
  write(iotape,'(a)') '</grid>'

! valence charge q(r) = cdv / r^2

  do i = 2,nr
    f(i) = cdv(i) / (r(i)*r(i))
  enddo
  call extrap0(f)
  qval = integ(cdv)

  write(iotape,'(a,a,a)') '<valence-charge total-charge="', trim(rfmt(qval)), '"'
  write(iotape,'(a)') ' is-unscreening-charge="yes" rescaled-to-z-pseudo="no">'
  call write_radfunc(f, .TRUE.)
  write(iotape,'(a)') '</valence-charge>'

! pseudo-core charge

  if(lcore) then
    do i = 2,nr
      f(i) = cdc(i) / (r(i)*r(i))
    enddo
    call extrap0(f)
    write(iotape,'(a)') '<pseudocore-charge>'
    call write_radfunc(f, .TRUE.)
    write(iotape,'(a)') '</pseudocore-charge>'
  endif

! semilocal potentials V = vionic / (2 r) in Hartree

  if(lso) then
    call write_slps('scalar_relativistic', 0)
    call write_slps('spin_orbit', 2)
    call write_slps('lj', 1)
  else
    call write_slps('non_relativistic', 0)
  endif

! local potential

  do i = 2,nr
    f(i) = vlocal(i) / (2*r(i))
  enddo
  call extrap0(f)
  if(llocal < 0) then
    ltype = 'maximum-of-semilocal'
  else
    write(ltype,'("l=",i0)') llocal
  endif
  write(iotape,'(a)') '<local-potential type="'//trim(ltype)//'">'
  call write_radfunc(f, .FALSE.)
  write(iotape,'(a)') '</local-potential>'

! Kleinman-Bylander projectors

  if(lso) then
    call write_proj('scalar_relativistic', 0)
    call write_proj('spin_orbit', 2)
    call write_proj('lj', 1)
  else
    call write_proj('non_relativistic', 0)
  endif

! pseudo wave-functions R = rpsi / r

  if(lso) then
    call write_pswf('scalar_relativistic', 0)
    call write_pswf('lj', 1)
  else
    call write_pswf('non_relativistic', 0)
  endif

  write(iotape,'(a)') '</psml>'

  close(unit = iotape)

  deallocate(nq, lhas)
  deallocate(f)

  return

  contains

! semilocal potentials of a set.  mode 0: average, 1: j components, 2: spin-orbit

  subroutine write_slps(cname, mode)

    character(len=*), intent(in)    ::  cname
    integer, intent(in)             ::  mode

    integer      ::  jj, js

    cset = cname
    write(iotape,'(a)') '<semilocal-potentials set="'//trim(cset)//'">'
    do l = 0,lmax
      if(.not. lhas(l)) cycle
      if(mode == 2 .and. l == 0) cycle
      do jj = 1,2
        if(mode /= 1 .and. jj == 2) cycle
        if(mode == 1 .and. l == 0 .and. jj == 2) cycle
        js = 3 - 2*jj
        do i = 2,nr
          if(mode == 0) then
            f(i) = vionic(i,l,0) / (2*r(i))
          elseif(mode == 1) then
            f(i) = vionic(i,l,js) / (2*r(i))
          else
            f(i) = 2*(vionic(i,l,1) - vionic(i,l,-1)) / ((2*l+1)*2*r(i))
          endif
        enddo
        call extrap0(f)
        if(mode == 1) then
          xj = l + js*UM/2
          write(iotape,'(a,i0,a,a,a,a,a,a,a)') '<slps n="', nq(l), '" l="', IL(l),        &
              '" j="', trim(rfmt(xj)), '" rc="', trim(rfmt(rc(l))), '">'
        else
          write(iotape,'(a,i0,a,a,a,a,a)') '<slps n="', nq(l), '" l="', IL(l),            &
              '" rc="', trim(rfmt(rc(l))), '">'
        endif
        call write_radfunc(f, .FALSE.)
        write(iotape,'(a)') '</slps>'
      enddo
    enddo
    write(iotape,'(a)') '</semilocal-potentials>'

  end subroutine write_slps

! projectors of a set.  mode 0: scalar (j components weighted if lso), 1: lj, 2: spin-orbit

  subroutine write_proj(cname, mode)

    character(len=*), intent(in)    ::  cname
    integer, intent(in)             ::  mode

    integer      ::  jj, js

    cset = cname
    write(iotape,'(a)') '<nonlocal-projectors set="'//trim(cset)//'">'
    do l = 0,lmax
      nseq = 0
      do jj = 1,2
        if(lso) then
          js = 3 - 2*jj
          if(l == 0 .and. js == -1) cycle
        else
          if(jj == 2) cycle
          if(l == llocal) cycle
          js = 0
        endif
        if(inorm(l,js) == 0) cycle

!       V_KB = inorm |vkbproj><vkbproj| in Rydberg, vkbproj = chi (no factor r)

        do i = 1,nr
          f(i) = r(i)*r(i)*vkbproj(i,l,js)*vkbproj(i,l,js)
        enddo
        pnorm = integ(f)
        ekb = (inorm(l,js)*pnorm)/2

!       a vanishing projector (the channel of the local potential) is
!       numerical noise and is not written

        if(abs(ekb) < EKBMIN) cycle

        if(.not. lso .or. l == 0) then
          fsr = UM
          fso = ZERO
        elseif(js == 1) then
          fsr = (l+1)/(2*l+UM)
          fso = 2/(2*l+UM)
        else
          fsr = l/(2*l+UM)
          fso = -2/(2*l+UM)
        endif

        if(mode == 2 .and. fso == ZERO) cycle

        do i = 2,nr
          f(i) = vkbproj(i,l,js) / sqrt(pnorm)
        enddo
        if(l == 0) then
          call extrap0(f)
        else
          f(1) = ZERO
        endif

!       in the lj set seq counts the projectors of each l, j

        nseq = nseq + 1
        if(mode == 1) nseq = 1
        if(mode == 1) then
          xj = l + js*UM/2
          write(iotape,'(a,a,a,a,a,i0,a,a,a)') '<proj l="', IL(l), '" j="', trim(rfmt(xj)),  &
              '" seq="', nseq, '" ekb="', trim(rfmt(ekb)), '" type="KB">'
        elseif(mode == 0) then
          write(iotape,'(a,a,a,i0,a,a,a)') '<proj l="', IL(l), '" seq="', nseq,            &
              '" ekb="', trim(rfmt(fsr*ekb)), '" type="KB">'
        else
          write(iotape,'(a,a,a,i0,a,a,a)') '<proj l="', IL(l), '" seq="', nseq,            &
              '" ekb="', trim(rfmt(fso*ekb)), '" type="KB">'
        endif
        call write_radfunc(f, .TRUE.)
        write(iotape,'(a)') '</proj>'
      enddo
    enddo
    write(iotape,'(a)') '</nonlocal-projectors>'

  end subroutine write_proj

! pseudo wave-functions of a set.  mode 0: average, 1: j components

  subroutine write_pswf(cname, mode)

    character(len=*), intent(in)    ::  cname
    integer, intent(in)             ::  mode

    integer      ::  jj, js

    cset = cname
    write(iotape,'(a)') '<pseudo-wave-functions set="'//trim(cset)//'">'
    do l = 0,lmax
      if(.not. lhas(l)) cycle
      do jj = 1,2
        if(mode == 0 .and. jj == 2) cycle
        if(mode == 1 .and. l == 0 .and. jj == 2) cycle
        js = 0
        if(mode == 1) js = 3 - 2*jj
        do i = 2,nr
          f(i) = rpsi(i,l,js) / r(i)
        enddo
        if(l == 0) then
          call extrap0(f)
        else
          f(1) = ZERO
        endif
        if(mode == 1) then
          xj = l + js*UM/2
          write(iotape,'(a,i0,a,a,a,a,a,a,a)') '<pswf n="', nq(l), '" l="', IL(l),        &
              '" j="', trim(rfmt(xj)), '" energy_level="', trim(rfmt(ev(l,js)/2)), '">'
        else
          write(iotape,'(a,i0,a,a,a,a,a)') '<pswf n="', nq(l), '" l="', IL(l),            &
              '" energy_level="', trim(rfmt(ev(l,js)/2)), '">'
        endif
        call write_radfunc(f, .TRUE.)
        write(iotape,'(a)') '</pswf>'
      enddo
    enddo
    write(iotape,'(a)') '</pseudo-wave-functions>'

  end subroutine write_pswf

! writes a radfunc.  If ltrim the trailing negligible values are not written.

  subroutine write_radfunc(g, ltrim)

    real(REAL64), intent(in)        ::  g(:)
    logical, intent(in)             ::  ltrim

    integer      ::  np, ii
    real(REAL64) ::  gmax

    np = nr
    if(ltrim) then
      gmax = maxval(abs(g(1:nr)))
      do while(np > 2 .and. abs(g(np)) <= 1.0E-14_REAL64*gmax)
        np = np - 1
      enddo
    endif

    write(iotape,'(a)') '<radfunc>'
    write(iotape,'(a,i0,a)') '<data npts="', np, '">'
    write(iotape,'(1p,4e21.13)') (clean(g(ii)),ii=1,np)
    write(iotape,'(a)') '</data>'
    write(iotape,'(a)') '</radfunc>'

  end subroutine write_radfunc

! linear extrapolation to r(1) = 0

  subroutine extrap0(g)

    real(REAL64), intent(inout)     ::  g(:)

    g(1) = g(2) - (g(3)-g(2))*(r(2)-r(1))/(r(3)-r(2))

  end subroutine extrap0

! integral on the grid (Simpson in the index variable)

  function integ(g) result(s)

    real(REAL64), intent(in)        ::  g(:)
    real(REAL64)                    ::  s

    integer      ::  m, ii

    m = nr
    if(mod(m,2) == 0) m = m - 1
    s = g(1)*drdi(1) + g(m)*drdi(m)
    do ii = 2,m-1,2
      s = s + 4*g(ii)*drdi(ii)
    enddo
    do ii = 3,m-2,2
      s = s + 2*g(ii)*drdi(ii)
    enddo
    s = s / 3
    if(m < nr) s = s + (g(nr)*drdi(nr) + g(m)*drdi(m))/2

  end function integ

! values too small for a 3 digit exponent are written as zero

  function clean(x) result(y)

    real(REAL64), intent(in)        ::  x
    real(REAL64)                    ::  y

    y = x
    if(abs(x) < 1.0E-99_REAL64) y = ZERO

  end function clean

! real number for an attribute

  function rfmt(x) result(c)

    real(REAL64), intent(in)        ::  x
    character(len=24)               ::  c

    write(c,'(es22.14)') x
    c = adjustl(c)

  end function rfmt

! concatenation of psdtitle

  function join(t) result(c)

    character(len=10), intent(in)   ::  t(20)
    character(len=200)              ::  c

    integer      ::  m

    c = ' '
    do m = 1,20
      c(10*m-9:10*m) = t(m)
    enddo

  end function join

! escapes the XML special characters of an attribute value

  function xmlesc(s) result(c)

    character(len=*), intent(in)    ::  s
    character(len=6*len(s))         ::  c

    integer      ::  m, p

    c = ' '
    p = 0
    do m = 1,len_trim(s)
      select case(s(m:m))
      case('&')
        c(p+1:p+5) = '&amp;'
        p = p + 5
      case('<')
        c(p+1:p+4) = '&lt;'
        p = p + 4
      case('>')
        c(p+1:p+4) = '&gt;'
        p = p + 4
      case('"')
        c(p+1:p+6) = '&quot;'
        p = p + 6
      case default
        c(p+1:p+1) = s(m:m)
        p = p + 1
      end select
    enddo

  end function xmlesc

end subroutine atom_kb_psd_out_psml


!>  random (version 4) uuid for the psml file

subroutine atom_kb_psml_uuid(uuid)

! written 28 October 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  character(len=36), intent(out)    ::  uuid                             !<  uuid

  character(len=16), parameter  ::  HEX = '0123456789abcdef'

  integer                 ::  nseed, ival(8)
  integer, allocatable    ::  iseed(:)
  real(REAL64)            ::  x
  integer                 ::  k, m, d

  call random_seed(size = nseed)
  allocate(iseed(nseed))
  call date_and_time(values = ival)
  do k = 1,nseed
    iseed(k) = 1000*ival(8) + 60000*ival(7) + 3600000*ival(6) + 7919*k + 104729*ival(3)
  enddo
  call random_seed(put = iseed)
  deallocate(iseed)

  uuid = 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'
  do m = 1,36
    if(uuid(m:m) == 'x' .or. uuid(m:m) == 'y') then
      call random_number(x)
      d = int(16*x)
      if(uuid(m:m) == 'y') d = 8 + mod(d,4)
      uuid(m:m) = HEX(d+1:d+1)
    endif
  enddo

  return

end subroutine atom_kb_psml_uuid
