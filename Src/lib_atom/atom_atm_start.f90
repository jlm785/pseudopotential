!>  Reads calculation parameters, by default from atom.dat,
!>  and sets the initial values of several quantities.
!>
!>  \author       Sverre Froyen, Norm Troullier, Carlos Balbas, Jose Luis Martins
!>  \version      6.0.3
!>  \date         1980s, 4 August 2021, 12 September 2021.
!>  \copyright    GNU Public License v2

subroutine atom_atm_start(nameat, ctype, ititle, kerker, icorr, ispp,    &
    znuc, rmax, aa, bb, ncore, nval, ni, li, zd, zu, evd,                &
    itype, ikerk, nr, a, b, r, drdi, d2rodr,                             &
    norb, no, lo, iso, zo, zel, evi,                                     &
    iowrite, mxdnr, mxdorb, mxdl)


! njtj ***  modifications  ***
!   The input and output variables passed have been changed.
!   There are five new pseudopotential generation options
!   The input variables znuc,zsh,rsh,rmax,aa,bb are
!   compared to a small positive value - eliminates
!   floating point comparisions errors(zero is
!   not always zero).
! njtj ***  modifications  ***

! lcb
!   modified for GGA by Carlos Balbas,   January 97.
! lcb

! cleanup and new interface, July 2019. JLM
! minor cleanup, new name, 22 June 2021. JLM
! separate read file from processing, 4 August 2021. JLM
! so->iso 12 September 2021. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

!  data for orbitals, 1s,2s,2p,3s,3p,3d,4s,4p,4d,5s,5p,4f,5d,6s,6p,5f,6d,7s,7p

  integer, parameter    ::  ncmax = 19                                   !  maximum number of core orbitals
  integer, parameter    ::  nc(ncmax) = (/1,2,2,3,3,3,4,4,4,5,5,4,5,6,6,5,6,7,7/)
  integer, parameter    ::  lc(ncmax) = (/0,0,1,0,1,2,0,1,2,0,1,3,2,0,1,3,2,0,1/)


! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdorb                           !<  dimension of the number of orbitals
  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points
  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum l

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbol of the atom

  character(len=2), intent(in)      ::  ctype                            !<  type of calculation flag
  character(len=10), intent(in)     ::  ititle(5)                        !<  title of calculation

  character(len=3), intent(in)      ::  kerker                           !<  type of pseudopotential flag

  character(len=2), intent(in)      ::  icorr                            !<  correlation type

  integer, intent(in)               ::  ni(mxdorb)                       !<  valence principal quantum number
  integer, intent(in)               ::  li(mxdorb)                       !<  valence angular quantum number
  real(REAL64), intent(in)          ::  zd(mxdorb), zu(mxdorb)           !<  occupation of valence down and up orbitals
  real(REAL64), intent(in)          ::  evd(mxdorb)                      !<  default eigenvalue



! output

  integer, intent(out)              ::  norb                             !<  dimension of the number of orbitals

  integer, intent(out)              ::  no(mxdorb)                       !<  principal quantum number n
  integer, intent(out)              ::  lo(mxdorb)                       !<  angular quantum number l
  integer, intent(out)              ::  iso(mxdorb)                      !<  2*spin or 2*(j-l)
  real(REAL64), intent(out)         ::  zo(mxdorb)                       !<  orbital occupation
  real(REAL64), intent(out)         ::  evi(mxdorb)                      !<  orbital default energy

  integer, intent(out)              ::  itype                            !<  type of calculation (-1 signals end of calculation)
  integer, intent(out)              ::  ikerk                            !<  type of pseudopotential

  real(REAL64), intent(out)         ::  zel                              !<  electron charge

! output and input

  real(REAL64), intent(inout)       ::  znuc                             !<  nuclear charge

  integer, intent(inout)            ::  ncore                            !<  number of orbitals treated as core
  integer, intent(inout)            ::  nval                             !<  number of valence orbitals

  character(len=1), intent(inout)   ::  ispp                             !<  spin polarization  ' ', 's', 'r'

  real(REAL64), intent(inout)       ::  rmax                             !<  maximum mesh radius
  real(REAL64), intent(inout)       ::  aa                               !<  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64), intent(inout)       ::  bb                               !<  a = exp(-aa)/znuc, b = 1/bb

  real(REAL64), intent(inout)       ::  a                                !<  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64), intent(inout)       ::  b                                !<  r(i) = a*(exp(b*(i-1))-1)

  integer, intent(inout)            ::  nr                               !<  dimension of the number of radial points

  real(REAL64), intent(inout)       ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(inout)       ::  drdi(mxdnr)                      !<  d r(i) / d i
  real(REAL64), intent(inout)       ::  d2rodr(mxdnr)                    !<  (d^2 r(i) / d i^2) / (d r / d i)

! local

  integer, allocatable           ::  nomin(:)

  real(REAL64)      ::  zcore                                            !  core charge
  real(REAL64)      ::  zval                                             !  valence charge
  real(REAL64)      ::  zion                                             !  ionic charge
  integer           ::  isc, isi                                         !  2*spin or 2*(j-l) in core or valence

  character(len=10) ::  iray(5)                                          !  for output
  character(len=3)  ::  name                                             !  prefix to relativistic
  real(REAL64)      ::  xji                                              !  0 or l+s

  integer           ::  iz
  character(len=9)  ::  bdate

  integer           ::  jmax

! counters

  integer     ::  i, j, l

! constants

  real(REAL64), parameter    ::  UM = 1.0_REAL64, ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  PFIVE = 0.5_REAL64
  real(REAL64), parameter    ::  EPS = 0.00001_REAL64


  allocate(nomin(0:mxdl))

  do l = 0,mxdl
    nomin(l) = 1000
  enddo

  do i = 1,mxdorb
    no(i) = 0
    lo(i) = 0
    iso(i) = 0
    zo(i) = ZERO
    evi(i) = ZERO
  enddo

!   converts the type of calculation
!   ctypr -> itype
!   ae -> 0 all electron calculation
!   pg -> 1 pseudopotential generation w/o core correction
!   pe -> 2 pseudopotential generation w/  core correction exchange
!   ph -> 3 pseudopotential generation w/  core correction hartree/exc
!   pt -> 4 pseudopotential test
!   pm -> 5 pseudopotential test + valence charge modify


!  if ctype = '  ' , no more data, program ends

  if (ctype == 'ae') then
    itype = 0
  elseif (ctype == 'pg') then
    itype = 1
  elseif (ctype == 'pe') then
    itype = 2
  elseif (ctype == 'ph') then
    itype = 3
  elseif (ctype == 'pt') then
    itype = 4
  elseif (ctype == 'pm') then
    itype = 5
  else
    itype =-1

    return

  endif

! njtj  ***  major modification  start  ***
! There are seven ways to generate the pseudopotential :
!    kerker = van Vanderbilt
!    kerker = tam Troullier and Martins
!    kerker = ker (yes) Kerker
!    kerker = hsc (no)  Hamann Schluter and Chiang
!    kerker = min (oth) datafile made for minimization
!    kerker = bhs Bachelet, Hamann and Schluter
!    kerker = tm2 Improved Troullier and Martins

  if (itype > 0) then

    if(kerker == 'tm2' .or. kerker == 'TM2') then
      ikerk = 6
    elseif(kerker == 'bhs' .or. kerker == 'BHS') then
       ikerk = 5
    elseif(kerker == 'oth' .or. kerker == 'OTH' .or.                     &
           kerker == 'min' .or. kerker == 'MIN') then
      ikerk = 4
    elseif (kerker == 'van' .or. kerker =='VAN') then
      ikerk = 3
    elseif (kerker == 'tbk' .or. kerker == 'TBK' .or.                     &
            kerker == 'tam' .or. kerker == 'TAM') then
      ikerk = 2
    elseif (kerker == 'yes' .or. kerker == 'YES' .or.                     &
            kerker == 'ker' .or. kerker == 'KER') then
      ikerk = 1
    elseif (kerker == 'no ' .or. kerker == ' no' .or.                     &
            kerker == 'NO ' .or. kerker == ' NO' .or.                     &
            kerker== 'hsc' .or. kerker == 'HSC') then
      ikerk = 0
    else
      write(iowrite,'(//,"error in input_start - kerker = ",a3,          &
        &      " unknown")')  kerker

      STOP                                                 !  exits program with error

    endif
  endif

! njtj  ***  major modification end  ***

!   processes spin/relativistic
!   ispp = ' ' - nonspin calculation
!   ispp = s  - spin polarized calculation
!   ispp = r  - relativistic calculation

       if (ispp /= 's' .and. ispp /= 'r') ispp = ' '

! jlm    options that are no longer available
!    if (ispp == 's' .and. icorr == 'xa') ispp = ' '
!    if (ispp == 's' .and. icorr == 'wi') ispp = ' '
!    if (ispp == 's' .and. icorr == 'hl') ispp = ' '

! njtj   ***  major modification start  ***
!   Floating point comparison error modification.

  if (abs(znuc) < EPS) then
    call atom_p_tbl_charge(nameat,iz)
    znuc = UM*iz
  endif

  if (itype < 4) then

! set up grid if you are not testing or modifying a pseudopotential

    if (abs(rmax) < EPS) rmax = 120*UM
    if (abs(aa) < EPS) aa = 6*UM
    if (abs(bb) < EPS) bb = 200*UM

!   bb = 40.0 standard
!   bb = 80.0 2 * standard
!   bb = 120.0 3 * standard

    call atom_atm_setmesh(iowrite, aa, bb, znuc, rmax, nr, r ,drdi, d2rodr, mxdnr)

    a = exp(-aa)/znuc
    b = UM/bb

  endif

! njtj  ***  major modification end  ***

! checks number of core and valence orbitals


  if (ncore > ncmax .or. ncore < 0) then
    write(iowrite,'(//," error in input_start - max number of core",     &
       &    " orbitals is ",i5, "input value is ",i5)')  ncmax, ncore

    stop                                                   !  exits program with error

  endif
  if (nval < 0) then
    write(iowrite,'(//," error in input_start -number of",               &
      &  " valence orbitals is ",i5)')  nval

    stop                                                   !  exits program with error

  endif

! compute occupation numbers and orbital quantum numbers for the core

  zcore = ZERO
  if (ncore /= 0 .and. itype < 4) then

    norb = 0
    do i = 1,ncore

      if (ispp == ' ') then
        jmax = 1
        isc = 0
      elseif (ispp == 'r' .and. lc(i) == 0) then
        jmax = 1
        isc = 1
      else
        jmax = 2
        isc =-1
      endif

      do j = 1,jmax

        norb = norb + 1
        no(norb) = nc(i)
        lo(norb) = lc(i)
        iso(norb) = isc

        zo(norb) = (2*lo(norb)+1)*UM
        if (ispp == ' ') zo(norb) = 2*zo(norb)
        if (ispp == 'r') zo(norb) = (2*(lo(norb)+isc*PFIVE)+1)*UM

        zcore = zcore + zo(norb)
        if (ispp /= ' ') isc = -isc
      enddo

    enddo
    ncore = norb
  endif

! compute occupation numbers and orbital quantum numbers for the valence

  if (itype >= 4) ncore = 0
  norb = ncore
  zval = ZERO

  if (nval /= 0) then

    do i = 1,nval

      if (ispp == ' ') then
        jmax = 1
        isi = 0
      elseif (ispp == 'r' .and. li(i) == 0) then
        jmax = 1
        isi = 1
      else
        jmax = 2
        isi =-1
      endif

      do j = 1,jmax

        norb = norb + 1
        no(norb) = ni(i)
        lo(norb) = li(i)
        iso(norb) = isi
        zo(norb) = zd(i) + zu(i)
        if (zo(norb) == ZERO) evi(norb) = evd(i)

        if (ispp == 's') then
          if (isi == 0 .or. isi == 1) then
            zo(norb) = zd(i)
          else
            zo(norb) = zu(i)
          endif
        elseif (ispp == 'r') then
          zo(norb) = zo(norb)*(2*(li(i)+isi*PFIVE)+1)/(4*li(i)+2)
        endif

        zval = zval + zo(norb)
        if (norb /= 0 .and. nomin(lo(norb)) > no(norb)) then
          nomin(lo(norb)) = no(norb)
        endif
        if (ispp /= ' ') isi = -isi

      enddo

    enddo

! abort if two orbitals are equal

    nval = norb - ncore
    do i = 1,norb
      do j = 1,norb
        if (i <= j) exit
        if (no(i) /= no(j)) exit
        if (lo(i) /= lo(j)) exit
        if (iso(i) /= iso(j)) exit
        write(iowrite,'(//," error in input - orbital ",i2,              &
             &   " is already occupied")') i

        stop                                                        !  exits program with error

      enddo
    enddo

! reduce n quantum number if pseudoatom

    if (itype >= 4) then
      do i = 1,nval
        no(i) = no(i) - nomin(lo(i)) + lo(i) + 1
      enddo
    endif

  endif

! end of loop over valence orbitals

  zion = znuc - zcore - zval
  zel = zval
  if (itype < 4) then
    zel = zel + zcore
  else
    znuc = znuc - zcore
  endif

! find jobname and date and printout, zedate is a machine dependent routine

  if (icorr == 'pb') then
    iray(1) = 'atom-GGA96'
  elseif (icorr== 'pw') then
    iray(1) = 'atom-LDApw'
  elseif (icorr == 'ca') then
    iray(1)='atom-LDAca'
  else
    stop 'unrecognized correlation choice'
  endif

  call zedate(bdate)
  iray(2) = bdate//' '

! printout

  write(iowrite,'(1x,a10,a10,5x,5a10,/,21("*"),/)') iray(1),iray(2),ititle

  if (itype == 0) then
    write(iowrite,'(1x,a2," all electron calculation ",/,1x,27("-"),/)') nameat
  elseif (itype < 4) then
    write(iowrite,'(1x,a2," pseudopotential generation",/,1x,29("-"),/)') nameat
  elseif (itype == 4) then
    write(iowrite,'(1x,a2," pseudopotential test",/,1x,23("-"),/)') nameat
  elseif (itype == 5) then
    write(iowrite,'(1x,a2," pseudo test + charge mod ",/,1x,27("-"),/)') nameat
  endif

  if (ispp == 'r') then
    write(iowrite,'(" R E L A T I V I S T I C ! !",/)')
    name = '   '
  elseif (ispp == ' ') then
    name = 'non'
  else
    name = '   '
  endif

  write(iowrite,'(" correlation = ",a2,3x,a3,"spin-polarized",/)') icorr,name

  write(iowrite,'(" nuclear charge             =",f10.6,/,               &
       & " number of core orbitals    =",i3,/,                           &
       & " number of valence orbitals =",i3,/,                           &
       & " electronic charge          =",f10.6,/,                        &
       & " ionic charge               =",f10.6,//)')                     &
         znuc,ncore,nval,zel,zion

!   if (zsh > EPS) write(iowrite,'(" shell charge =",f6.2," at radius =",  &
!        &  f6.2,//)') zsh,rsh

  write(iowrite,'(" input data for orbitals",//,                         &
       & "  i    n    l    s     j     occ",/)')
  xji = ZERO
  do i=1,norb
    if (ispp == 'r') xji = lo(i) + iso(i)*PFIVE
    write(iowrite,'(1x,i2,2i5,2f6.1,f10.4)') i, no(i), lo(i),            &
        &      iso(i)*PFIVE, xji, zo(i)
  enddo

  if (itype < 4) then
    write(iowrite,'(//," radial grid parameters",//," r(1) = .0 , r(2) =",    &
         & e11.3," , ... , r(",i4,") =",f8.3,/," a =",f7.3,"  b =",f8.3,/)')  &
         r(2), nr, r(nr), aa, bb
  endif

  deallocate(nomin)

  return
end subroutine atom_atm_start
