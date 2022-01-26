!>  reads calculation parameters, by default from atom.dat
!>
!>  \author       Sverre Froyen, Norm Troullier, Carlos Balbas, Jose Luis Martins
!>  \version      6.013
!>  \date         1980s, 4 August 2021
!>  \copyright    GNU Public License v2

subroutine atom_kb_test_start(nameat, ititle, icorr, ispp,               &
    nval, ni, li, zd, zu, evd,                                           &
    norb, no, lo, iso, zo, zel, evi,                                     &
    iowrite, mxdorb, mxdl)


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


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdorb                           !<  dimension of orbitals
  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbol of the atom

  character(len=10), intent(in)     ::  ititle(5)                        !<  title of calculation

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

  real(REAL64), intent(out)         ::  zel                              !<  electron charge

! output and input

  integer, intent(inout)            ::  nval                             !<  number of valence orbitals

  character(len=1), intent(inout)   ::  ispp                             !<  spin polarization  ' ', 's', 'r'


! local

  integer           ::  nomin(0:mxdl)
  real(REAL64)      ::  zval                                             !  valence charge
!  real(REAL64)      ::  zion                                             !  ionic charge
  real(REAL64)      ::  znuc                                             !  nuclear charge
  integer           ::  isi                                              !  spin in core or valence
  character(len=10) ::  iray(5)                                          !  for output
  character(len=3)  ::  name                                             !  prefix to relativistic
  real(REAL64)      ::  xji                                              !  0 or l+s

  integer           ::  iz
  character(len=9)  ::  bdate

  integer           ::  jmax

! counters

  integer     ::  i, j

! constants

  real(REAL64), parameter    ::  UM = 1.0_REAL64, ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  PFIVE = 0.5_REAL64
  real(REAL64), parameter    ::  EPS = 0.00001_REAL64

  do i = 0,mxdl
    nomin(i) = 100
  enddo

  do i = 1,mxdorb
    no(i) = 0
    lo(i) = 0
    iso(i) = 0
    zo(i) = ZERO
    evi(i) = ZERO
  enddo

  call atom_p_tbl_charge(nameat,iz)
  znuc = UM*iz

! compute occupation numbers and orbital quantum numbers for the valence


  norb = 0
  zval = ZERO

  if (nval > 0) then

    do i = 1,nval

      if (ispp == ' ') then
        jmax = 1
        isi = ZERO
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

    nval = norb
    do i = 1,norb
      do j = 1,norb
        if (i <= j) exit
        if (no(i) /= no(j)) exit
        if (lo(i) /= lo(j)) exit
        if (iso(i) /= iso(j)) exit
        write(iowrite,'(//," error in input - orbital ",i2,              &
             & " is already occupied")') i

        stop                                                        !  exits program with error

      enddo
    enddo

! reduce n quantum number if pseudoatom

      do i = 1,nval
        no(i) = no(i) - nomin(lo(i)) + lo(i) + 1
      enddo

  endif

! end of loop over valence orbitals

  zel = zval

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

  write(iowrite,'(1x,a2," KB pseudopotential test",/,1x,23("-"),/)') nameat

  if (ispp == ' ') then
    name = 'non'
  else
    name = '   '
  endif

  write(iowrite,'(" correlation = ",a2,3x,a3,"spin-polarized",/)') icorr,name

  write(iowrite,'(" nuclear charge             =",f10.6,/,               &
       & " number of valence orbitals =",i3,/,                           &
       & " electronic charge          =",f10.6,/)')                      &
         znuc,nval,zel

!
  write(iowrite,'(" input data for orbitals",//,                         &
       & "  i    n    l    s     occ",/)')
  xji = ZERO
  do i=1,norb
    write(iowrite,'(1x,i2,2i5,f6.1,f10.4)') i, no(i), lo(i), iso(i)*PFIVE, zo(i)
  enddo

  return
end subroutine atom_kb_test_start
