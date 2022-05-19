!>  reads by default from atom.dat and obtain the required array dimensions
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.8
!>  \date         22 June 2021, 21 April 2022.
!>  \copyright    GNU Public License v2

subroutine atom_atm_input_size(ioread, iopsd, filepsd, mxdnr, mxdorb, mxdl)

! calls default mesh. 21 April 2022. JLM
! psdtitle, 18 May 2022. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

!  data for orbitals, 1s,2s,2p,3s,3p,3d,4s,4p,4d,5s,5p,4f,5d,6s,6p,5f,6d,7s,7p

  integer, parameter    ::  ncmax = 19                                   !  maximum number of core orbitals
  integer, parameter    ::  nc(ncmax) = (/1,2,2,3,3,3,4,4,4,5,5,4,5,6,6,5,6,7,7/)
  integer, parameter    ::  lc(ncmax) = (/0,0,1,0,1,2,0,1,2,0,1,3,2,0,1,3,2,0,1/)


! input

  integer, intent(in)               ::  ioread                           !<  default tape for reading

  integer, intent(in)               ::  iopsd                            !<  default tape for pseudopotential
  character(len=*), intent(in)      ::  filepsd                          !<  name of default tape for reading pseudopotential in parsec format

! output

  integer, intent(out)              ::  mxdorb                           !<  dimension of the number of orbitals
  integer, intent(out)              ::  mxdnr                            !<  dimension of the number of radial points
  integer, intent(out)              ::  mxdl                             !<  dimension for angular momentum l

! local variables

  integer                            ::  norb                            !  number of orbitals


  integer                           ::  itype                            !  type of calculation (-1 signals end of calculation)
  character(len=2)                  ::  nameat                           !  chemical symbol of the atom
  character(len=2)                  ::  icorr                            !  correlation type
  character(len=1)                  ::  ispp                             !  spin polarization  ' ', 's', 'r'

  real(REAL64)                      ::  znuc                             !  nuclear charge

  real(REAL64)                      ::  zsh                              !  shell charge (seldom used)
  real(REAL64)                      ::  rsh                              !  shell radius (seldom used)

  integer                           ::  ncore                            !  number of orbitals treated as core

  real(REAL64)                      ::  a                                !  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64)                      ::  b                                !  r(i) = a*(exp(b*(i-1))-1)

  real(REAL64)                      ::  ri                               !  r(i) = a*(exp(b*(i-1))-1)

! local

  character(len=2)           ::  ctype                                   !  type of calculation flag
  character(len=3)           ::  kerker                                  !  type of pseudopotential flag
  real(REAL64)               ::  rmax                                    !  maximum mesh radius
  real(REAL64)               ::  aa                                      !  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64)               ::  bb                                      !  a = exp(-aa)/znuc, b = 1/bb
  integer                    ::  nval                                    !  number of valence orbitals
  integer                    ::  ni, li                                  !  principal and angular quantum number
  integer                    ::  iz                                      !  atomic number
  integer                    ::  lmax                                    !  maximum angular momentum l

  integer                    ::  imax

  character(len=2)           ::  icorrt, namet
  character(len=3)           ::  irel
  character(len=4)           ::  nicore
  character(len=10)          ::  iray(6), psdtitle(20)

  integer                    ::  nrm
  integer                    ::  npotd, npotu

  real(REAL64)               ::  rmax_def, aa_def, bb_def

! counters

  integer     ::  i

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  EPS = 0.00001_REAL64

  read(ioread,'(3x,a2)') ctype

!  if ctype = ' ' , no more data, program ends

  if (ctype == 'ae') then
    itype=0
  elseif (ctype == 'pg') then
    itype=1
  elseif (ctype == 'pe') then
    itype=2
  elseif (ctype == 'ph') then
    itype=3
  elseif (ctype == 'pt') then
    itype=4
  elseif (ctype == 'pm') then
    itype=5
  else
    itype=-1

    return

  endif

  if (itype > 0) then
    read(ioread,*) kerker
  endif

  read(ioread,'(3x,a2,3x,a2,a1)') nameat,icorr,ispp

  if (ispp /= 's' .and. ispp /= 'r') ispp = ' '

  read(ioread,*) znuc,zsh,rsh,rmax,aa,bb

  if (abs(znuc) < EPS) then
    call atom_p_tbl_charge(nameat,iz)
    znuc = ONE*iz
  endif

  if (itype < 4) then

! set up grid if you are not testing or modifying a pseudopotential

    call atom_atm_default_mesh(rmax_def, aa_def, bb_def)

    if (abs(rmax) < EPS) rmax = rmax_def
    if (abs(aa) < EPS) aa = aa_def
    if (abs(bb) < EPS) bb = bb_def

!   bb = 40.0 standard
!   bb = 80.0 2 * standard
!   bb = 120.0 3 * standard

!   call atom_atm_setmesh(iowrite, aa, bb, znuc, rmax, nr, r ,drdi, d2rodr, mxdnr)

    a = exp(-aa)/znuc
    b = ONE/bb

    imax = nint( log(rmax/a+ONE) / b + ONE )

    do i = min(2,imax-5),imax+5
      mxdnr = i-1
      ri = a*(exp(b*(i-1))-1)
      if (ri > rmax) exit
    enddo

  else

    open(unit=iopsd, file=trim(filepsd), form='unformatted', status='unknown')

    read(iopsd) namet, icorrt, irel ,nicore, (iray(i),i=1,6),            &
      (psdtitle(i),i=1,7), npotd, npotu, nrm

    mxdnr = nrm + 1

    close(unit=iopsd)

  endif

! njtj  ***  major modification end  ***

! read the number of core and valence orbitals


  read(ioread,*) ncore,nval

  if (ncore > ncmax .or. ncore < 0) then
    write(6,*)
    write(6,*)
    write(6,'(" error in atm_input_size:   max number of core",          &
         & " orbitals is ",i5, "input value is ",i5)')  ncmax, ncore

    STOP                                                   !  exits program with error

  endif

! compute occupation numbers and orbital quantum numbers for the core

  norb = 0
  lmax = 0
  if (ncore /= 0) then
    do i = 1,ncore

      if(ispp == 's') then
        norb = norb + 2
      elseif(ispp == 'r') then
        if(lc(i) == 0) then
          norb = norb + 1
        else
          norb = norb + 2
        endif
      else
        norb = norb + 1
      endif

      if (lmax < lc(i)) lmax = lc(i)

    enddo

  endif

! compute occupation numbers and orbital quantum numbers for the valence

  if (itype >= 4) then
    lmax = 0
    norb = 0
  endif

  if (nval /= 0) then

    do i = 1,nval

      read(ioread,'(2i5)') ni,li

      if(ispp == 's') then
        norb = norb + 2
      elseif(ispp == 'r') then
        if(li == 0) then
          norb = norb + 1
        else
          norb = norb + 2
        endif
      else
        norb = norb + 1
      endif

      if (lmax < li) lmax = li

    enddo

  endif

  mxdorb = norb

! allow excitation in tests

  mxdl = lmax+2

  return

end subroutine atom_atm_input_size
