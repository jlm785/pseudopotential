!>  Reads the results of the all-electron calculation
!>  preparation of the pseudopotential construction.
!>  Reverse of atom_atm_datout
!>
!>  \author       Sverre Froyen,  Norm Troullier, Jose Luis Martins
!>  \version      6.0.3
!>  \date         1980s and 1990s, April 2018, 23 September 2021
!>  \copyright    GNU Public License v2

subroutine atom_psd_dat_in(ioae, fileae, itype, icorr, ispp, lmax,       &
     & nr, a, b, r, drdi, d2rodr,                                        &
     & nameat, norb, ncore, no, lo, iso, zo, znuc, zel,                  &
     & cdc, vionic, vhxc, ev,                                            &
     & mxdnr, mxdorb, mxdl)

! converted to f90, April 2018
! vionic, so->iso, vhxc. 14, 23 September 2021. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  ioae                             !<  default tape for writing
  character(len=*), intent(in)      ::  fileae                           !<  name of default tape for all-electron results


  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points
  integer, intent(in)               ::  mxdorb                           !<  dimension of the number of orbitals
  integer, intent(in)               ::  mxdl                             !<  dimension maximum angular momentum

! output

  integer, intent(out)              ::  itype                            !<  type of calculation (-1 signals end of calculation)
  character(len=2), intent(out)     ::  icorr                            !<  correlation type
  character(len=1), intent(out)     ::  ispp                             !<  spin polarization  ' ', 's', 'r'

  real(REAL64), intent(out)         ::  a                                !<  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64), intent(out)         ::  b                                !<  r(i) = a*(exp(b*(i-1))-1)

  integer, intent(out)              ::  nr                               !<  number of radial points
  real(REAL64), intent(out)         ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(out)         ::  drdi(mxdnr)                      !<  d r(i) / d i
  real(REAL64), intent(out)         ::  d2rodr(mxdnr)                    !<  (d^2 r(i) / d i^2) /  (d r(i) / d i)

  character(len=2), intent(out)     ::  nameat                           !<  chemical symbot of the atom

  integer, intent(out)              ::  norb                             !<  number of orbitals
  integer, intent(out)              ::  ncore                            !<  number of orbitals treated as core

  integer, intent(out)              ::  lmax                             !<  maximum angular momentum+1

  integer, intent(out)              ::  no(mxdorb)                       !<  principal quantum number n
  integer, intent(out)              ::  lo(mxdorb)                       !<  angular quantum number l
  integer, intent(out)              ::  iso(mxdorb)                      !<  2*spin or 2*(j-l)
  real(REAL64), intent(out)         ::  zo(mxdorb)                       !<  orbital occupation

  real(REAL64), intent(out)         ::  znuc                             !<  nuclear charge
  real(REAL64), intent(out)         ::  zel                              !<  electron charge

  real(REAL64), intent(out)         ::  cdc(mxdnr)                       !<  core charge density (total)

  real(REAL64), intent(out)         ::  vionic(mxdnr,0:mxdl,-1:1)        !<  r*ionic potential in Rydberg (down or total)

  real(REAL64), intent(out)         ::  vhxc(mxdnr,-1:1)                 !<  screening potential (Ry)

  real(REAL64), intent(out)         ::  ev(mxdorb)                       !<  orbital energy

! local variables

  real(REAL64), allocatable         ::  so(:)

  real(REAL64)                      ::  rtry, drtry
  integer                           ::  lmaxp1
  integer                           ::  id

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  SMALL = 1.0E-8_REAL64

! counter

  integer       ::  i, l


  allocate(so(mxdorb))

  vionic = ZERO
  vhxc = ZERO

! Open and read out data (default: from  datafile.dat).

  open (unit=ioae, file=trim(fileae), status='old', form='unformatted')

  read(ioae) itype, icorr, ispp, nr, a, b

  if(nr > mxdnr) then
    write(6,*) '  Stopped in atom_psd_dat_in:'
    write(6,*) '  nr = ',nr,'  mxdnr = ',mxdnr

    STOP

  endif

  if(ispp == ' ') then
    id = 0
  else
    id = 1
  endif

  read(ioae) (r(i),i=1,nr)
  read(ioae) (drdi(i),i=1,nr)

! paranoid check

  do i = 1,nr
    rtry = a*(exp(b*(i-1))-1)

    if(abs(rtry - r(i)) > SMALL) then
      write(6,*) '  Stopped in atom_psd_dat_in:'
      write(6,'("   INCOMPATIBLE  a, b r(i) ")')

      STOP

    endif

    drtry = (r(i)+a)*b

    if(abs(drtry - drdi(i)) > SMALL) then
      write(6,*) '  Stopped in atom_psd_dat_in:'
      write(6,'("   INCOMPATIBLE  a, b drdi(i) ")')

      STOP

    endif

    d2rodr(i) = b
  enddo


  read(ioae) lmaxp1, nameat, norb, ncore
! old convention
  lmax = lmaxp1-1

  if(lmax > mxdl) then
    write(6,*) '  Stopped in atom_psd_dat_in:'
    write(6,*) '  lmax = ',lmax,'  mxdl = ',mxdl

    STOP

  endif

  if(norb > mxdorb) then
    write(6,*) '  Stopped in atom_psd_dat_in:'
    write(6,*) '  norb = ',norb,'  mxdorb = ',mxdorb

    STOP

  endif

  read(ioae) (no(i),i=1,norb)
  read(ioae) (lo(i),i=1,norb)
  read(ioae) (so(i),i=1,norb)

  do i = 1,norb
    if(so(i) > 0.1) then
      iso(i) = 1
    elseif(so(i) < -0.1) then
      iso(i) =-1
    else
      iso(i) = 0
    endif
  enddo

  read(ioae) (zo(i),i=1,norb)
  read(ioae) znuc, zel
  read(ioae) (cdc(i),i=1,nr)

  do l = 0,lmax
    read(ioae) (vionic(i,l,id),i=1,nr)
  enddo
  do l = 0,lmax
    read(ioae) (vionic(i,l,-1),i=1,nr)
  enddo

  read(ioae) (vhxc(i,id),i=1,nr)
  read(ioae) (vhxc(i,-1),i=1,nr)
  read(ioae) (ev(i),i=1,norb)

  close (unit=ioae)


  deallocate(so)

  return

end subroutine atom_psd_dat_in
