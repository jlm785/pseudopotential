!>  sets the quasi-logarithmic atomic mesh
!>
!>  reasonable values  aa ~ 1,...,6,  bb ~ 40,...,200, rmax ~ 100,...,200
!>  accuracy (mesh size) is mostly controled by bb.
!>  r(1) = 0, r(2) ~ exp(-aa) / bb,    r(j+1) ~ r(j)*(1+1/bb)
!>
!>  r(i) = a*(exp(b*(i-1))-1)
!>
!>  \author       Jose Luis Martins
!>  \version      6.013
!>  \date         22 June 2021
!>  \copyright    GNU Public License v2

subroutine atom_atm_setmesh(iowrite, aa, bb, znuc, rmax,                 &
           nr, r, drdi, d2rodr,                                          &
           mxdnr)

! written 18 July 2019. JLM
! jlm  version 6.013


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points

  real(REAL64), intent(in)          ::  aa, bb                           !<  grid parameters, r(i) =(exp(-aa)/znuc)*(exp((i-1)/bb)-1)
  real(REAL64), intent(in)          ::  znuc                             !<  nuclear charge
  real(REAL64), intent(in)          ::  rmax                             !<  r(i) does not exceed rmax

! output

  integer, intent(out)              ::  nr                               !<  dimension of the number of radial points

  real(REAL64), intent(out)         ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(out)         ::  drdi(mxdnr)                      !<  d r(i) / d i
  real(REAL64), intent(out)         ::  d2rodr(mxdnr)                    !<  (d^2 r(i) / d i^2) / (d r / d i)

! local variables

  real(REAL64)       ::  a, b                                            !<  r(i) = a*(exp(b*(i-1))-1)
  integer            ::  imax

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64

! counters

  integer     ::  i

! sanity checks

  if(bb <= ZERO) then
    write(iowrite,'("   atm_setmesh, bb = ",g12.3," < 0")') bb

    stop

  endif

  if(znuc <= ZERO) then
    write(iowrite,'("   atm_setmesh, znuc = ",g12.3," < 0")') znuc

    stop

  endif
  if(rmax < ZERO) then
    write(iowrite,'("   atm_setmesh, rmax = ",g12.3," < 0")') rmax

    stop

  endif

  a = exp(-aa)/znuc
  b = 1/bb

  do i = 1,mxdnr
    r(i) = a*(exp(b*(i-1))-1)
    drdi(i) = (r(i)+a)*b
    d2rodr(i) = b
  enddo

  imax = nint(log(rmax/a+ONE)/b+ONE)

  if(r(mxdnr) <= rmax) then
    if(a*(exp(b*mxdnr)-1) > rmax) then
      nr = mxdnr
    else
      write(iowrite,*)
      write(iowrite,'(" Error in atm_setmesh - arraylimits",        &
          & " for radial array exceeded",/)')
      write(iowrite,'(" Try mxdnr = ",i10)') imax + 10

      stop                                                          !  exits program with error

    endif
  else
    do i = min(2,imax-2),max(nr,imax+2)
      nr = i-1
      if (r(i) > rmax) exit
    enddo
  endif

  return

end subroutine atom_atm_setmesh
