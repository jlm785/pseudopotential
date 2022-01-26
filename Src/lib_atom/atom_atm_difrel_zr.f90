!>  Finds the eigenvalue ev, the  minor and major
!>  components of the relativistic wavefunction ar and br (r u(r))
!>  for the boundary condition ar(revi) = 0.
!>  It is only reliable small revi (single digit, low tenths),
!>  or for scattering states ev ~ 0.
!>
!>  \author       J.L.Martins
!>  \version      6.0.7
!>  \date         September 2021.
!>  \copyright    GNU Public License v2

subroutine atom_atm_difrel_zr(nr, r, drdi, v, ar, br,                    &
        no, lo, iso, znuc, vzero, ev, iflag, revi, nrevi, tol,           &
        iowrite, mxdnr)

! uses a different method (secant+bissection) than old difnrl
! exits if bracketing energy range is small. 11 December 2021. JLM

  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdnr                            !<  maximum number of radial grid points.

  integer, intent(in)               ::  nr                               !<  number of points in the radial grid
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points   r(i) = a*(exp(b*i)-1), i=1,...,nr
  real(REAL64), intent(in)          ::  drdi (mxdnr)                     !<  (d r(i) / d i)

  real(REAL64), intent(in)          ::  v(mxdnr)                         !<  effective potential (RYDBERG) includes centrifugal potential

  integer, intent(in)               ::  no                               !<  principal quantum number of orbital
  integer, intent(in)               ::  lo                               !<  angular momentum of orbital
  integer, intent(in)               ::  iso                              !<  2*spin or 2*(j-l)

  real(REAL64), intent(in)          ::  znuc                             !<  nuclear charge
  real(REAL64), intent(in)          ::  vzero                            !<  effective potential in Rydberg

  real(REAL64), intent(in)          ::  revi                             !<  wave-function calculated up to revi

  real(REAL64), intent(in)          ::  tol                              !<  precision of the solution tol ~1.0e-10

! input and output

  real(REAL64), intent(inout)       ::  ev                               !<  guess eigenvalue on input eigenvalue on output

! output

  real(REAL64), intent(out)         ::  ar(mxdnr)                        !<  major component radial wave-function u = rR(r)
  real(REAL64), intent(out)         ::  br(mxdnr)                        !<  minor component (integral ar**2 + br**2 = 1)

  integer, intent(out)              ::  iflag                            !<  iflag = 0: success; iflag = 1: failed to converge, iflag > 3: major error

  integer, intent(out)              ::  nrevi                            !<  wave-functions was calculated up to r(nrevi)

! local variables

  real(REAL64)             ::  evi                                       !  only one eigenvalue needed

  real(REAL64)             ::  dev, dev0
  real(REAL64)             ::  x, x1, x2, f, f1, f2

  real(REAL64)             ::  xmin, xmax
  real(REAL64)             ::  fmin, fmax

  real(REAL64)             ::  temp
  integer                  ::  nrv
  integer                  ::  nodes
  logical                  ::  lfind1, lfind2, lfind3
  logical                  ::  lup

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64

! counters

  integer               ::  j, k


  iflag = 0
  lfind1 = .FALSE.
  lfind2 = .FALSE.
  lfind3 = .FALSE.


  evi =  ev
  dev0 = max(-evi/10, ONE/100)

! loop to find wave-function with no-lo-1 nodes

  dev = dev0

  do j = 1,100
    call atom_atm_difrel_one(nr, r, drdi, v, ar, br,                     &
        lo, iso, znuc, vzero, evi, revi, nrv, iflag,                     &
        iowrite, mxdnr)

    if(iflag /= 0) exit

    nrevi = nrv

    nodes = 0
    do k = 3,nrevi
      if(ar(k)*ar(k-1) < ZERO) nodes = nodes+1
    enddo

    if(nodes == no-lo-1) then

      x1 = evi
      f1 = ar(nrevi)
      lfind1 = .TRUE.
      xmin = x1
      fmin = f1

      exit

    elseif(nodes < no-lo-1) then

      if(j > 1) then
        if(lup) then
          dev = dev + dev/5
        else
          dev = dev / 2
        endif
      endif
      evi = evi + dev
      lup = .TRUE.

    else

      if(j > 1) then
        if(lup) then
          dev = dev / 2
        else
          dev = dev + dev / 5
        endif
      endif
      evi = evi - dev
      lup = .FALSE.

    endif
  enddo

  dev = dev0

  evi = evi + dev


  if(lfind1) then

!   loop to find wave-function with no-lo-1 nodes

    do j = 1,100
      call atom_atm_difrel_one(nr, r, drdi, v, ar, br,                   &
          lo, iso, znuc, vzero, evi, revi, nrv, iflag,                   &
          iowrite, mxdnr)

      if(iflag /= 0) exit

      nrevi = nrv

      nodes = 0
      do k = 2,nrevi
        if(ar(k)*ar(k-1) < ZERO) nodes = nodes+1
      enddo

      if(nodes == no-lo) then

        x2 = evi
        f2 = ar(nrevi)
        lfind2 = .TRUE.
        xmax = x2
        fmax = f2

        exit

      elseif(nodes < no-lo) then

        if(j > 1) then
          if(lup) then
            dev = dev + dev/5
          else
            dev = dev / 2
          endif
        endif
        evi = evi + dev
        lup = .TRUE.

      else

        if(j > 1) then
          if(lup) then
            dev = dev / 2
          else
            dev = dev + dev / 5
          endif
        endif
        evi = evi - dev
        lup = .FALSE.

      endif
    enddo

  endif


  if(lfind2) then

!   x2 closest to solution

    if(abs(f1) < abs(f2)) then
      temp = x1
      x1 = x2
      x2 = temp
      temp = f1
      f1 = f2
      f2 = temp
    endif

! secant method with some bissection

    lfind3 = .FALSE.
    do j = 1,100

      x = x2 - f2*(x1 - x2) / (f1 -f2)
      if(x < xmin .or. x > xmax) then
        if(f2*fmin < ZERO) then
          x = (xmin + x2) / 2
        else
          x = (xmax + x2) / 2
        endif
      endif

      call atom_atm_difrel_one(nr, r, drdi, v, ar, br,                   &
          lo, iso, znuc, vzero, x, revi, nrv, iflag,                     &
          iowrite, mxdnr)

      if(iflag /= 0) exit

      nrevi = nrv

      nodes = 0
      do k = 2,nrevi
        if(ar(k)*ar(k-1) < ZERO) nodes = nodes+1
      enddo

      f = ar(nrevi)

      if(nodes == no-lo-1) then
        xmin = x
        fmin = f
      elseif(nodes == no-lo) then
        xmax = x
        fmax = f
      endif

      if(abs(f) < tol .OR. abs(xmin - xmax)/max(abs(xmin),abs(xmax)) < tol) then
        lfind3 = .TRUE.
        exit
      endif

      x1 = x2
      f1 = f2
      x2 = x
      f2 = f
    enddo
  endif

  if(.not. lfind3) iflag = 5

  ev = x

  return

end subroutine atom_atm_difrel_zr
