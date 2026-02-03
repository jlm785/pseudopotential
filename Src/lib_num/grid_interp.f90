!>  Interpolates from one grid to another grid using Lagrange interpolation.
!>  Returns also the derivatives of the function on the new grid
!>  and an estimate of the respective errors.
!>  Based on D. B. Hunter, The Computer Journal, 3, 270 (1961)
!>  with the small difference algorithm of Numerical Recipes Eq. 3.1.5
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.6 of atom
!>  \date         15 July 2021.
!>  \copyright    GNU Public License v3


subroutine grid_interp(n, nd, nin, xin, fin, nout, xout, fout, dymax, ifail, mxdout)


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter          ::  iowrite = 6                            !  default tape for writing

! input

  integer, intent(in)               ::  mxdout                           !<  dimension of points in the output (new) grid

  integer, intent(in)               ::  n                                !<  order of Lagrange interpolation
  integer, intent(in)               ::  nd                               !<  order of derivatives

  integer, intent(in)               ::  nin                              !<  number of points in the input (old) grid
  real(REAL64), intent(in)          ::  xin(nin)                         !<  points on the input grid
  real(REAL64), intent(in)          ::  fin(nin)                         !<  f(x(i))

  integer, intent(in)               ::  nout                             !<  number of points in the output (new) grid
  real(REAL64), intent(in)          ::  xout(mxdout)                     !<  points on the output (new) grid

! output

  real(REAL64), intent(out)         ::  fout(mxdout,0:nd)                !<  fout(:,0): value of the function on the new grid; fout(:,n): n-th derivative
  real(REAL64), intent(out)         ::  dymax(0:nd)                      !<  estimate of the accuracy (from n-1 interpolation)
  integer, intent(out)              ::  ifail                            !<  ifail =/= 0 indicates failure of the process.

! allocatable arrays

  real(REAL64), allocatable         ::  xs(:)                            !  local old shifted grid
  real(REAL64), allocatable         ::  ys(:)                            !  function on old shifted grid

  real(REAL64), allocatable         ::  y(:)                             !  f(0):  interpolated value at x=0;  f(n):  interpolated derivative at x=0.
  real(REAL64), allocatable         ::  dy(:)                            !  error estimate (last correction)

  integer, allocatable              ::  isin(:), isout(:)                !  indexation of grids

! local variables

  logical        ::  linup, lindw
  logical        ::  loutup, loutdw

  real(REAL64)   ::  xdif
  integer        ::  near
  integer        ::  nn
  real(REAL64)   ::  xgmin, xgmax
  integer        ::  ninf, nsup

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  EPS = 1.0E-9_REAL64


! counters

  integer     ::  i, j


! check grid sizes and if original grid is ordered

  ifail = 0
  if(nin < 2 .or. nin < n+1) then
    write(iowrite,*) '  grid_interp,  nin = ', nin, '  n = ',n

    ifail = 1
    return

  endif

  if(n < 1) then
    write(iowrite,*) '  grid_interp,  n = ', n

    ifail = 2
    return

  endif

  if(nout > mxdout) then
    write(iowrite,*) '  grid_interp,  nout = ', nout, '  mxdout = ', mxdout

    ifail = 3
    return

  endif

  if(n > 10) then
    write(iowrite,*) '  grid_interp,  n = ', n
    write(iowrite,*) '  Do you know what you are doing? '
  endif

  xgmin = xin(1)
  xgmax = xin(1)
  do i = 2,nin
    if(xgmin > xin(i)) xgmin = xin(i)
    if(xgmax < xin(i)) xgmax = xin(i)
  enddo

  xdif = EPS*(xgmax - xgmin)
  do j = 1,nout
    if(xout(j) < xgmin-xdif .or. xout(j) > xgmax+xdif) then

    write(iowrite,*) '  grid_interp,  extrapolation not allowed'
    write(iowrite,*) '  xout(j),xgmin,xgmax = ', xout(j),xgmin,xgmax

    ifail = 4
    return

    endif
  enddo

! check if the input grid is ordered.  Close points are not allowed

  linup = .TRUE.
  do i = 1,nin-1
    if(xin(i+1) < xin(i)+xdif) then
      linup = .FALSE.

      exit

    endif
  enddo

  lindw = .TRUE.
  do i = 1,nin-1
    if(xin(i+1)+xdif > xin(i)) then
      lindw = .FALSE.

      exit

    endif
  enddo

! check if the output grid is ordered

  loutup = .TRUE.
  do i = 1,nout-1
    if(xout(i+1) < xout(i)) then
      loutup = .FALSE.

      exit

    endif
  enddo

  loutdw = .TRUE.
  do i = 1,nout-1
    if(xout(i+1) > xout(i)) then
      loutdw = .FALSE.

      exit

    endif
  enddo

! allocate arrays

  allocate(xs(0:n),ys(0:n))
  allocate(y(0:nd),dy(0:nd))

  do i = 0,nd
    dymax(i) = ZERO
  enddo

  if( (linup .or. lindw) .and. (loutup .or. loutdw) ) then

    if(  (linup .and. loutup) .or. (lindw .and. loutdw) ) then
      near = 1
    else
      near = nin
    endif

    do j = 1,nout

!     find nearest point below

      if(linup .and. loutup) then
        nn = near
        do i = near,nin
          if(xin(i) > xout(j)) exit
          nn = i
        enddo
        near = nn
      elseif(lindw .and. loutup) then
        nn = near
        do i = near,1,-1
          nn = i
          if(xin(i) > xout(j)) exit
        enddo
        near = nn
      elseif(linup .and. loutdw) then
        nn = near
        do i = near,1,-1
          nn = i
          if(xin(i) < xout(j)) exit
        enddo
        near = nn
      elseif(lindw .and. loutdw) then
        nn = near
        do i = near,nin
          if(xin(i) < xout(j)) exit
          nn = i
        enddo
        near = nn
      endif
      ninf = near - n/2
      if(ninf < 1) ninf = 1
      nsup = ninf + n
      if(nsup > nin) then
        nsup = nin
        ninf = nsup - n
      endif

      if(linup) then
        do i = 0,n
          xs(i) = xin(ninf+i) - xout(j)
          ys(i) = fin(ninf+i)
        enddo
      else
        do i = 0,n
          xs(i) = xin(nsup-i) - xout(j)
          ys(i) = fin(nsup-i)
        enddo
      endif

      call poly_interp(y, dy, xs, ys, n, nd)

      do i = 0,nd
        fout(j,i) = y(i)
        if(dymax(i) < dy(i)) dymax(i) = dy(i)
      enddo

    enddo

  else

!   one or both of the grids is not ordered

    allocate(isin(nin),isout(nout))

    if(linup) then
      do j = 1,nin
        isin(j) = j
      enddo
    else
      call sort(nin,xin,isin)
    endif

    if(loutup) then
      do j = 1,nout
        isout(j) = j
      enddo
    else
      call sort(nout,xout,isout)
    endif

!   check if the input grid has close points

    do i = 1,nin-1
      if( xin(isin(i+1)) < xin(isin(i))+xdif ) then

        write(iowrite,*) '  grid_interp,  duplicate poins in input grid'
        write(6,*) 'xin(',isin(i+1),') = ',xin(isin(i+1))
        write(6,*) 'xin(',isin(i),') = ',xin(isin(i))
        ifail = 5

        return

      endif
    enddo

    near = 1

    do j = 1,nout

!     find nearest point below

      nn = near
      do i = near,nin
        if(xin(isin(i)) > xout(isout(j))) exit
        nn = i
      enddo
      near = nn

      ninf = near - n/2
      if(ninf < 1) ninf = 1
      nsup = ninf + n
      if(nsup > nin) then
        nsup = nin
        ninf = nsup - n
      endif

      do i = 0,n
        xs(i) = xin(isin(ninf+i)) - xout(isout(j))
        ys(i) = fin(isin(ninf+i))
      enddo

      call poly_interp(y, dy, xs, ys, n, nd)

      do i = 0,nd
        fout(isout(j),i) = y(i)
        if(dymax(i) < dy(i)) dymax(i) = dy(i)
      enddo

    enddo

  endif

  deallocate(xs,ys)
  deallocate(y,dy)

  return
  end subroutine grid_interp
