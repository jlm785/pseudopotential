!>  Fourier transforms radial functions
!>  Tail is integrated with a spline interpolation.
!>
!>  \author       Sverre Froyen, Norm Troullier, JL Martins
!>  \version      6.0.4
!>  \date         80s, April 2012, 28 September 2021.
!>  \copyright    GNU Public License v2

subroutine ft_radial(nr, nq, r, drdi, qp, fin, fout, l)

! included the (2 pi)**3 factor here.

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  nr                               !<  number of points in the radial grid
  integer, intent(in)               ::  nq                               !<  number of points for the non-local Fourier grid

  real(REAL64), intent(in)          ::  r(0:nr)                          !<  radial grid points   r(i) = a*(exp(b*i)-1)
  real(REAL64), intent(in)          ::  drdi(0:nr)                       !<  d r(i) / d i
  real(REAL64), intent(in)          ::  qp(0:nq)                         !<  fourier grid points

  real(REAL64), intent(in)          ::  fin(0:nr)                        !<  function to be transformed
  integer, intent(in)               ::  l                                !<  angular momentum

! output

  real(REAL64), intent(out)         :: fout(0:nq)                        !<  Fourier transform of function

! allocatable arrays

  real(REAL64), allocatable  ::  a(:)
  real(REAL64), allocatable  ::  y(:), yp(:), ypp(:)                     !  function and derivatives in spline original grid
  real(REAL64), allocatable  ::  r2(:)                                   !  interpolated coordinates
  real(REAL64), allocatable  ::  w(:,:)                                  !  spline work (LU decomp)
  real(REAL64), allocatable  ::  v(:), vp(:), vpp(:)                     !  interpolated function and derivatives

! local variables

  real(REAL64)  ::  deltar
  real(REAL64)  ::  fout2j
  real(REAL64)  ::  fmax

  integer       ::  ierr, kerr
  integer       ::  nrmax, nr2, nr4

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  real(REAL64), parameter    ::  PI = 4*atan(ONE)
  real(REAL64), parameter    ::  EPS = 1.0E-12_REAL64
  real(REAL64), parameter    ::  XN = SQRT(2*ONE / PI)

! counters

  integer     :: j, k

! functions

  real(REAL64), external   :: sbessj


  do j = 0,nq
    fout(j) = ZERO
  enddo

! finds nrmax

  fmax = abs(fin(1)*drdi(1))
  do k = 1,nr
    if(abs(fin(k)*drdi(k)) > fmax) fmax = abs(fin(k)*drdi(k))
  enddo

  do k = nr,1,-1
    nrmax = k

    if(abs(fin(k)*drdi(k)) > fmax*EPS) exit

  enddo

! finds nr2

  deltar = 0.1/qp(nq)

  nr2 = nrmax
  do k = 1,nrmax
    if (r(k)-r(k-1) > deltar) then
      nr2 = k

      exit

    endif
  enddo

  nr2 = 7*(nr2/7)

! finds nr4

  if(nrmax <= nr2) then
    nr4 = 0
  else
    nr4 = int((r(nrmax)-r(nr2))/deltar)
    nr4 = 7*(nr4/7+1)+1
  endif

  if(nr2 > 0) then

    allocate(a(0:nr2))

!   Due to the high number of oscilations in the integrand,
!   an eight point Newton-Cotes intagration method is used.
!   See  Abramowitz and Stegun Eq. 25.4.17

    do j = 0,nq
      a(0) = ZERO
      do k = 1,nr2
        a(k) = fin(k)*r(k)*drdi(k)*sbessj(l,qp(j)*r(k))
      enddo
      do k = 0,nr2-7,7
        fout(j) = fout(j) + 751*(a(k)+a(k+7)) + 3577*(a(k+1)+a(k+6)) +   &
                            1323*(a(k+2)+a(k+5)) + 2989*(a(k+3)+a(k+4))
      enddo
      fout(j) = (7*XN*fout(j))/17280
    enddo

    deallocate(a)

  endif

! second part using splines

  if(nr4 > 7) then

    allocate(y(0:nr),yp(0:nr),ypp(0:nr))
    allocate(w(0:nr,3))

    y(0) = zero
    do k = 1,nr
      y(k) = fin(k)*r(k)
    enddo

    ierr = 0
    call splift(r, y, yp, ypp, nr+1, w, ierr, 0,-ONE/2,-ONE/2,ZERO,ZERO)
    if(ierr /= 1) then
      write(6,*) ' error in splift in ft_radial '

      stop

    endif

    deallocate(w)


    allocate(r2(nr4))
    allocate(v(nr4),vp(nr4),vpp(nr4))

    do k = 1,nr4
      r2(k) = r(nr2) + (k-1)*deltar
    enddo

    call splint(r, y, ypp, nr, r2, v, vp, vpp, nr4, kerr)

    deallocate(y,yp,ypp)
    deallocate(vp,vpp)

    allocate(a(nr4))

    do j = 0,nq
      fout2j = ZERO
      do k = 1,nr4
        a(k) = v(k)*sbessj(l,qp(j)*r2(k))
      enddo
      do k = 1,nr4-7,7
        fout2j = fout2j + 751*(a(k)+a(k+7)) + 3577*(a(k+1)+a(k+6)) +     &
                          1323*(a(k+2)+a(k+5)) + 2989*(a(k+3)+a(k+4))
      enddo
      fout(j) = fout(j) + (7*deltar*XN*fout2j) / 17280
    enddo

  endif

  deallocate(r2)
  deallocate(a)
  deallocate(v)

  return

end subroutine ft_radial
