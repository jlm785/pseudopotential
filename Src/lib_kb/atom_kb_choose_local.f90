!>  Chooses the local potential in KB transformation
!>
!>  \author       Norm Troullier, J.L.Martins
!>  \version      6.0.7
!>  \date         November 90, May 2012, July 2021, 8 December 2021. JLM
!>  \copyright    GNU Public License v2

subroutine atom_kb_choose_local(npot, lo, llocal ,nr, r, vionic, vlocal,   &
          mxdl, mxdnr)

! adapted from the old program jlm 22/5/2012 to 9/6/2012
! mxdl, mxdnr,  17 September 2021. JLM
! smoothing of local. 8 December 2021. JLM
! bug in smoothing of local potential. 26 January 2022. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum
  integer, intent(in)               ::  mxdnr                            !<  dimension of radial grid points.

  integer, intent(in)               ::  npot(-1:1)                       !<  number of orbitals to be calculated
  integer, intent(in)               ::  lo(mxdl+1,-1:1)                  !<  angular momentum of orbital

  integer, intent(in)               ::  nr                               !<  number of points in the radial grid
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points   r(i) = a*(exp(b*i)-1), i=1,...,nr
  real(REAL64), intent(in)          ::  vionic(mxdnr,0:mxdl,-1:1)        !<  r*pseudopotential
  integer, intent(in)               ::  llocal                           !<  choice of local potential (-1 is supremum)

! output

  real(REAL64), intent(out)         ::  vlocal(mxdnr)                    !<  r*local pseudopotential

! local variables

  real(REAL64), allocatable         ::  vk(:)
  real(REAL64), allocatable         ::  xi(:), yi(:), ypi(:), yppi(:)

  real(REAL64)               ::  vmax
  integer                    ::  jc
  logical                    ::  lbreak
  integer                    ::  iold, kold, inew, knew
  integer                    ::  ncross, jcross(20)
  logical                    ::  lcross(20)

  integer                    ::  jmin, jmax

  real(REAL64)               ::  x(4), y(4), yp(4), ypp(4), w(4,3)
  integer                    ::  ierr

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter    ::  EPS = 1.0E-6_REAL64
  integer, parameter         ::  JJ = 6                                  !  width of smoothing around a kink

! counters

  integer                  :: i, j, k

  do j = 1,nr
    vlocal(j) = ZERO
  enddo

  if (llocal >= 0 ) then

    do j = 2,nr
      vlocal(j) = vionic(j,llocal,0)
    enddo

  else

    allocate(vk(nr))

!   outside core radius for all

    lbreak = .FALSE.
    do j = nr,2,-1
      jc = j
      vmax = vionic(j,lo(npot(0),0),0)
      do k = -1,1
        if(npot(k) > 0) then
          do i = 1,npot(k)
            if (abs(vionic(j,lo(i,k),k) - vmax) > EPS) lbreak = .TRUE.
            if(lbreak) exit
          enddo
        endif
        if(lbreak) exit
      enddo
      if(lbreak) exit
    enddo

    jc = min(jc+5,nr-1)

    do j = jc+1,nr
      vk(j) = vionic(j,lo(npot(0),0),0) / r(j)
    enddo

!   inside core radius.  Identify crossings

    iold = 0
    kold = 0
    vmax = vionic(2,lo(npot(0),0),0) / r(2)
    do k = -1,1
      if(npot(k) > 0) then
        do i = 1,npot(k)
          if (vionic(2,lo(i,k),k) / r(2) > vmax + EPS) then
            vmax = vionic(2,lo(i,k),k) / r(2)
            iold = i
            kold = k
          endif
        enddo
      endif
    enddo
    vk(2) = vmax

    ncross = 0
    do j = 3,jc

      vmax = vionic(j,lo(iold,kold),kold) / r(j)
      inew = iold
      knew = kold

      do k = -1,1
        if(npot(k) > 0) then
          do i = 1,npot(k)
            if (vionic(j,lo(i,k),k) / r(j) > vmax + EPS) then
              vmax = vionic(j,lo(i,k),k)  / r(j)
              inew = i
              knew = k
            endif
          enddo
        endif
      enddo
      vk(j) = vmax

      if(inew /= iold .or. knew /= kold) then
        ncross = ncross+1
        if(ncross > 20) then
          write(6,*)
          write(6,*) '   WARNING:   Too many crossings in atom_kb_choose_local.'
          write(6,*) '   Check the local potential.'
          write(6,*)
        else
          jcross(ncross) = j
          iold = inew
          kold = knew
        endif
      endif

    enddo

!   smoothing of local potential with spline

    do i = 1,ncross
      lcross(i) = .TRUE.
    enddo

    do i = 1,ncross
      if(lcross(i)) then

        jmin = jcross(i) - JJ
        jmax = jcross(i) + JJ-1
        if(i /= ncross) then
          do k = i+1,ncross
            if(jcross(k) - jcross(k-1) > 2*JJ-3) then
              exit
            else
              lcross(k) = .FALSE.
              jmax = jcross(k) + JJ-1
            endif
          enddo
        endif

        x(1) = r(jmin)
        y(1) = vk(jmin)
        x(2) = r(jmin+1)
        y(2) = vk(jmin+1)
        x(3) = r(jmax-1)
        y(3) = vk(jmax-1)
        x(4) = r(jmax)
        y(4) = vk(jmax)

        call splift(x,y,yp,ypp,4,w,ierr,0,ZERO,UM/2,ZERO,UM/2)

        allocate(xi(jmax-jmin+1),yi(jmax-jmin+1))
        allocate(ypi(jmax-jmin+1),yppi(jmax-jmin+1))

        do j = jmin,jmax
          xi(j-jmin+1) = r(j)
        enddo

        call splint(x,y,ypp,4,xi,yi,ypi,yppi,jmax-jmin+1,ierr)

        do j = jmin,jmax
          vk(j) = yi(j-jmin+1)
        enddo

        deallocate(xi,yi,ypi,yppi)

      endif
    enddo

!   store result

    do j = 2,nr
      vlocal(j) = vk(j) * r(j)
    enddo

    deallocate(vk)

  endif

  return

end subroutine atom_kb_choose_local
