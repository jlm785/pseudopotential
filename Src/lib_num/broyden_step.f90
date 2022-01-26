!>  Solution of a vectorial equation fvec = 0 using a Broyden method with a Jacobian.
!>  It is a reverse communication interface (RCI). Given an argument x_old and fvec(x_old)
!>  suggests a new argument x_new.
!>  Use for a moderate number of unknowns n.

!>  Rewrite with other values of EPS, TOL, STPMAX if not using natural units.

!   The typical use of this RCI for a function "funcv" would be

!   imax is the maximum number of iterations.  
!
!   do i = 1,imax
! 
!     call funcv(n,x,fvec)
! 
!     if(lnewjac) then
! 
!!      gets new jacobian
! -------------------------------------------------------------------------
!       call fdjac(n,x,fvec,NP,r)      !  jacobian by finite differences
!
!--------------------       OR  (but not both)   --------------------------
!
!       call funcv_jac(r,n,x)          !  analytical jacobian
!--------------------------------------------------------------------------       
!     endif
!     
!     call broyden_step(n,x,fvec,r,np,lnewjac,ireturn)
! 
!     if(ireturn /= 0 .and. ireturn /= 9) exit
! 
!   end do


subroutine broyden_step(n, x, fvec, tolf, r, np, lnewjac, ireturn)

! written May/June 2018. JLM
! Documentation 3 November 2020. JLM
! Copyright Jose Luis Martins/INESC-MN

! version 4.98

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)                ::  n                               !<  active space dimension
  integer, intent(in)                ::  np                              !<  leading dimension of r, np >= n

  real(REAL64), intent(in)           ::  fvec(n)                         !<  vectorial function

  real(REAL64), intent(in)           ::  tolf                            !<  convergence if abs(fvec(i)) < tolf for all i

! input and output

  real(REAL64), intent(inout)        ::  x(n)                            !<  argument of function
  real(REAL64), intent(inout)        ::  r(np,np)                        !<  original jacobian/right matrix of QR
  logical, intent(inout)             ::  lnewjac                         !<  new jacobian calculated (input)/ needed (output)

! output

  integer, intent(out)               ::  ireturn                         !<  status of calculation

! if(ireturn == 0) then
!   write(6,*) '  maximum number of iterations was exceeded'
! elseif(ireturn == 1) then
!   write(6,*) '  solution found after jacobian step'
!   write(6,*) '  guess is solution? number of iterations = ', its
! elseif(ireturn == 2) then
!   write(6,*) '  solution found after linesearch step'
!   write(6,*) '  number of iterations = ', its
! elseif(ireturn == 7) then
!   write(6,*) '  convergence in delta x'
! elseif(ireturn == 8) then
!   write(6,*) '  linesearch failed'
! elseif(ireturn == 9) then
!   write(6,*) '  needs new jacobian. its = ', its
!   write(6,*) '  maximum number of iterations was exceeded?'
! elseif(ireturn == 11) then
!   write(6,*) '  new jacobian did not solve the problem'
! elseif(ireturn == 12) then
!   write(6,*) '  spurious convergence. grad f = 0'
! elseif(ireturn == 20) then
!   write(6,*) '  singular Jacobian in broyden_step'
! elseif(ireturn == 21) then
!   write(6,*) '  r singular in broydn'
! endif

  
! allocatable local arrays.  It is not threadsafe.
  
  real(REAL64), allocatable, save    ::  fvcold(:)                       !  old fvec
  real(REAL64), allocatable, save    ::  xold(:)                         !  old x
  
  real(REAL64), allocatable, save    ::  qt(:,:)                         !  Q^T
 
  real(REAL64), allocatable, save    ::  g(:)                            !  gradient direction
  real(REAL64), allocatable, save    ::  p(:)                            !  linear search direction

  real(REAL64), allocatable, save    ::  c(:), d(:)                      !  denominator and diagonal of QR
  real(REAL64), allocatable, save    ::  s(:), t(:), w(:)                !  work arrays

! other local variables.  It is not threadsafe.

  logical                 ::  restrt, sing, skip
  logical                 ::  check
 
  integer, save           ::  istatus                                    !  linesearch status
  logical                 ::  lexitline

  real(REAL64)            ::  fmin                                       !  half norm squared of fvec
  real(REAL64), save      ::  fold                                       !  old value of half norm squared of fvec

  real(REAL64), save      ::  stpmax                                     !  maximum step in linesearch

  real(REAL64)            ::  den, xsum, temp, test

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64

  real(REAL64), parameter    ::  EPS = 3.0E-16_REAL64
  real(REAL64), parameter    ::  TOLMIN = 1.E-12_REAL64, TOLX = EPS
  real(REAL64), parameter    ::  STPMX = 100.0_REAL64

! counters

  integer         ::  i, j, k

! half norm squared of fvec

  fmin = ZERO
  do i = 1,n
    fmin = fmin + fvec(i)*fvec(i)
  enddo
  fmin = fmin/2


  if(lnewjac) then

    lexitline = .TRUE.

!   a new jacobian was calculated
    

    test = ZERO
    do i = 1,n
      if(abs(fvec(i)) > test) test = abs(fvec(i))
    enddo
    if(test < 0.01*TOLF) then
      check = .FALSE.

!     solution found

      ireturn = 1

      return

    endif

!   allocate arrays

    if(allocated(fvcold)) deallocate(fvcold)
    allocate(fvcold(np))
    if(allocated(xold)) deallocate(xold)
    allocate(xold(np))

    if(allocated(qt)) deallocate(qt)
    allocate(qt(np,np))

    if(allocated(g)) deallocate(g)
    allocate(g(np))
    if(allocated(p)) deallocate(p)
    allocate(p(np))

    if(allocated(c)) deallocate(c)
    allocate(c(np))
    if(allocated(d)) deallocate(d)
    allocate(d(np))
    if(allocated(s)) deallocate(s)
    allocate(s(np))
    if(allocated(t)) deallocate(t)
    allocate(t(np))
    if(allocated(w)) deallocate(w)
    allocate(w(np))

!   estimates stpmax

    xsum = ZERO
    do i = 1,n
      xsum = xsum + x(i)*x(i)
    enddo
    stpmax = STPMX*max(sqrt(xsum), n*UM)

!   QR decomposition

    call broyden_qr_decomp(r,n,NP,c,d,sing)

    if(sing) then

!     singular Jacobian in broyden_step

      ireturn = 20

      return

    endif

!   calculates Q^T

    do i = 1,n
      do j = 1,n
        qt(i,j) = ZERO
      enddo
      qt(i,i) = UM
    enddo

    do k = 1,n-1
      if(c(k) /= ZERO)then
        do j = 1,n
          xsum = ZERO
          do i = k,n
            xsum = xsum + r(i,k)*qt(i,j)
          enddo
          xsum = xsum/c(k)
          do i = k,n
            qt(i,j) = qt(i,j) - xsum*r(i,k)
          enddo
        enddo
      endif
    enddo

    do i = 1,n
      r(i,i) = d(i)
      do j = 1,i-1
        r(i,j) = ZERO
      enddo
    enddo


  else


!   jacobian already exists continue linesearch

 
    call broyden_linesearch(n,xold,fold,g,p,x,fmin,stpmax,istatus)


    if(istatus < 3) then
      lexitline = .TRUE.
    else
      lexitline = .FALSE.
    endif

    if(lexitline) then

!     linesearch step terminated successfuly or failed completely

      check = .FALSE.
      if(istatus < 0) check = .TRUE.

      test = ZERO
      do i = 1,n
        if(abs(fvec(i)) > test) test = abs(fvec(i))
      enddo

      if(test < TOLF)then
        check = .FALSE.

!       solution found

        ireturn = 2

        return

      endif

      if(check)then

!       linesearch had problems

        if(lnewjac) then

!         new jacobian didn't solve the problem.  This should never happen!  

          ireturn = 11

          return

        else

!         tests for new jacobian

          test = ZERO
          den = max(fmin, (UM*n)/2 )

          do i = 1,n
            temp = abs(g(i))*max(abs(x(i)),UM)/den
            if(temp > test) test = temp
          enddo

          if(test < TOLMIN) then

!           spurious convergence

            ireturn = 12

            return

          else

!           needs new jacobian

            lnewjac = .TRUE.
            ireturn = 9
           
            return

          endif

        endif

      else

!       linesearch was successful

        restrt = .FALSE.

!       test for same point

        test = ZERO
        do i = 1,n
          temp = (abs(x(i)-xold(i)))/max(abs(x(i)),UM)
          if(temp > test) test = temp
        enddo

        if(test < TOLX) then

!         convergence in delta x

          ireturn = 7

          return

        endif

      endif
               
!     update the jacobian


      do i = 1,n
        s(i) = x(i) - xold(i)
      enddo

      do i = 1,n
        xsum = ZERO
        do j = i,n
          xsum = xsum + r(i,j)*s(j)
        enddo
        t(i) = xsum
      enddo

      skip = .TRUE.

      do i = 1,n
        xsum = ZERO
        do j = 1,n
          xsum = xsum + qt(j,i)*t(j)
        enddo
        w(i) = fvec(i) - fvcold(i) - xsum

        if(abs(w(i)) >= EPS*(abs(fvec(i))+abs(fvcold(i)))) then
          skip = .FALSE.
        else
          w(i) = ZERO
        endif

      enddo

      if(.NOT. skip)then

!       QR update

        do i = 1,n
          xsum = ZERO
          do j = 1,n
            xsum = xsum + qt(i,j)*w(j)
          enddo
          t(i)=xsum
        enddo

        den = ZERO
        do i = 1,n
          den = den + s(i)*s(i)
        enddo
        do i = 1,n
          s(i) = s(i)/den
        enddo

        call broyden_qr_update(r,qt,n,NP,t,s)

        do i = 1,n
          if(r(i,i) == ZERO) then

!            r singular in broyden

             ireturn = 21

             return

          endif
          d(i) = r(i,i)
        enddo

      endif

    endif                                                                !   end of lexitline if

  endif                                                                  ! end of lnewjac if/then/else

  if(lexitline) then

    lnewjac = .FALSE.

!   calculates new gradient

    do i = 1,n
      xsum = ZERO
      do j = 1,n
        xsum = xsum + qt(i,j)*fvec(j)
      enddo
      g(i) = xsum
    enddo

    do i = n,1,-1
      xsum = ZERO
      do j = 1,i
        xsum = xsum + r(j,i)*g(j)
      enddo
      g(i) = xsum
    enddo

!   copies input results to "old" variables

    do i = 1,n
      xold(i) = x(i)
      fvcold(i) = fvec(i)
    enddo
    fold = fmin

!   calculates search direction

    do i = 1,n
      xsum = ZERO
      do j = 1,n
        xsum = xsum + qt(i,j)*fvec(j)
      enddo
      p(i) = -xsum
    enddo

    call broyden_rsolv(r,n,NP,d,p)

!   do full step with new g, p 

    istatus = 0


    call broyden_linesearch(n,xold,fold,g,p,x,fmin,stpmax,istatus)


    if(istatus == -1) then

!     linesearch failed

      check = .TRUE.
      ireturn = 8

      return

    endif

  endif

! continue the search, which seems to be working...

  ireturn = 0

  return
end subroutine broyden_step

