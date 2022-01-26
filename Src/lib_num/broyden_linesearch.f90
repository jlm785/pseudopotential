!>  Reverse communication interface for the linesearch algorithm A6.3.1 of
!>  Numerical Methods for Unconstrained Optimization and Nonlinear Equations
!>  JEDennis and RBSchnabel.  See also Numerical Recipes lnsrch

subroutine broyden_linesearch(n,xold,fold,g,p,xnew,fnew,stpmax,istatus)

! written May/June 2018. JLM
! Documentation 3 November 2020. JLM
! Copyright Jose Luis Martins/INESC-MN

! version 4.98

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  n                                !<  space dimension

  real(REAL64), intent(in)          ::  xold(n)                          !<  x_in
  real(REAL64), intent(in)          ::  fold                             !<  f(x_in)

  real(REAL64), intent(in)          ::  g(n)                             !<  grad f(x_in)

  real(REAL64), intent(in)          ::  fnew                             !<  f(x_out)

  real(REAL64), intent(in)          ::  stpmax                           !<  maximum allowed step

! output

  real(REAL64), intent(out)         ::  xnew(n)                          !<  x_out

! input and output

  real(REAL64), intent(inout)       ::  p(n)                             !<  search direction (may be rescaled on output)
  integer, intent(inout)            ::  istatus                          !<  status of the calculation
!                                                                        !   1: success on first attempt; 2 success after backtracking;
!                                                                        !   3: try again; 0: setup run ; < 0 : failure

! local variables that must be saved (subroutine is NOT THREADSAFE!)

  real(REAL64),save      ::  slope
  real(REAL64),save      ::  alam, alam2, alamin
  real(REAL64),save      ::  fnew2

! local variables

  real(REAL64)     ::  a, b, disc, rhs1, rhs2
  real(REAL64)     ::  xsum, temp, test,tmplam

! external function

  real(REAL64), external            ::  fmin

! constant

  real(REAL64), parameter   ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  real(REAL64), parameter   ::  ALF = 1.0E-4_REAL64, TOLX = 3.0E-16_REAL64

! counters

  integer         ::  i

  if(istatus == 0) then

    istatus = 3

    xsum = ZERO
    do i = 1,n
      xsum = xsum + p(i)*p(i)
    enddo
    xsum = sqrt(xsum)
    if(xsum > stpmax)then
      do i = 1,n
        p(i) = p(i)*stpmax/xsum
      enddo
    endif

    slope = ZERO
    do i = 1,n
      slope = slope + g(i)*p(i)
    enddo

    if(slope >= ZERO) then

!     roundoff error or wrong g or p, exit with warning

      istatus = -1

!     exit

    else

      test = ZERO
      do i = 1,n
        temp = abs(p(i))/max(abs(xold(i)),UM)
        if(temp > test) test = temp
      enddo

      alamin = TOLX/test
      alam = UM

      do i = 1,n
        xnew(i) = xold(i) + alam*p(i)
      enddo

    endif

  else

!   fnew = fmin(xnew)
  
    if(alam < alamin)then

!     alam is too small, exit with warning

      do i = 1,n
        xnew(i) = xold(i)
      enddo

      istatus = -2

!      exit

    elseif(fnew <= fold + ALF*alam*slope) then

!     f decreased by reasonable amount, successful exit

      if(istatus == 3) then
        istatus = 1
      else
        istatus = 2
      endif
      

!      exit

    else

!     find lambda for backtrack

      if(istatus == 3) then

!       first backtrack tmplam ~ 1/2

        tmplam = -slope/(2*(fnew-fold-slope))

      else

        rhs1 = fnew - fold - alam*slope
        rhs2 = fnew2 - fold - alam2*slope
        a = (rhs1/(alam*alam) - rhs2/(alam2*alam2)) / (alam-alam2)
        b = (-alam2*rhs1/(alam*alam) + alam*rhs2/(alam2*alam2))     &
               / (alam-alam2)

        if(abs(a) < TOLX/3) then
!         if(a==0.d0)then
          tmplam = -slope/(2*b)
        else
          disc = b*b - 3*a*slope
          if(disc <= ZERO) then
            write(6,*) 'roundoff problem in linesearch'
            tmplam = alam / 2
          elseif(b <=0) then
            tmplam = (-b+sqrt(disc))/(3*a)
          else
            tmplam = -slope / (b+sqrt(disc))
          endif
        endif

        if(tmplam > alam/2 ) tmplam = alam/2
      endif

      istatus = 4

!     endif

      alam2 = alam
      fnew2 = fnew
      alam = max(tmplam,alam/10)

      do i = 1,n
        xnew(i) = xold(i) + alam*p(i)
      enddo
    
    endif

  endif

! nddo

! should never occur

!  if(istatus == -3) then
!    write(6,*) "  Max iterations exceeded in linesearch"
!  endif

  return

  end subroutine broyden_linesearch
