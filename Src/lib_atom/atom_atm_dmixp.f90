!>  Computes a new vector in a iterative scheme
!>  using Anderson's extrapolation scheme
!>  eqs 4.1-4.9,4.15-4.18 of
!>  D.G.Anderson J.Assoc.Computing Machinery,12,547(1965)
!>
!>  \author       Jose Luis Martins
!>  \version      6.013
!>  \date         1980s, 22 June 2021
!>  \copyright    GNU Public License v2

subroutine atom_atm_dmixp(a, b, beta, icy, id, nr, vmem, mxdnr)

!  *************************************************************
!  *
!  *    adapted from K.C.Pandey
!  *    input a=newpot b=oldpot
!  *    output a=a-b b=newpot
!  *    beta=mixing,in=iter. number
!  *    id=1,2 or 3 diff conv meth.
!  *    icy cycle number ,icy=1 on first/zeroth call
!  *    c,d work arrays of size nr
!  *    vmem storage arrays of size (nmsh,4)
!  *
!  *************************************************************

! converted to fortran 90, March/April 2018
! cleanup and new interface, July 2019. JLM
! jlm  version 6.00


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points

  integer, intent(in)               ::  nr                               !<  dimension of the number of radial points

  real(REAL64), intent(in)          ::  beta                             !<  mixing coefficient
  integer, intent(in)               ::  icy                              !<  cycle number, icy = 1 on first/zeroth call
  integer, intent(in)               ::  id                               !<  method for convergence (id = 1,2,3).

! input and output

  real(REAL64), intent(inout)       ::  a(mxdnr)                         !<  input: calculated potential, output:  correction
  real(REAL64), intent(inout)       ::  b(mxdnr)                         !<  input: old potential, output:  new potential

  real(REAL64), intent(inout)       ::  vmem(mxdnr,4)                    !<  save data for later.  Do not change outside subroutine!!!!

! work arrays

  real(REAL64), allocatable         ::  c(:)
  real(REAL64), allocatable         ::  d(:)
  real(REAL64), allocatable         ::  vn1(:)
  real(REAL64), allocatable         ::  vn12(:)
  real(REAL64), allocatable         ::  vn2(:)
  real(REAL64), allocatable         ::  vn22(:)

! local variables

  integer                    ::  in
  real(REAL64)               ::  r2
  real(REAL64)               ::  d11, d12, d22
  real(REAL64)               ::  rd1m, rd2m
  real(REAL64)               ::  t1, t2
  real(REAL64)               ::  x
  real(REAL64)               ::  bt1, bt2
  real(REAL64)               ::  a2
  real(REAL64)               ::  det, dett

! external

  real(REAL64), external            ::  ddot

! counter

  integer                    ::  i

! parameters

  real(REAL64), parameter           ::  ZERO = 0.0_REAL64
  real(REAL64), parameter           ::  UM = 1.0_REAL64
  real(REAL64), parameter           ::  DETOL = 1.0 / 10.0**9


  allocate(c(nr),d(nr))
  allocate(vn1(nr),vn2(nr),vn12(nr),vn22(nr))

  do i = 1,nr
    vn1(i)  = vmem(i,1)
    vn12(i) = vmem(i,2)
    vn2(i)  = vmem(i,3)
    vn22(i) = vmem(i,4)
  enddo

  in = icy - 1

  if(in == 0) then

    call daxpy(nr,UM,a,1,b,1)

  else

    call daxpy(nr,-UM,b,1,a,1)
    r2 = ddot(nr,a,1,a,1)

    if(id == 1) then

      call daxpy(nr,beta,a,1,b,1)

    else

      if(in == 1) then

        call dcopy(nr,a,1,vn1,1)
        call dcopy(nr,b,1,vn2,1)
        call daxpy(nr,beta,a,1,b,1)

      else

        call dcopy(nr,vn1,1,c,1)

        if(id == 3 .and. in > 2) then
          call dcopy(nr,vn12,1,d,1)
        endif

        call dcopy(nr,a,1,vn1,1)

        if(id > 2 .and. in > 1) then
          call dcopy(nr,c,1,vn12,1)
        endif

        call daxpy(nr,-UM,a,1,c,1)
        d11 = ddot(nr,c,1,c,1)
        rd1m = ddot(nr,a,1,c,1)

        if(in <=2 .or. id <=2) then

          t1 = -rd1m/d11
          x = UM - t1
          bt1 = beta*t1
          call dscal(nr,beta,a,1)
          call daxpy(nr,bt1,c,1,a,1)
          call dcopy(nr,vn2,1,d,1)
          call daxpy(nr,t1,d,1,a,1)
          call dcopy(nr,b,1,vn2,1)

          if(id > 2 .and. in == 2) then
            call dcopy(nr,d,1,vn22,1)
          endif

          do i = 1,nr
            b(i) = x*b(i) + a(i)
          enddo

        else

          call daxpy(nr,-UM,a,1,d,1)
          d22 = ddot(nr,d,1,d,1)
          d12 = ddot(nr,c,1,d,1)
          rd2m = ddot(nr,a,1,d,1)
          a2 = d11*d22
          det = a2 - d12*d12
          dett = det/a2

          if(abs(dett) > DETOL) then
            t1 = (-rd1m*d22 + rd2m*d12) / det
            t2 = ( rd1m*d12 - rd2m*d11) / det
          else
            t1 = -rd1m/d11
            t2 = ZERO
          endif

          x = UM - t1 - t2
          bt1 = beta*t1
          bt2 = beta*t2
          call dscal(nr,beta,a,1)
          call daxpy(nr,bt1,c,1,a,1)
          call daxpy(nr,bt2,d,1,a,1)
          call daxpy(nr,t1,vn2,1,a,1)
          call daxpy(nr,t2,vn22,1,a,1)
          call dcopy(nr,vn2,1,vn22,1)
          call dcopy(nr,b,1,vn2,1)
          do i = 1,nr
            b(i) = x*b(i) + a(i)
          enddo

        endif

      endif

    endif

  endif
  do i = 1,nr
    vmem(i,1) = vn1(i)
    vmem(i,2) = vn12(i)
    vmem(i,3) = vn2(i)
    vmem(i,4) = vn22(i)
  enddo

  deallocate(c,d)
  deallocate(vn1,vn2,vn12,vn22)

  return

end subroutine atom_atm_dmixp
