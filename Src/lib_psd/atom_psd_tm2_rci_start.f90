!>  Calculates the error and jacobian of the equations of the
!>  improved scheme of N. Troullier and J. L. Martins
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.0.6
!>  \date         3 November 2021.
!>  \copyright    GNU Public License v2

subroutine atom_psd_tm2_rci_start(rjrc, polydrc, bkrk)


 implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter          ::  ntm2 = 6                              !  order of tm2 polynomial

! input

  real(REAL64), intent(in)          ::  rjrc                             !<  r(jrc)
  real(REAL64), intent(in)          ::  polydrc(0:ntm2-2)                !<  n-th derivatives of Kerker polynomial ar rc

! output

  real(REAL64), intent(out)         ::  bkrk(0:ntm2)                     !<  coefficient of polynomial

! local variables

! local variables

  integer          ::  info

! local allocatable arrays

  real(REAL64), allocatable         ::  rcpn(:)                          !  rc^n
  real(REAL64), allocatable         ::  aj(:,:)                          !  linear equation matrix
  integer, allocatable              ::  ipvt(:)

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64

! counter

  integer      ::  j, k


! setup matrix of linear equation

  allocate(rcpn(0:2*ntm2))
  allocate(aj(0:ntm2-2,0:ntm2-2))
  allocate(ipvt(0:ntm2-2))

  rcpn(0) = ONE
  rcpn(1) = rjrc
  do j = 2,2*ntm2
    rcpn(j) = rcpn(j-1)*rjrc
  enddo

  do j = 0,ntm2-2
    aj(0,j) = ONE
    do k = 1,ntm2-2
      aj(k,j) = ZERO
    enddo
  enddo
  do k = 1,ntm2-2
    do j = (k+1)/2,ntm2-2
      aj(k,j) = (2*j-k+1)*aj(k-1,j)
    enddo
  enddo

  do k = 0,ntm2-2
    do j = (k+1)/2,ntm2-2
      aj(k,j) = aj(k,j)*rcpn(2*j-k)
    enddo
  enddo

! solve linear system

  do k = 0,ntm2-2
    bkrk(k) = polydrc(k)
  enddo

  call dgefa(aj, ntm2-1, ntm2-1, ipvt, info)
  if(info /= 0) stop
  call dgesl(aj, ntm2-1, ntm2-1, ipvt, bkrk, 0)

! last two coefficients are zero

  bkrk(ntm2-1) = ZERO
  bkrk(ntm2) = ZERO

  deallocate(rcpn)
  deallocate(aj)
  deallocate(ipvt)

  return

end subroutine atom_psd_tm2_rci_start
