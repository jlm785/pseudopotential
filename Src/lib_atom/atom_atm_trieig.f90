!>  computes the eigenvalues and eigenvectors of a tridiagonal matrix
!>  by calling LAPACK subroutines
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.4
!>  \date         12 October 2021.
!>  \copyright    GNU Public License v2

subroutine atom_atm_trieig(n, d, sd, m, e, z, iowrite, ldz)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  ldz                              !<  leading dimension of z

  integer, intent(in)               ::  n                                !<  dimension of matrix
  real(REAL64), intent(in)          ::  d(n)                             !<  diagonal of matrix
  real(REAL64), intent(in)          ::  sd(n)                            !<  sub-diagonal of matrix sd(1) arbitrary
  integer, intent(in)               ::  m                                !<  number of lowest eigenvalues

! output

  real(REAL64), intent(out)         ::  e(m)                             !<  eigenvalues
  real(REAL64), intent(out)         ::  z(ldz,m)                         !<  eigenvectors

! allocatable arrays

  integer, allocatable              ::  iblock(:)
  integer, allocatable              ::  isplit(:)
  real(REAL64), allocatable         ::  work(:)
  integer, allocatable              ::  iwork(:)
  integer, allocatable              ::  ifail(:)

! local variables

  real(REAL64)   ::  dum
  real(REAL64)   ::  eps                              !<  requested precision (if negative-> machine precision)
  integer        ::  mm, nsplit, info

  allocate(iblock(n))
  allocate(isplit(n))
  allocate(work(5*n))
  allocate(iwork(3*n))
  allocate(ifail(n))

! negative eps makes LAPACK use machine precision.  sd convention
! is different between eispack and LAPACK

  eps = -1.0

  call dstebz('I', 'B', n, dum, dum, 1, m, eps, d, sd(2:n),              &
      mm, nsplit, e, iblock, isplit, work, iwork, info )

  if (info /= 0) write(iowrite,'(/,"  atom_atm_trieig:  error in",       &
      &       " dstebz       info =",i3,/)') info

  call dstein( n, d, sd(2:n), m, e, iblock, isplit, z, ldz, work,        &
      iwork, ifail, info )

  if (info /= 0) write(iowrite,'(/,"  atom_atm_trieig:  error in",       &
      &       " dstein       info =",i3,/)') info

  deallocate(iblock)
  deallocate(isplit)
  deallocate(work)
  deallocate(iwork)
  deallocate(ifail)

  return

end subroutine atom_atm_trieig
