!>  Solves the equations of the
!>  improved scheme of N. Troullier and J. L. Martins
!>  using a reverse communication interface
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.8
!>  \date         4 November 2021. 19 May 2022.
!>  \copyright    GNU Public License v2

subroutine atom_psd_mrpp_eqsolve_rci(r, drdi, d2rodr, jrc, lo,           &
     polydrc, cdrc, cdrc2, ar2, ar2p, ev1, ev2, bkrk,                    &
     iowrite, mxdnr)

! extracted from the atom_psd_tm2 code. 29 October 2021. JLM
! deallocation, 19 May 2022. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter          ::  MRPP = 8                              !  order of tm2 polynomial

! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdnr                            !<  dimension  of the number of radial points

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i
  real(REAL64), intent(in)          ::  d2rodr(mxdnr)                    !<  (d^2 r(i) / d i^2) / (d r / d i)

  integer, intent(in)               ::  jrc                              !<  r(jrc) = r_cut

  real(REAL64), intent(in)          ::  polydrc(0:MRPP-4)                !<  n-th derivatives of Kerker polynomial ar rc
  real(REAL64), intent(in)          ::  cdrc                             !<  int_0^r_cut rho(r) dr
  real(REAL64), intent(in)          ::  cdrc2                            !<  int_0^r_cut rho(r) dr second orbital

  real(REAL64), intent(in)          ::  ar2                              !<  r*psi(jrc) for second orbital
  real(REAL64), intent(in)          ::  ar2p                             !<  d r*psi / d r (jrc) for second orbital

  real(REAL64), intent(in)          ::  ev1                              !<  orbital energy of main orbital
  real(REAL64), intent(in)          ::  ev2                              !<  orbital energy of second orbital

  integer, intent(in)               ::  lo                               !<  angular quantum number l

! input and output

  real(REAL64), intent(inout)         ::  bkrk(0:MRPP)                   !<  coefficient of polynomial

! local variables

  real(REAL64)                      ::  errsq                            !  sum of square errors

  logical         ::  lnewjac
  integer         ::  njac
  integer         ::  ireturn

! allocatable arrays

  real(REAL64), allocatable         ::  fvec(:)                          !  error in the solution of equation
  real(REAL64), allocatable         ::  xj(:,:)                          !  jacobian
  real(REAL64), allocatable         ::  xjac(:,:)                        !  jacobian

! constants

  real(REAL64), parameter    ::  TOLF = 1.0E-8_REAL64

  integer, parameter ::   MAXITS = 100

! counter

  integer      ::  j, k, ITS



  allocate(fvec(0:MRPP))
  allocate(xj(0:MRPP,0:MRPP))
  allocate(xjac(0:MRPP,0:MRPP))

  lnewjac = .TRUE.
  njac = 0

  do its = 1,MAXITS

    call atom_psd_mrpp_func(r, drdi, d2rodr, jrc, lo,                    &
        polydrc, cdrc, cdrc2, ar2, ar2p, ev1, ev2, bkrk,                 &
        errsq, fvec, xj,                                                 &
        iowrite, mxdnr)

    if(lnewjac) then
      do j = 0,MRPP
      do k = 0,MRPP
        xjac(k,j) = xj(k,j)
      enddo
      enddo
      njac = njac+1
    endif

    call broyden_step(MRPP+1, bkrk, fvec, TOLF, xjac, MRPP+1, lnewjac, ireturn)

    if(ireturn /= 0 .and. ireturn /= 9) exit

  enddo

  if(ireturn == 0) then

    write(6,*)
    write(6,*) '  stopped in atom_psd_mrpp_eqsolve_rci:'
    write(6,*) '  maximum number of iterations was exceeded', its, MAXITS

    STOP

  endif

  if(ireturn == 8 .or. ireturn > 9) then

    write(6,*)
    write(6,*) '  stopped in atom_psd_mrpp_eqsolve_rci:'
    write(6,*) '  problem in broyden_step  ireturn = ', ireturn

    STOP

  endif

  deallocate(fvec)
  deallocate(xj)
  deallocate(xjac)

  return

end subroutine atom_psd_mrpp_eqsolve_rci
