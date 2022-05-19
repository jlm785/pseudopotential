!>  Solves the equations of the
!>  improved scheme of N. Troullier and J. L. Martins
!>  using a reverse communication interface
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.8
!>  \date         4 November 2021. 19 May 2022.
!>  \copyright    GNU Public License v2

subroutine atom_psd_tm2_eqsolve_rci(r, drdi, jrc, lo,                    &
     polydrc, cdrc, bkrk,                                                &
     mxdnr)

! extracted from the atom_psd_tm2 code. 29 October 2021. JLM
! deallocate, 19 May 2022. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter          ::  ntm2 = 6                              !  order of tm2 polynomial

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension  of the number of radial points

  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  integer, intent(in)               ::  jrc                              !<  r(jrc) = r_cut

  real(REAL64), intent(in)          ::  cdrc                             !<  int_0^r_cut rho(r) dr
  real(REAL64), intent(in)          ::  polydrc(0:ntm2-2)                !<  n-th derivatives of Kerker polynomial ar rc

  integer, intent(in)               ::  lo                               !<  angular quantum number l

! output

  real(REAL64), intent(out)         ::  bkrk(0:ntm2)                     !<  coefficient of polynomial

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

  real(REAL64), parameter    ::  TOLF = 1.0E-12_REAL64

  integer, parameter ::   MAXITS = 500

! counter

  integer      ::  j, k, ITS


! Find starting point  without norm-conservation or V''(0) = 0

  call atom_psd_tm2_rci_start(r(jrc), polydrc, bkrk)

  allocate(fvec(0:ntm2))
  allocate(xj(0:ntm2,0:ntm2))
  allocate(xjac(0:ntm2,0:ntm2))

  lnewjac = .TRUE.
  njac = 0

  do its = 1,MAXITS

    call atom_psd_tm2_func(r, drdi, jrc, lo,                             &
        polydrc, cdrc, bkrk,                                             &
        errsq, fvec, xj,                                                 &
        mxdnr)

    if(lnewjac) then
      do j = 0,ntm2
      do k = 0,ntm2
        xjac(k,j) = xj(k,j)
      enddo
      enddo
      njac = njac+1
    endif

    call broyden_step(ntm2+1, bkrk, fvec, TOLF, xjac, ntm2+1, lnewjac, ireturn)

    if(ireturn /= 0 .and. ireturn /= 9) exit

  enddo

  if(ireturn == 0) then

    write(6,*)
    write(6,*) '  stopped in atom_psd_tm2_eqsolve_rci:'
    write(6,*) '  maximum number of iterations was exceeded', its, MAXITS
    write(6,*)
    write(6,'("   errsq = ",e13.3)') errsq

    STOP

  endif

  if(ireturn == 8 .or. ireturn > 9) then

    write(6,*)
    write(6,*) '  stopped in atom_psd_tm2_eqsolve_rci:'
    write(6,*) '  problem in broyden_step  ireturn = ', ireturn

    STOP

  endif

  allocate(fvec)
  allocate(xj)
  allocate(xjac)

  return

end subroutine atom_psd_tm2_eqsolve_rci
