!>  Finds total exchange-correlation energy and potential for a
!>  spherical electron density distribution.
!>  this version implements the local (spin) density approximation and
!>  the generalized-gradient-aproximation with the 'explicit mesh
!>  functional' method of White & Bird, prb 50, 4954 (1994).
!>
!>  \author       L.C.Balbas, J.M.Soler, Jose Luis Martins
!>  \version      6.013
!>  \date         22 June 2021
!>  \copyright    GNU Public License v2

subroutine atom_xc(functl, author, irel,                                 &
                  nr, maxr, rmesh, nspin, dens,                          &
                  ex, ec, dx, dc, vxc,                                   &
                  iowrite)

! gradients are 'defined' by numerical derivatives, using 2*NN+1 mesh
! points, where NN is a parameter defined below
! coded by L.C.Balbas and J.M.Soler. December 1996. version 0.5.

! character*(*) functl : functional to be used:
!              'lda' or 'lsd' => local (spin) density approximation
!                       'gga' => generalized gradient corrections
!                                uppercase is optional
! character*(*) author : parametrization desired:
!     'ca' or 'pz' => lsd perdew & zunger, prb 23, 5075 (1981)
!           'pw92' => lsd perdew & wang, prb, 45, 13244 (1992). this is
!                     the local density limit of the next:
!            'pbe' => gga perdew, burke & ernzerhof, prl 77, 3865 (1996)
!                     uppercase is optional

! ************************ units ************************************
! distances in atomic units (bohr).
! densities in atomic units (electrons/bohr**3)
! energy unit depending of parameter EUNIT below
! ********* routines called *****************************************
! xc_gg, xc_lda
! *******************************************************************

! converted to f90, April 2018. JLM
! cleanup, new interface, July 2019. JLM
! jlm  version 6.00


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! fix the order of the numerical derivatives: the number of radial
! points used is 2*NN+1

  integer, parameter   ::   NN = 5

! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  maxr                             !<  physical first dimension of rmesh, dens and vxcmaximum number of points
  integer, intent(in)               ::  nspin                            !<  nspin=1 => unpolarized; nspin=2 => polarized

  character(len=*), intent(in)      ::  functl                           !<  functional to be used
  character(len=*), intent(in)      ::  author                           !<  parametrization desired
  integer, intent(in)               ::  irel                             !<  relativistic exchange? (0 => no, 1 => yes)
  integer, intent(in)               ::  nr                               !<  number of radial mesh points

  real(REAL64), intent(in)          ::  rmesh(maxr)                      !<  radial mesh points
  real(REAL64), intent(in)          ::  dens(maxr,nspin)                 !<  total (nspin=1) or spin (nspin=2) electron density at mesh points

! output

  real(REAL64), intent(out)         ::  ex                               !<  total exchange energy
  real(REAL64), intent(out)         ::  ec                               !<  total correlation energy
  real(REAL64), intent(out)         ::  dx                               !<  integralof( rho * (eps_x - v_x) )
  real(REAL64), intent(out)         ::  dc                               !<  integralof( rho * (eps_c - v_c) )
  real(REAL64), intent(out)         ::  vxc(maxr,nspin)                  !<  (spin) exch-corr potential

! local variables and arrays

  logical              ::  gga
  real(REAL64)         ::  d(2), dexdd(2), decdd(2)
  real(REAL64)         ::  dexdgd(3,2), decdgd(3,2)
  real(REAL64)         ::  dgdm(-NN:NN), dgidfj(-NN:NN)
  real(REAL64)         ::  drdm, dvol
  real(REAL64)         ::  epsc, epsx
  real(REAL64)         ::  f1, f2, gd(3,2)
  real(REAL64)         ::  pi

  external ggaxc, ldaxc

! allocatable arrays

  real(REAL64), allocatable ::  aux(:)

! constants

  real(REAL64), parameter   ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64

! fix energy unit:  EUNIT=1.0 => hartrees,
!                   EUNIT=0.5 => rydbergs,
!                   EUNIT=0.03674903 => ev

  real(REAL64), parameter   ::  EUNIT = 0.5_REAL64

! dvmin is added to differential of volume to avoid division by zero

  real(REAL64), parameter   ::  DVMIN = 1.d-12

! counters

  integer         ::  in, in1, in2, ir, is, jn


! set gga switch

  if ( functl == 'lda' .or. functl == 'LDA' .or.                    &
       functl == 'lsd' .or. functl == 'LSD' ) then
    gga = .false.
  elseif ( functl == 'gga' .or. functl == 'GGA') then
    gga = .true.
  else
    write(iowrite,*) 'atomxc: unknown functional ', functl

    stop

  endif

! allocate aux

  if (gga) then
    allocate(aux(nr))
  endif

! initialize output

  ex = ZERO
  ec = ZERO
  dx = ZERO
  dc = ZERO

  do is = 1,nspin
    do ir = 1,nr
      vxc(ir,is) = ZERO
    enddo
  enddo

! get number pi

  pi = 4 * atan(ONE)

! loop on mesh points

  do ir = 1,nr

!   find interval of neighbour points to calculate derivatives

    in1 = max(  1, ir-NN ) - ir
    in2 = min( nr, ir+NN ) - ir

!   find weights of numerical derivation from lagrange
!   interpolation formula

    do in = in1,in2
      if (in  ==  0) then
        dgdm(in) = 0
        do jn = in1,in2
          if (jn /= 0) dgdm(in) = dgdm(in) + ONE / (0 - jn)
        enddo
      else
        f1 = 1
        f2 = 1
        do jn = in1,in2
          if (jn /= in .and. jn /= 0) f1 = f1 * (0  - jn)
          if (jn /= in) f2 = f2 * (in - jn)
        enddo
        dgdm(in) = f1 / f2
      endif
    enddo

!   find dr/dmesh

    drdm = 0
    do in = in1,in2
      drdm = drdm + rmesh(ir+in) * dgdm(in)
    enddo

!   find differential of volume. use trapezoidal integration rule

    dvol = 4 * pi * rmesh(ir)*rmesh(ir) * drdm

!   dvmin is a small number added to avoid a division by zero

    dvol = dvol + dvmin
    if (ir == 1 .or. ir == nr) dvol = dvol / 2
    if (gga) aux(ir) = dvol

!   find the weights for the derivative d(gradf(i))/d(f(j)), of
!   the gradient at point i with respect to the value at point j

    if (gga) then
      do in = in1,in2
        dgidfj(in) = dgdm(in) / drdm
      enddo
    endif

!   find density and gradient of density at this point

    do is = 1,nspin
      d(is) = dens(ir,is)
    enddo
    if (gga) then
      do is = 1,nspin
        gd(1,is) = ZERO
        gd(2,is) = ZERO
        gd(3,is) = ZERO
        do in = in1,in2
          gd(3,is) = gd(3,is) + dgidfj(in) * dens(ir+in,is)
        enddo
      enddo
    endif

!   find exchange and correlation energy densities and their
!   derivatives with respect to density and density gradient

    if (gga) then
      call atom_xc_gga( author, irel, nspin, d, gd,                      &
                 epsx, epsc, dexdd, decdd, dexdgd, decdgd )
    else
      call atom_xc_lda( author, irel, nspin, d, epsx, epsc, dexdd, decdd)
    endif

!   add contributions to exchange-correlation energy and its
!   derivatives with respect to density at all points

    do is = 1,nspin
      ex = ex + dvol * d(is) * epsx
      ec = ec + dvol * d(is) * epsc
      dx = dx + dvol * d(is) * (epsx - dexdd(is))
      dc = dc + dvol * d(is) * (epsc - decdd(is))
      if (gga) then
        vxc(ir,is) = vxc(ir,is) + dvol * ( dexdd(is) + decdd(is) )
        do in = in1,in2
          dx= dx - dvol * dens(ir+in,is) * dexdgd(3,is) * dgidfj(in)
          dc= dc - dvol * dens(ir+in,is) * decdgd(3,is) * dgidfj(in)
          vxc(ir+in,is) = vxc(ir+in,is) + dvol *                         &
                (dexdgd(3,is) + decdgd(3,is)) * dgidfj(in)
        enddo
      else
        vxc(ir,is) = dexdd(is) + decdd(is)
      endif
    enddo

  enddo

! divide by volume element to obtain the potential (per electron)

  if (gga) then
    do is = 1,nspin
      do ir = 1,nr
        dvol = aux(ir)
        vxc(ir,is) = vxc(ir,is) / dvol
      enddo
    enddo
  endif

! divide by energy unit

  ex = ex / EUNIT
  ec = ec / EUNIT
  dx = dx / EUNIT
  dc = dc / EUNIT
  do is = 1,nspin
    do ir = 1,nr
      vxc(ir,is) = vxc(ir,is) / EUNIT
    enddo
  enddo

  if (gga) then
    deallocate(aux)
  endif

  end subroutine atom_xc
