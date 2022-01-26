!>  Finds the exchange and correlation energies at a point, and their
!>  derivatives with respect to density and density gradient, in the
!>  Generalized Gradient Correction approximation.
!>  Lengths in Bohr, energies in Hartrees
!>
!>  \author       L.C.Balbas, J.M.Soler, Jose Luis Martins
!>  \version      6.013
!>  \date         December 96, 22 June 2021
!>  \copyright    GNU Public License v2

subroutine atom_xc_gga( author, irel, nspin, dens, gdens,                &
                   ex, ec, dexdd, decdd, dexdgd, decdgd )

! Written by L.C.Balbas and J.M.Soler, Dec'96. Version 0.5.
! density should always be positive   jlm  August 2001

! converted to f90, July 2018
! cleanup and new interface, July 2019. JLM
! jlm  version 6.00

  implicit          none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter          ::  iowrite = 6                            !  default tape for writing

! input

  character(len=*), intent(in)      ::  author                           !<  parametrization desired
  integer, intent(in)               ::  irel                             !<  relativistic exchange? (0=>no, 1=>yes)

  integer, intent(in)               ::  nspin                            !<  nspin=1 => unpolarized; nspin=2 => polarized

  real(REAL64), intent(in)          ::  dens(nspin)                      !<  charge density
  real(REAL64), intent(in)          ::  gdens(3,nspin)                   !<  charge density

! output

  real(REAL64), intent(out)         ::  ex                               !<  exchange energy density
  real(REAL64), intent(out)         ::  ec                               !<  correlation energy density
  real(REAL64), intent(out)         ::  dexdd(nspin)                     !<  d rho*ex / d dens(spin)
  real(REAL64), intent(out)         ::  decdd(nspin)                     !<  d rho*ec / d dens(spin)
  real(REAL64), intent(out)         ::  dexdgd(3,nspin)                  !<  d rho*ex / d dgens(i,spin)
  real(REAL64), intent(out)         ::  decdgd(3,nspin)                  !<  d rho*ec / d dgens(i,spin)

! local

  real(REAL64)    ::  rho(nspin)

! constants

  real(REAL64), parameter   ::  ZERO = 0.0_REAL64

! counters

  integer         ::  i

! jlm
  do i = 1,nspin
    rho(i) = dens(i)
    if(rho(i) < ZERO) rho(i) = ZERO
  enddo
! jlm

  if (author == 'pbe' .or. author == 'PBE') then
    call atom_xc_pbe( irel, nspin, rho, gdens,                           &
               ex, ec, dexdd, decdd, dexdgd, decdgd )
  else
    write(iowrite,*) 'xc_gga: unknown author ', author

    stop

  endif

  return

end subroutine atom_xc_gga
