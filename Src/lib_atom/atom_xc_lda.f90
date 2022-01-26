!>  Finds the exchange and correlation energies and potentials, in the
!>  local (spin) density approximation.
!>  Lengths in bohr, energies in hartrees.
!>
!>  \author       L.C.Balbas and J.M.Soler
!>  \version      6.013
!>  \date         December 96, 22 June 2021
!>  \copyright    GNU Public License v2

subroutine atom_xc_lda( author, irel, nspin, dens, epsx, epsc, vx, vc )

! Written by L.C.Balbas and J.M.Soler, dec'96. version 0.5.
! density should always be positive.   jlm  august  2001


! converted to f90, May 2018
! cleanup and new interface, July 2019. JLM
! jlm  version 6.00


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter          ::  iowrite = 6                            !  default tape for writing

! input

  character(len=*), intent(in)      ::  author                           !<  parametrization desired
  integer, intent(in)               ::  irel                             !<  relativistic exchange? (0 => no, 1 => yes)

  integer, intent(in)               ::  nspin                            !<  nspin=1 => unpolarized; nspin=2 => polarized

  real(REAL64), intent(in)          ::  dens(nspin)                      !<  charge density

! output

  real(REAL64), intent(out)         ::  epsx                             !<  exchnange  energy densities
  real(REAL64), intent(out)         ::  epsc                             !<  correlation energy densities
  real(REAL64), intent(out)         ::  vx(nspin)                        !<  exchnange potential
  real(REAL64), intent(out)         ::  vc(nspin)                        !<  correlation potential

! local

  real(REAL64)    ::  rho(nspin)

! constants

  real(REAL64), parameter   ::  ZERO = 0.0_REAL64

! counters

  integer         ::  i

! jlm
  do i = 1,nspin
    rho(i) = dens(i)
    if(dens(i) < ZERO) rho(i) = ZERO
  enddo
! jlm

  if ( author == 'ca' .or. author == 'CA' .or.                           &
      author == 'pz' .or. author == 'PZ') then

    call atom_xc_x( irel, nspin, rho, epsx, vx )
    call atom_xc_pzc( nspin, rho, epsc, vc )

  elseif ( author == 'pw92' .or. author == 'PW92' ) then

    call atom_xc_x( irel, nspin, rho, epsx, vx )
    call atom_xc_pw92c( nspin, rho, epsc, vc )

  else
    write(iowrite,*) '  xc_lda: unknown author ', author

    stop

  endif

  return
  end subroutine atom_xc_lda

