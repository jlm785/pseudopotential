!>  Calculates the charge inside the core radius: int_0^r_c |psi|^2 dr
!>
!>  \author       Sverre Froyen,  Norm Troullier, Jose Luis Martins
!>  \version      6.013
!>  \date         1980s and 1990s, 30 June 2021
!>  \copyright    GNU Public License v2

subroutine atom_psd_psi_charge_nrl(jrc, drdi, rpsi, cdrc,                &
     mxdnr)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension  of the number of radial points

  integer, intent(in)               ::  jrc                              !<  r(jrc) = r_cut

  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  real(REAL64), intent(in)          ::  rpsi(mxdnr)                      !<  r*psi

! output

  real(REAL64), intent(out)         ::  cdrc                             !<  int_0^r_cut rho(r) dr

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64

! local variables

  integer            ::  ll

! counters

  integer  ::  k

  ll = 2

  cdrc = -rpsi(jrc)*rpsi(jrc) * drdi(jrc)
  if (jrc /= 2*(jrc/2)) then
    do k = jrc,1,-1
      cdrc = cdrc +  ll * rpsi(k)*rpsi(k) * drdi(k)
      ll = 6 - ll
    enddo
  else
    do k = jrc,4,-1
      cdrc = cdrc +  ll * rpsi(k)*rpsi(k) * drdi(k)
      ll = 6 - ll
    enddo
    cdrc = cdrc - rpsi(4)*rpsi(4) * drdi(4)
    cdrc = cdrc + 9*(  rpsi(1)*rpsi(1) * drdi(1) +                       &
                     3*rpsi(2)*rpsi(2) * drdi(2) +                       &
                     3*rpsi(3)*rpsi(3) * drdi(3) +                       &
                       rpsi(4)*rpsi(4) * drdi(4)) / (8*ONE)
  endif
  cdrc = cdrc/3

  return

end subroutine atom_psd_psi_charge_nrl
