!>  Finds the core radii and sets xc partial core correction
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.6
!>  \date         27 August 2021. 2 November 2021.
!>  \copyright    GNU Public License v2

subroutine atom_psd_core_radii(rc, cfac, rcfac, rc_tab, rz, rx, lpmax,   &
         iowrite, ioae, ioread, filein, fileae)

! mxdlc in atom_psd_max_zero.  2 November 2021. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter                ::  mxdlc = 4                        !  hard coded maximum valence angular momentum

! input

  integer, intent(in)               ::  iowrite                          !<  default output tape

  integer, intent(in)               ::  ioae                             !<  default tape for writing
  character(len=*), intent(in)      ::  fileae                           !<  name of default tape for all-electron results

  integer, intent(in)               ::  ioread                           !<  default tape for reading
  character(len=*), intent(in)      ::  filein                           !<  name of default tape for reading

! output

  real(REAL64), intent(out)         ::  rc(0:mxdlc)                      !<  core radius r_c(l)
  real(REAL64), intent(out)         ::  cfac                             !<  criteria for pseudo-core charge
  real(REAL64), intent(out)         ::  rcfac                            !<  pseudo-core radius

  real(REAL64), intent(out)         ::  rc_tab(0:mxdlc)                  !<  tabulated core radii
  real(REAL64), intent(out)         ::  rz(0:mxdlc), rx(0:mxdlc)         !<  avearge zero and extrma for l

  integer, intent(out)              ::  lpmax                            !<  maximum value of l in pseudopotentials

! local variables

  character(len=2)        ::  nameat                                     !  chemical symbol of the atom
  logical                 ::  lex, lres
  character(len=30)       ::  status                                     !  quality of psedopotential

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64

! counter

  integer     ::  l

  call atom_psd_print_info(iowrite, ioae, fileae, nameat)

! tabulated and default values

  call atom_p_tbl_psd_tm2(nameat, rc_tab, status)

! reads core radii fromatom.dat

  inquire(file=filein, exist = lex)

  lres = .FALSE.
  if(lex) then
    call atom_psd_read_atom(ioread, filein, nameat, rc, cfac, rcfac, lres)
  endif

  if(.not. lres .or. .not. lex) then
    cfac = ONE
    rcfac = ZERO

    do l = 0,mxdlc
      rc(l) = rc_tab(l)
    enddo

    write(iowrite,*)
    write(iowrite,*) '  Core radii obtained from a table with reliability identified as: '
    write(iowrite,*) '  ',status
    write(iowrite,*)

  endif

  call atom_psd_max_zero(rz, rx, mxdlc, iowrite, ioae, fileae)

  do l = 0,mxdlc
    if(rc(l) > 0.0001) lpmax = l
  enddo

  return

end subroutine atom_psd_core_radii
