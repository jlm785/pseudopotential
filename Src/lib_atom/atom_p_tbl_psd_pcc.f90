!>  Indicates if partial core correction is recommended
!>  for troullier Martins pseudopotential generation.
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.7
!>  \date         22 June 2021, 21 December 2021.
!>  \copyright    GNU Public License v2

subroutine atom_p_tbl_psd_pcc(name, pcc, jhard)

!  ***********************************************************
!  *                                                         *
!  *   All elements from H to Og are included.               *
!  *                                                         *
!  *  Version dated May 1, 1991                              *
!  *  njtj                                                   *
!  *                                                         *
!  *  Modified for integer. July 2013                        *
!  *                                                         *
!  *  Added new elements 2 January 2018.  JLM                *
!  ***********************************************************

! cleanup and new interface, July 2019. JLM
! jhard 21 December 2021. JLM

  implicit none

! input

  character(len=2), intent(in)      ::  name                             !<  chemical symbol of the atom
  integer, intent(in)               ::  jhard                            !<  flag for accuracy/speed compromise (-1 no pcc)

! output

  character(len=2), intent(out)     ::  pcc                              !<  no pcc: pg; pcc: pe.

  if (name == 'H ' .or. name == ' H') then
    pcc = 'pg'
  elseif (name == 'He') then
    pcc = 'pg'
  elseif (name == 'Li') then
    pcc = 'pe'
  elseif (name == 'Be') then
    pcc = 'pe'
  elseif (name == 'B ' .or. name == ' B') then
    pcc = 'pe'
  elseif (name == 'C ' .or. name == ' C') then
    pcc = 'pg'
  elseif (name == 'N ' .or. name == ' N') then
    pcc = 'pg'
  elseif (name == 'O ' .or. name == ' O') then
    pcc = 'pg'
  elseif (name == 'F ' .or. name == ' F') then
    pcc = 'pg'
  elseif (name == 'Ne') then
    pcc = 'pg'
  elseif (name == 'Na') then
    pcc = 'pe'
  elseif (name == 'Mg') then
    pcc = 'pe'
  elseif (name == 'Al') then
    pcc = 'pg'
  elseif (name == 'Si') then
    pcc = 'pg'
  elseif (name == 'P ' .or. name == ' P') then
    pcc = 'pg'
  elseif (name == 'S ' .or. name == ' S') then
    pcc = 'pg'
  elseif (name == 'Cl') then
    pcc = 'pg'
  elseif (name == 'Ar') then
    pcc = 'pg'
  elseif (name == 'K ' .or. name == ' K') then
    pcc = 'pe'
  elseif (name == 'Ca') then
    pcc = 'pe'
  elseif (name == 'Sc') then
    pcc = 'pg'
  elseif (name == 'Ti') then
    pcc = 'pg'
  elseif (name == 'V ' .or. name == ' V') then
    pcc = 'pg'
  elseif (name == 'Cr') then
    pcc = 'pg'
  elseif (name == 'Mn') then
    pcc = 'pg'
  elseif (name == 'Fe') then
    pcc = 'pg'
  elseif (name == 'Co') then
    pcc = 'pg'
  elseif (name == 'Ni') then
    pcc = 'pg'
  elseif (name == 'Cu') then
    pcc = 'pg'
  elseif (name == 'Zn') then
    pcc = 'pg'
  elseif (name == 'Ga') then
    pcc = 'pe'
  elseif (name == 'Ge') then
    pcc = 'pe'
  elseif (name == 'As') then
    pcc = 'pg'
  elseif (name == 'Se') then
    pcc = 'pg'
  elseif (name == 'Br') then
    pcc = 'pg'
  elseif (name == 'Kr') then
    pcc = 'pg'
  elseif (name == 'Rb') then
    pcc = 'pe'
  elseif (name == 'Sr') then
    pcc = 'pe'
  elseif (name == 'Y ' .or. name == ' Y') then
    pcc = 'pg'
  elseif (name == 'Zr') then
    pcc = 'pg'
  elseif (name == 'Nb') then
    pcc = 'pg'
  elseif (name == 'Mo') then
    pcc = 'pg'
  elseif (name == 'Tc') then
    pcc = 'pg'
  elseif (name == 'Ru') then
    pcc = 'pg'
  elseif (name == 'Rh') then
    pcc = 'pg'
  elseif (name == 'Pd') then
    pcc = 'pg'
  elseif (name == 'Ag') then
    pcc = 'pg'
  elseif (name == 'Cd') then
    pcc = 'pg'
  elseif (name == 'In') then
    pcc = 'pg'
  elseif (name == 'Sn') then
    pcc = 'pg'
  elseif (name == 'Sb') then
    pcc = 'pg'
  elseif (name == 'Te') then
    pcc = 'pg'
  elseif (name == 'I ' .or. name == ' I') then
    pcc = 'pg'
  elseif (name == 'Xe') then
    pcc = 'pg'
  elseif (name == 'Cs') then
    pcc = 'pe'
  elseif (name == 'Ba') then
    pcc = 'pe'
  elseif (name == 'La') then
    pcc = 'pe'
  elseif (name == 'Ce') then
    pcc = 'pe'
  elseif (name == 'Pr') then
    pcc = 'pe'
  elseif (name == 'Nd') then
    pcc = 'pe'
  elseif (name == 'Pm') then
    pcc = 'pe'
  elseif (name == 'Sm') then
    pcc = 'pe'
  elseif (name == 'Eu') then
    pcc = 'pe'
  elseif (name == 'Gd') then
    pcc = 'pe'
  elseif (name == 'Tb') then
    pcc = 'pe'
  elseif (name == 'Dy') then
    pcc = 'pe'
  elseif (name == 'Ho') then
    pcc = 'pe'
  elseif (name == 'Er') then
    pcc = 'pe'
  elseif (name == 'Tm') then
    pcc = 'pe'
  elseif (name == 'Yb') then
    pcc = 'pe'
  elseif (name == 'Lu') then
    pcc = 'pe'
  elseif (name == 'Hf') then
    pcc = 'pg'
  elseif (name == 'Ta') then
    pcc = 'pg'
  elseif (name == 'W ' .or. name == ' W') then
    pcc = 'pg'
  elseif (name == 'Re') then
    pcc = 'pg'
  elseif (name == 'Os') then
    pcc = 'pg'
  elseif (name == 'Ir') then
    pcc = 'pg'
  elseif (name == 'Pt') then
    pcc = 'pg'
  elseif (name == 'Au') then
    pcc = 'pg'
  elseif (name == 'Hg') then
    pcc = 'pg'
  elseif (name == 'Tl') then
    pcc = 'pg'
  elseif (name == 'Pb') then
    pcc = 'pg'
  elseif (name == 'Bi') then
    pcc = 'pg'
  elseif (name == 'Po') then
    pcc = 'pg'
  elseif (name == 'At') then
    pcc = 'pg'
  elseif (name == 'Rn') then
    pcc = 'pg'
  elseif (name == 'Fr') then
    pcc = 'pg'
  elseif (name == 'Ra') then
    pcc = 'pg'
  elseif (name == 'Ac') then
    pcc = 'pg'
  elseif (name == 'Th') then
    pcc = 'pg'
  elseif (name == 'Pa') then
    pcc = 'pg'
  elseif (name == ' U' .or. name == 'U ') then
    pcc = 'pg'
  elseif (name == 'Np') then
    pcc = 'pg'
  elseif (name == 'Pu') then
    pcc = 'pg'
  elseif (name == 'Am') then
    pcc = 'pg'
  elseif (name == 'Cm') then
    pcc = 'pg'
  elseif (name == 'Bk') then
    pcc = 'pg'
  elseif (name == 'Cf') then
    pcc = 'pg'
  elseif (name == 'Es') then
    pcc = 'pg'
  elseif (name == 'Fm') then
    pcc = 'pg'
  elseif (name == 'Md') then
    pcc = 'pg'
  elseif (name == 'No') then
    pcc = 'pg'
  elseif (name == 'Lr') then
    pcc = 'pg'
  elseif (name == 'Rf') then
    pcc = 'pg'
  elseif (name == 'Db') then
    pcc = 'pg'
  elseif (name == 'Sg') then
    pcc = 'pg'
  elseif (name == 'Bh') then
    pcc = 'pg'
  elseif (name == 'Hs') then
    pcc = 'pg'
  elseif (name == 'Mt') then
    pcc = 'pg'
  elseif (name == 'Ds') then
    pcc = 'pg'
  elseif (name == 'Rg') then
    pcc = 'pg'
  elseif (name == 'Cn') then
    pcc = 'pg'
  elseif (name == 'Nh') then
    pcc = 'pg'
  elseif (name == 'Fl') then
    pcc = 'pg'
  elseif (name == 'Mc') then
    pcc = 'pg'
  elseif (name == 'Lv') then
    pcc = 'pg'
  elseif (name == 'Ts') then
    pcc = 'pg'
  elseif (name == 'Og') then
    pcc = 'pg'
  else
    write(6,*)
    write(6,*)
    write(6,'("  Element ",a2," unknown")') name
    write(6,'("  Using charge 200")')
    write(6,*)
    pcc = 'pg'
  endif

  if(jhard == -1) pcc = 'pg'

  return

end subroutine atom_p_tbl_psd_pcc
