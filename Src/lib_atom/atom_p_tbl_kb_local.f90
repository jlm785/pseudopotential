!>  Function determines the default KB local potential channel
!>  for a TM pseudopotential
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.7
!>  \date         July 2021. 21 December 2021.
!>  \copyright    GNU Public License v2

subroutine atom_p_tbl_kb_local(name, llocal, jhard)

! jhard replace nval2. 21 December 2021. JLM

  implicit none

! input

  character(len=2), intent(in)      ::  name                             !<  chemical symbol of the atom
  integer, intent(in)               ::  jhard                            !<  flag for accuracy/speed compromise (rarely used here)

! output

  integer, intent(out)              ::  llocal                           !<  angular momentum for local potential (negative: maximum of l-dependent)

  if (name == 'H ' .or. name == ' H') then
    llocal = 0
  elseif (name == 'He') then
    llocal = 0
  elseif (name == 'Li') then
    llocal = -1
  elseif (name == 'Be') then
    llocal = -1
  elseif (name == 'B ' .or. name == ' B') then
    llocal = -1
  elseif (name == 'C ' .or. name == ' C') then
    llocal = -1
  elseif (name == 'N ' .or. name == ' N') then
    llocal = -1
  elseif (name == 'O ' .or. name == ' O') then
    llocal = -1
  elseif (name == 'F ' .or. name == ' F') then
    llocal = -1
  elseif (name == 'Ne') then
    llocal = -1
  elseif (name == 'Na') then
    llocal = -1
  elseif (name == 'Mg') then
    llocal = -1
  elseif (name == 'Al') then
    llocal = -1
  elseif (name == 'Si') then
    llocal = -1
  elseif (name == 'P ' .or. name == ' P') then
    llocal = -1
  elseif (name == 'S ' .or. name == ' S') then
    llocal = -1
  elseif (name == 'Cl') then
    llocal = -1
  elseif (name == 'Ar') then
    llocal = -1
  elseif (name == 'K ' .or. name == ' K') then
    llocal = -1
  elseif (name == 'Ca') then
    llocal = -1
  elseif (name == 'Sc') then
    llocal = -1
  elseif (name == 'Ti') then
    llocal = -1
  elseif (name == 'V ' .or. name == ' V') then
    llocal = -1
  elseif (name == 'Cr') then
    llocal = -1
  elseif (name == 'Mn') then
    llocal = -1
  elseif (name == 'Fe') then
    llocal = -1
  elseif (name == 'Co') then
    llocal = -1
  elseif (name == 'Ni') then
    llocal = -1
  elseif (name == 'Cu') then
    llocal = -1
  elseif (name == 'Zn') then
    llocal = -1
  elseif (name == 'Ga') then
    llocal = -1
  elseif (name == 'Ge') then
    llocal = -1
  elseif (name == 'As') then
    llocal = -1
  elseif (name == 'Se') then
    llocal = -1
  elseif (name == 'Br') then
    llocal = -1
  elseif (name == 'Kr') then
    llocal = -1
  elseif (name == 'Rb') then
    llocal = -1
  elseif (name == 'Sr') then
    llocal = -1
  elseif (name == 'Y ' .or. name == ' Y') then
    llocal = -1
  elseif (name == 'Zr') then
    llocal = -1
  elseif (name == 'Nb') then
    llocal = -1
  elseif (name == 'Mo') then
    llocal = -1
  elseif (name == 'Tc') then
    llocal = -1
  elseif (name == 'Ru') then
    llocal = -1
  elseif (name == 'Rh') then
    llocal = -1
  elseif (name == 'Pd') then
    llocal = -1
  elseif (name == 'Ag') then
    llocal = -1
  elseif (name == 'Cd') then
    llocal = -1
  elseif (name == 'In') then
    llocal = -1
  elseif (name == 'Sn') then
    llocal = -1
  elseif (name == 'Sb') then
    llocal = -1
  elseif (name == 'Te') then
    llocal = -1
  elseif (name == 'I ' .or. name == ' I') then
    llocal = -1
  elseif (name == 'Xe') then
    llocal = -1
  elseif (name == 'Cs') then
    llocal = -1
  elseif (name == 'Ba') then
    llocal = -1
  elseif (name == 'La') then
    llocal = -1
  elseif (name == 'Ce') then
    llocal = -1
  elseif (name == 'Pr') then
    llocal = -1
  elseif (name == 'Nd') then
    llocal = -1
  elseif (name == 'Pm') then
    llocal = -1
  elseif (name == 'Sm') then
    llocal = -1
  elseif (name == 'Eu') then
    llocal = -1
  elseif (name == 'Gd') then
    llocal = -1
  elseif (name == 'Tb') then
    llocal = -1
  elseif (name == 'Dy') then
    llocal = -1
  elseif (name == 'Ho') then
    llocal = -1
  elseif (name == 'Er') then
    llocal = -1
  elseif (name == 'Tm') then
    llocal = -1
  elseif (name == 'Yb') then
    llocal = -1
  elseif (name == 'Lu') then
    llocal = -1
  elseif (name == 'Hf') then
    llocal = -1
  elseif (name == 'Ta') then
    llocal = -1
  elseif (name == 'W ' .or. name == ' W') then
    llocal = -1
  elseif (name == 'Re') then
    llocal = -1
  elseif (name == 'Os') then
    llocal = -1
  elseif (name == 'Ir') then
    llocal = -1
  elseif (name == 'Pt') then
    llocal = -1
  elseif (name == 'Au') then
    llocal = -1
  elseif (name == 'Hg') then
    llocal = -1
  elseif (name == 'Tl') then
    llocal = -1
  elseif (name == 'Pb') then
    llocal = -1
  elseif (name == 'Bi') then
    llocal = -1
  elseif (name == 'Po') then
    llocal = -1
  elseif (name == 'At') then
    llocal = -1
  elseif (name == 'Rn') then
    llocal = -1
  elseif (name == 'Fr') then
    llocal = -1
  elseif (name == 'Ra') then
    llocal = -1
  elseif (name == 'Ac') then
    llocal = -1
  elseif (name == 'Th') then
    llocal = -1
  elseif (name == 'Pa') then
    llocal = -1
  elseif (name == ' U' .or. name == 'U ') then
    llocal = -1
  elseif (name == 'Np') then
    llocal = -1
  elseif (name == 'Pu') then
    llocal = -1
  elseif (name == 'Am') then
    llocal = -1
  elseif (name == 'Cm') then
    llocal = -1
  elseif (name == 'Bk') then
    llocal = -1
  elseif (name == 'Cf') then
    llocal = -1
  elseif (name == 'Es') then
    llocal = -1
  elseif (name == 'Fm') then
    llocal = -1
  elseif (name == 'Md') then
    llocal = -1
  elseif (name == 'No') then
    llocal = -1
  elseif (name == 'Lr') then
    llocal = -1
  elseif (name == 'Rf') then
    llocal = -1
  elseif (name == 'Db') then
    llocal = -1
  elseif (name == 'Sg') then
    llocal = -1
  elseif (name == 'Bh') then
    llocal = -1
  elseif (name == 'Hs') then
    llocal = -1
  elseif (name == 'Mt') then
    llocal = -1
  elseif (name == 'Ds') then
    llocal = -1
  elseif (name == 'Rg') then
    llocal = -1
  elseif (name == 'Cn') then
    llocal = -1
  elseif (name == 'Nh') then
    llocal = -1
  elseif (name == 'Fl') then
    llocal = -1
  elseif (name == 'Mc') then
    llocal = -1
  elseif (name == 'Lv') then
    llocal = -1
  elseif (name == 'Ts') then
    llocal = -1
  elseif (name == 'Og') then
    llocal = -1
  else
    write(6,*)
    write(6,*)
    write(6,'("  Element ",a2," unknown")') name
    write(6,'("  Using llocal = -2")')
    write(6,*)
    llocal = -2
  endif

  return

end subroutine atom_p_tbl_kb_local
