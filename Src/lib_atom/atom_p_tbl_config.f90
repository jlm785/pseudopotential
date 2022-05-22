!>  Function determines the canonical number of core orbitals of an element
!>  and the atomic ground state valence configuration.
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.7
!>  \date         12 July 2021, 21 December 2021
!>  \copyright    GNU Public License v2

subroutine atom_p_tbl_config(name, ncore, nval, no, lo, zo, jhard)

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
! jhard replace nval2. 21 December 2021. JLM

  implicit none

  integer, parameter          ::  REAL64 = selected_real_kind(12)
  integer, parameter          ::  lc = 4

! input

  character(len=2), intent(in)      ::  name                             !<  chemical symbol of the atom
  integer, intent(in)               ::  jhard                            !<  flag for accuracy/speed compromise

! output

  integer, intent(out)              ::  ncore                            !<  canonical number of core orbitals
  integer, intent(out)              ::  nval                             !<  canonical number of interesting valence orbitals
  integer, intent(out)              ::  no(lc+1), lo(lc+1)               !<  configuration
  real(REAL64), intent(out)         ::  zo(lc+1)                         !<  occupation

  real(REAL64), parameter    ::  ONE = 1.0_REAL64

  if (name == 'H ' .or. name == ' H') then
    ncore = 0
    nval = 1
    no(1) = 1
    lo(1) = 0
    zo(1) = 1*ONE
  elseif (name == 'He') then
    ncore = 0
    nval = 1
    no(1) = 1
    lo(1) = 0
    zo(1) = 2*ONE
  elseif (name == 'Li') then
    ncore = 1
    nval = 2
    no(1) = 2
    lo(1) = 0
    zo(1) = 1*ONE
    no(2) = 2
    lo(2) = 1
    zo(2) = 0*ONE
    if(jhard == 1) then
      nval = 3
      no(3) = 3
      lo(3) = 2
      zo(3) = 0*ONE
    endif
 elseif (name == 'Be') then
    ncore = 1
    nval = 2
    no(1) = 2
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 2
    lo(2) = 1
    zo(2) = 0*ONE
    if(jhard == 1) then
      nval = 3
      no(3) = 3
      lo(3) = 2
      zo(3) = 0*ONE
    endif
  elseif (name == 'B ' .or. name == ' B') then
    ncore = 1
    nval = 2
    no(1) = 2
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 2
    lo(2) = 1
    zo(2) = 1*ONE
    if(jhard == 1) then
      nval = 3
      no(3) = 3
      lo(3) = 2
      zo(3) = 0*ONE
    endif
  elseif (name == 'C ' .or. name == ' C') then
    ncore = 1
    nval = 2
    no(1) = 2
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 2
    lo(2) = 1
    zo(2) = 2*ONE
    if(jhard == 1) then
      nval = 3
      no(3) = 3
      lo(3) = 2
      zo(3) = 0*ONE
    endif
  elseif (name == 'N ' .or. name == ' N') then
    ncore = 1
    nval = 2
    no(1) = 2
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 2
    lo(2) = 1
    zo(2) = 3*ONE
    if(jhard == 1) then
      nval = 3
      no(3) = 3
      lo(3) = 2
      zo(3) = 0*ONE
    endif
  elseif (name == 'O ' .or. name == ' O') then
    ncore = 1
    nval = 2
    no(1) = 2
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 2
    lo(2) = 1
    zo(2) = 4*ONE
    if(jhard == 1) then
      nval = 3
      no(3) = 3
      lo(3) = 2
      zo(3) = 0*ONE
    endif
  elseif (name == 'F ' .or. name == ' F') then
    ncore = 1
    nval = 2
    no(1) = 2
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 2
    lo(2) = 1
    zo(2) = 5*ONE
    if(jhard == 1) then
      nval = 3
      no(3) = 3
      lo(3) = 2
      zo(3) = 0*ONE
    endif
  elseif (name == 'Ne') then
    ncore = 1
    nval = 2
    no(1) = 2
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 2
    lo(2) = 1
    zo(2) = 6*ONE
  elseif (name == 'Na') then
    ncore = 3
    nval = 2
    no(1) = 3
    lo(1) = 0
    zo(1) = 1*ONE
    no(2) = 3
    lo(2) = 1
    zo(2) = 0*ONE
    if(jhard == 1) then
      nval = 3
      no(3) = 3
      lo(3) = 2
      zo(3) = 0*ONE
    endif
  elseif (name == 'Mg') then
    ncore = 3
    nval = 3
    no(1) = 3
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 3
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 3
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Al') then
    ncore = 3
    nval = 3
    no(1) = 3
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 3
    lo(2) = 1
    zo(2) = 1*ONE
    no(3) = 3
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Si') then
    ncore = 3
    nval = 3
    no(1) = 3
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 3
    lo(2) = 1
    zo(2) = 2*ONE
    no(3) = 3
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'P ' .or. name == ' P') then
    ncore = 3
    nval = 3
    no(1) = 3
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 3
    lo(2) = 1
    zo(2) = 3*ONE
    no(3) = 3
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'S ' .or. name == ' S') then
    ncore = 3
    nval = 3
    no(1) = 3
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 3
    lo(2) = 1
    zo(2) = 4*ONE
    no(3) = 3
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Cl') then
    ncore = 3
    nval = 3
    no(1) = 3
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 3
    lo(2) = 1
    zo(2) = 5*ONE
    no(3) = 3
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Ar') then
    ncore = 3
    nval = 3
    no(1) = 3
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 3
    lo(2) = 1
    zo(2) = 6*ONE
    no(3) = 3
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'K ' .or. name == ' K') then
    ncore = 5
    nval = 3
    no(1) = 4
    lo(1) = 0
    zo(1) = 1*ONE
    no(2) = 4
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 3
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Ca') then
    ncore = 5
    nval = 3
    no(1) = 4
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 4
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 3
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Sc') then
    ncore = 5
    nval = 3
    no(1) = 4
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 4
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 3
    lo(3) = 2
    zo(3) = 1*ONE
  elseif (name == 'Ti') then
    ncore = 5
    nval = 3
    no(1) = 4
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 4
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 3
    lo(3) = 2
    zo(3) = 2*ONE
  elseif (name == 'V ' .or. name == ' V') then
    ncore = 5
    nval = 3
    no(1) = 4
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 4
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 3
    lo(3) = 2
    zo(3) = 3*ONE
  elseif (name == 'Cr') then
    ncore = 5
    nval = 3
    no(1) = 4
    lo(1) = 0
    zo(1) = 1*ONE
    no(2) = 4
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 3
    lo(3) = 2
    zo(3) = 5*ONE
  elseif (name == 'Mn') then
    ncore = 5
    nval = 3
    no(1) = 4
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 4
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 3
    lo(3) = 2
    zo(3) = 5*ONE
  elseif (name == 'Fe') then
    ncore = 5
    nval = 3
    no(1) = 4
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 4
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 3
    lo(3) = 2
    zo(3) = 6*ONE
  elseif (name == 'Co') then
    ncore = 5
    nval = 3
    no(1) = 4
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 4
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 3
    lo(3) = 2
    zo(3) = 7*ONE
  elseif (name == 'Ni') then
    ncore = 5
    nval = 3
    no(1) = 4
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 4
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 3
    lo(3) = 2
    zo(3) = 8*ONE
  elseif (name == 'Cu') then
    ncore = 5
    nval = 3
    no(1) = 4
    lo(1) = 0
    zo(1) = 1*ONE
    no(2) = 4
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 3
    lo(3) = 2
    zo(3) = 10*ONE
  elseif (name == 'Zn') then
    ncore = 5
    nval = 3
    no(1) = 4
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 4
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 3
    lo(3) = 2
    zo(3) = 10*ONE
  elseif (name == 'Ga') then
    ncore = 6
    nval = 3
    no(1) = 4
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 4
    lo(2) = 1
    zo(2) = 1*ONE
    no(3) = 4
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Ge') then
    ncore = 6
    nval = 3
    no(1) = 4
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 4
    lo(2) = 1
    zo(2) = 2*ONE
    no(3) = 4
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'As') then
    ncore = 6
    nval = 3
    no(1) = 4
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 4
    lo(2) = 1
    zo(2) = 3*ONE
    no(3) = 4
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Se') then
    ncore = 6
    nval = 3
    no(1) = 4
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 4
    lo(2) = 1
    zo(2) = 4*ONE
    no(3) = 4
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Br') then
    ncore = 6
    nval = 3
    no(1) = 4
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 4
    lo(2) = 1
    zo(2) = 5*ONE
    no(3) = 4
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Kr') then
    ncore = 6
    nval = 3
    no(1) = 4
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 4
    lo(2) = 1
    zo(2) = 6*ONE
    no(3) = 4
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Rb') then
    ncore = 8
    nval = 3
    no(1) = 5
    lo(1) = 0
    zo(1) = 1*ONE
    no(2) = 5
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 4
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Sr') then
    ncore = 8
    nval = 3
    no(1) = 5
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 5
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 4
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Y ' .or. name == ' Y') then
    ncore = 8
    nval = 3
    no(1) = 5
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 5
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 4
    lo(3) = 2
    zo(3) = 1*ONE
  elseif (name == 'Zr') then
    ncore = 8
    nval = 3
    no(1) = 5
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 5
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 4
    lo(3) = 2
    zo(3) = 2*ONE
  elseif (name == 'Nb') then
    ncore = 8
    nval = 3
    no(1) = 5
    lo(1) = 0
    zo(1) = 1*ONE
    no(2) = 5
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 4
    lo(3) = 2
    zo(3) = 4*ONE
  elseif (name == 'Mo') then
    ncore = 8
    nval = 3
    no(1) = 5
    lo(1) = 0
    zo(1) = 1*ONE
    no(2) = 5
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 4
    lo(3) = 2
    zo(3) = 5*ONE
  elseif (name == 'Tc') then
    ncore = 8
    nval = 3
    no(1) = 5
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 5
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 4
    lo(3) = 2
    zo(3) = 5*ONE
  elseif (name == 'Ru') then
    ncore = 8
    nval = 3
    no(1) = 5
    lo(1) = 0
    zo(1) = 1*ONE
    no(2) = 5
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 4
    lo(3) = 2
    zo(3) = 7*ONE
  elseif (name == 'Rh') then
    ncore = 8
    nval = 3
    no(1) = 5
    lo(1) = 0
    zo(1) = 1*ONE
    no(2) = 5
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 4
    lo(3) = 2
    zo(3) = 8*ONE
  elseif (name == 'Pd') then
    ncore = 8
    nval = 3
    no(1) = 5
    lo(1) = 0
    zo(1) = 0*ONE
    no(2) = 5
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 4
    lo(3) = 2
    zo(3) = 10*ONE
  elseif (name == 'Ag') then
    ncore = 8
    nval = 3
    no(1) = 5
    lo(1) = 0
    zo(1) = 1*ONE
    no(2) = 5
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 4
    lo(3) = 2
    zo(3) = 10*ONE
  elseif (name == 'Cd') then
    ncore = 8
    nval = 3
    no(1) = 5
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 5
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 4
    lo(3) = 2
    zo(3) = 10*ONE
  elseif (name == 'In') then
    ncore = 9
    nval = 3
    no(1) = 5
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 5
    lo(2) = 1
    zo(2) = 1*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Sn') then
    ncore = 9
    nval = 3
    no(1) = 5
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 5
    lo(2) = 1
    zo(2) = 2*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Sb') then
    ncore = 9
    nval = 3
    no(1) = 5
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 5
    lo(2) = 1
    zo(2) = 3*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Te') then
    ncore = 9
    nval = 3
    no(1) = 5
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 5
    lo(2) = 1
    zo(2) = 4*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'I ' .or. name == ' I') then
    ncore = 9
    nval = 3
    no(1) = 5
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 5
    lo(2) = 1
    zo(2) = 5*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Xe') then
    ncore = 9
    nval = 3
    no(1) = 5
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 5
    lo(2) = 1
    zo(2) = 6*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Cs') then
    ncore = 11
    nval = 3
    no(1) = 6
    lo(1) = 0
    zo(1) = 1*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Ba') then
    ncore = 11
    nval = 3
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'La') then
    ncore = 11
    nval = 3
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 1*ONE
  elseif (name == 'Ce') then
    ncore = 11
    nval = 4
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 1*ONE
    no(4) = 4
    lo(4) = 3
    zo(4) = 1*ONE
  elseif (name == 'Pr') then
    ncore = 11
    nval = 4
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 0*ONE
    no(4) = 4
    lo(4) = 3
    zo(4) = 3*ONE
  elseif (name == 'Nd') then
    ncore = 11
    nval = 4
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 0*ONE
    no(4) = 4
    lo(4) = 3
    zo(4) = 4*ONE
  elseif (name == 'Pm') then
    ncore = 11
    nval = 4
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 0*ONE
    no(4) = 4
    lo(4) = 3
    zo(4) = 5*ONE
  elseif (name == 'Sm') then
    ncore = 11
    nval = 4
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 0*ONE
    no(4) = 4
    lo(4) = 3
    zo(4) = 6*ONE
  elseif (name == 'Eu') then
    ncore = 11
    nval = 4
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 0*ONE
    no(4) = 4
    lo(4) = 3
    zo(4) = 7*ONE
  elseif (name == 'Gd') then
    ncore = 11
    nval = 4
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 1*ONE
    no(4) = 4
    lo(4) = 3
    zo(4) = 7*ONE
  elseif (name == 'Tb') then
    ncore = 11
    nval = 4
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 0*ONE
    no(4) = 4
    lo(4) = 3
    zo(4) = 9*ONE
  elseif (name == 'Dy') then
    ncore = 11
    nval = 4
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 0*ONE
    no(4) = 4
    lo(4) = 3
    zo(4) = 10*ONE
  elseif (name == 'Ho') then
    ncore = 11
    nval = 4
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 0*ONE
    no(4) = 4
    lo(4) = 3
    zo(4) = 11*ONE
  elseif (name == 'Er') then
    ncore = 11
    nval = 4
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 0*ONE
    no(4) = 4
    lo(4) = 3
    zo(4) = 12*ONE
  elseif (name == 'Tm') then
    ncore = 11
    nval = 4
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 0*ONE
    no(4) = 4
    lo(4) = 3
    zo(4) = 13*ONE
  elseif (name == 'Yb') then
    ncore = 11
    nval = 4
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 0*ONE
    no(4) = 4
    lo(4) = 3
    zo(4) = 14*ONE
  elseif (name == 'Lu') then
    ncore = 11
    nval = 4
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 1*ONE
    no(4) = 4
    lo(4) = 3
    zo(4) = 14*ONE
  elseif (name == 'Hf') then
    ncore = 12
    nval = 3
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 2*ONE
  elseif (name == 'Ta') then
    ncore = 12
    nval = 3
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 3*ONE
  elseif (name == 'W ' .or. name == ' W') then
    ncore = 12
    nval = 3
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 4*ONE
  elseif (name == 'Re') then
    ncore = 12
    nval = 3
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 5*ONE
  elseif (name == 'Os') then
    ncore = 12
    nval = 3
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 6*ONE
  elseif (name == 'Ir') then
    ncore = 12
    nval = 3
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 7*ONE
  elseif (name == 'Pt') then
    ncore = 12
    nval = 3
    no(1) = 6
    lo(1) = 0
    zo(1) = 1*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 9*ONE
  elseif (name == 'Au') then
    ncore = 12
    nval = 3
    no(1) = 6
    lo(1) = 0
    zo(1) = 1*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 10*ONE
  elseif (name == 'Hg') then
    ncore = 12
    nval = 3
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 5
    lo(3) = 2
    zo(3) = 10*ONE
  elseif (name == 'Tl') then
    ncore = 13
    nval = 3
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 1*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Pb') then
    ncore = 13
    nval = 3
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 2*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Bi') then
    ncore = 13
    nval = 3
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 3*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Po') then
    ncore = 13
    nval = 3
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 4*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'At') then
    ncore = 13
    nval = 3
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 5*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Rn') then
    ncore = 13
    nval = 3
    no(1) = 6
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 6
    lo(2) = 1
    zo(2) = 6*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Fr') then
    ncore = 15
    nval = 2
    no(1) = 7
    lo(1) = 0
    zo(1) = 1*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
  elseif (name == 'Ra') then
    ncore = 15
    nval = 3
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Ac') then
    ncore = 15
    nval = 4
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 1*ONE
    no(4) = 5
    lo(4) = 3
    zo(4) = 0*ONE
  elseif (name == 'Th') then
    ncore = 15
    nval = 4
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 1*ONE
    no(4) = 5
    lo(4) = 3
    zo(4) = 1*ONE
  elseif (name == 'Pa') then
    ncore = 15
    nval = 4
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 1*ONE
    no(4) = 5
    lo(4) = 3
    zo(4) = 2*ONE
  elseif (name == ' U' .or. name == 'U ') then
    ncore = 15
    nval = 4
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 1*ONE
    no(4) = 5
    lo(4) = 3
    zo(4) = 3*ONE
  elseif (name == 'Np') then
    ncore = 15
    nval = 4
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 0*ONE
    no(4) = 5
    lo(4) = 3
    zo(4) = 5*ONE
  elseif (name == 'Pu') then
    ncore = 15
    nval = 4
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 0*ONE
    no(4) = 5
    lo(4) = 3
    zo(4) = 6*ONE
  elseif (name == 'Am') then
    ncore = 15
    nval = 4
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 0*ONE
    no(4) = 5
    lo(4) = 3
    zo(4) = 7*ONE
  elseif (name == 'Cm') then
    ncore = 15
    nval = 4
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 1*ONE
    no(4) = 5
    lo(4) = 3
    zo(4) = 7*ONE
  elseif (name == 'Bk') then
    ncore = 15
    nval = 4
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 0*ONE
    no(4) = 5
    lo(4) = 3
    zo(4) = 9*ONE
  elseif (name == 'Cf') then
    ncore = 15
    nval = 4
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 0*ONE
    no(4) = 5
    lo(4) = 3
    zo(4) = 10*ONE
  elseif (name == 'Es') then
    ncore = 15
    nval = 4
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 0*ONE
    no(4) = 5
    lo(4) = 3
    zo(4) = 11*ONE
  elseif (name == 'Fm') then
    ncore = 15
    nval = 4
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 0*ONE
    no(4) = 5
    lo(4) = 3
    zo(4) = 12*ONE
  elseif (name == 'Md') then
    ncore = 15
    nval = 4
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 0*ONE
    no(4) = 5
    lo(4) = 3
    zo(4) = 13*ONE
  elseif (name == 'No') then
    ncore = 15
    nval = 4
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 0*ONE
    no(4) = 5
    lo(4) = 3
    zo(4) = 14*ONE
  elseif (name == 'Lr') then
    ncore = 15
    nval = 4
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 1*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 0*ONE
    no(4) = 5
    lo(4) = 3
    zo(4) = 14*ONE
  elseif (name == 'Rf') then
    ncore = 16
    nval = 3
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 2*ONE
  elseif (name == 'Db') then
    ncore = 16
    nval = 3
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 3*ONE
  elseif (name == 'Sg') then
    ncore = 16
    nval = 3
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 4*ONE
  elseif (name == 'Bh') then
    ncore = 16
    nval = 3
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 5*ONE
  elseif (name == 'Hs') then
    ncore = 16
    nval = 3
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 6*ONE
  elseif (name == 'Mt') then
    ncore = 16
    nval = 3
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 7*ONE
  elseif (name == 'Ds') then
    ncore = 16
    nval = 3
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 8*ONE
  elseif (name == 'Rg') then
    ncore = 16
    nval = 3
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 9*ONE
  elseif (name == 'Cn') then
    ncore = 16
    nval = 3
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 0*ONE
    no(3) = 6
    lo(3) = 2
    zo(3) = 10*ONE
  elseif (name == 'Nh') then
    ncore = 17
    nval = 3
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 1*ONE
    no(3) = 7
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Fl') then
    ncore = 17
    nval = 3
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 2*ONE
    no(3) = 7
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Mc') then
    ncore = 17
    nval = 3
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 3*ONE
    no(3) = 7
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Lv') then
    ncore = 17
    nval = 3
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 4*ONE
    no(3) = 7
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Ts') then
    ncore = 17
    nval = 3
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 5*ONE
    no(3) = 7
    lo(3) = 2
    zo(3) = 0*ONE
  elseif (name == 'Og') then
    ncore = 17
    nval = 3
    no(1) = 7
    lo(1) = 0
    zo(1) = 2*ONE
    no(2) = 7
    lo(2) = 1
    zo(2) = 6*ONE
    no(3) = 7
    lo(3) = 2
    zo(3) = 0*ONE
  else
    write(6,*)
    write(6,*)
    write(6,'("  Element ",a2," unknown")') name
    write(6,'("  Using ncore = 19, nval = 1")')
    write(6,*)
    ncore = 19
    nval = 1
    no(1) = 8
    lo(1) = 0
    zo(1) = 1*ONE

  endif

  return

end subroutine atom_p_tbl_config
