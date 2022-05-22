!>  Suggests core radii for the Troullier-Martins pseudopotential.
!>  This are on the harder (small radii) side.
!>
!>  \author       Jose Luis Martins
!>  \version      6.013
!>  \date         22 June 2021
!>  \copyright    GNU Public License v2

subroutine atom_p_tbl_psd_tm2(name, rc, status)

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
! jlm  version 6.013

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)! input

  integer, parameter                ::  lc = 4                           !  maximum valence angular momentum

! input

  character(len=2), intent(in)      ::  name                             !<  chemical symbol of the atom

! output

  real(REAL64), intent(out)         ::  rc(0:lc)                         !<  core radii for angular momentum l
  character(len=30), intent(out)    ::  status                           !<  quality of pseudopotential

! constant

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64

! counter

  integer        ::  l


  do l = 0,lc
    rc(l) = ZERO
  enddo

  status = '                              '

  if (name == 'H ' .or. name == ' H') then
    rc(0) = 0.50
    status = 'fully tested                  '
  elseif (name == 'He') then
    rc(0) = 0.80
    status = 'hard tested                   '
  elseif (name == 'Li') then
    rc(0) = 2.00
    rc(1) = 2.00
    rc(2) = 2.00
    status = 'hard tested                   '
 elseif (name == 'Be') then
    rc(0) = 1.70
    rc(1) = 1.70
    rc(2) = 1.70
    status = 'hard tested                   '
  elseif (name == 'B ' .or. name == ' B') then
    rc(0) = 1.40
    rc(1) = 1.40
    rc(2) = 1.40
    status = 'hard tested                   '
  elseif (name == 'C ' .or. name == ' C') then
    rc(0) = 1.10
    rc(1) = 1.10
    rc(2) = 1.10
    status = 'hard tested   (triple bond)   '
  elseif (name == 'N ' .or. name == ' N') then
    rc(0) = 1.00
    rc(1) = 1.00
    rc(2) = 1.00
    status = 'hard used      (triple bond)  '
  elseif (name == 'O ' .or. name == ' O') then
    rc(0) = 1.10
    rc(1) = 1.10
    rc(2) = 1.10
    status = 'hard used                     '
  elseif (name == 'F ' .or. name == ' F') then
    rc(0) = 1.15
    rc(1) = 1.15
    rc(2) = 1.15
    status = 'hard guess                    '
  elseif (name == 'Ne') then
    rc(0) = 2.0
    rc(1) = 2.0
    status = 'hard who cares                '
  elseif (name == 'Na') then
    rc(0) = 2.30
    rc(1) = 2.50
    rc(2) = 2.50
    status = 'hard used      (oxides)       '
  elseif (name == 'Mg') then
    rc(0) = 2.80
    rc(1) = 2.80
    rc(2) = 2.80
    status = 'hard used      (oxides)       '
  elseif (name == 'Al') then
    rc(0) = 2.28
    rc(1) = 2.28
    rc(2) = 2.28
    status = 'hard used      (ionic)        '
  elseif (name == 'Si') then
    rc(0) = 1.80
    rc(1) = 1.80
    rc(2) = 1.80
    status = 'hard used      (ionic)        '
  elseif (name == 'P ' .or. name == ' P') then
    rc(0) = 1.75
    rc(1) = 1.75
    rc(2) = 1.75
    status = 'hard guess                    '
  elseif (name == 'S ' .or. name == ' S') then
    rc(0) = 1.65
    rc(1) = 1.65
    rc(2) = 1.65
    status = 'hard guess                    '
  elseif (name == 'Cl') then
    rc(0) = 1.65
    rc(1) = 1.65
    rc(2) = 1.65
    status = 'hard guess                    '
  elseif (name == 'Ar') then
    rc(0) = 3.0
    rc(1) = 3.0
    rc(2) = 3.0
    status = 'hard who cares                '
  elseif (name == 'K ' .or. name == ' K') then
    rc(0) = 3.50
    rc(1) = 3.70
    rc(2) = 3.50
    status = 'hard guess                    '
  elseif (name == 'Ca') then
    rc(0) = 3.10
    rc(1) = 3.10
    rc(2) = 3.10
    status = 'hard guess                    '
  elseif (name == 'Sc') then
    rc(0) = 2.60
    rc(1) = 2.60
    rc(2) = 2.60
    status = 'hard guess                    '
  elseif (name == 'Ti') then
    rc(0) = 2.45
    rc(1) = 2.85
    rc(2) = 2.15
    status = 'hard based on used            '
  elseif (name == 'V ' .or. name == ' V') then
    rc(0) = 2.25
    rc(1) = 2.60
    rc(2) = 2.25
    status = 'hard guess                    '
  elseif (name == 'Cr') then
    rc(0) = 2.25
    rc(1) = 2.60
    rc(2) = 2.25
    status = 'hard guess                    '
  elseif (name == 'Mn') then
    rc(0) = 2.25
    rc(1) = 2.50
    rc(2) = 2.25
    status = 'hard guess                    '
  elseif (name == 'Fe') then
    rc(0) = 2.10
    rc(1) = 2.30
    rc(2) = 2.10
    status = 'hard guess                    '
  elseif (name == 'Co') then
    rc(0) = 2.10
    rc(1) = 2.25
    rc(2) = 2.10
    status = 'hard guess                    '
  elseif (name == 'Ni') then
    rc(0) = 2.10
    rc(1) = 2.25
    rc(2) = 2.10
    status = 'hard guess                    '
  elseif (name == 'Cu') then
    rc(0) = 2.00
    rc(1) = 2.20
    rc(2) = 2.00
    status = 'hard based on used            '
  elseif (name == 'Zn') then
    rc(0) = 2.20
    rc(1) = 2.20
    rc(2) = 2.20
    status = 'hard guess                    '
  elseif (name == 'Ga') then
    rc(0) = 2.00
    rc(1) = 2.40
    rc(2) = 2.00
    status = 'hard guess  try 3d in valence!'
  elseif (name == 'Ge') then
    rc(0) = 2.20
    rc(1) = 2.20
    rc(2) = 2.20
    status = 'hard based on used            '
  elseif (name == 'As') then
    rc(0) = 2.00
    rc(1) = 2.00
    rc(2) = 2.00
    status = 'hard guess                    '
  elseif (name == 'Se') then
    rc(0) = 1.80
    rc(1) = 1.80
    rc(2) = 1.80
    status = 'hard guess                    '
  elseif (name == 'Br') then
    rc(0) = 1.85
    rc(1) = 1.85
    rc(2) = 1.85
    status = 'hard guess                    '
  elseif (name == 'Kr') then
    rc(0) = 1.85
    rc(1) = 1.85
    rc(2) = 1.85
    status = 'hard who cares                '
  elseif (name == 'Rb') then
    rc(0) = 3.80
    rc(1) = 4.00
    rc(2) = 3.80
    status = 'hard used                     '
  elseif (name == 'Sr') then
    rc(0) = 3.50
    rc(1) = 3.70
    rc(2) = 3.50
    status = 'hard guess                    '
  elseif (name == 'Y ' .or. name == ' Y') then
    rc(0) = 3.00
    rc(1) = 3.40
    rc(2) = 2.50
    status = 'hard used                     '
  elseif (name == 'Zr') then
    rc(0) = 2.75
    rc(1) = 2.95
    rc(2) = 2.75
    status = 'hard guess                    '
  elseif (name == 'Nb') then
    rc(0) = 2.40
    rc(1) = 2.75
    rc(2) = 2.25
    status = 'hard guess                    '
  elseif (name == 'Mo') then
    rc(0) = 2.35
    rc(1) = 2.65
    rc(2) = 2.25
    status = 'hard guess                    '
  elseif (name == 'Tc') then
    rc(0) = 2.30
    rc(1) = 2.60
    rc(2) = 2.30
    status = 'hard who cares                '
  elseif (name == 'Ru') then
    rc(0) = 2.30
    rc(1) = 2.50
    rc(2) = 2.30
    status = 'hard guess                    '
  elseif (name == 'Rh') then
    rc(0) = 2.30
    rc(1) = 2.50
    rc(2) = 2.30
    status = 'hard guess                    '
  elseif (name == 'Pd') then
    rc(0) = 2.30
    rc(1) = 2.50
    rc(2) = 2.30
    status = 'hard guess                    '
  elseif (name == 'Ag') then
    rc(0) = 2.40
    rc(1) = 2.50
    rc(2) = 2.40
    status = 'hard guess                    '
  elseif (name == 'Cd') then
    rc(0) = 2.45
    rc(1) = 2.45
    rc(2) = 2.45
    status = 'hard guess                    '
  elseif (name == 'In') then
    rc(0) = 2.70
    rc(1) = 2.70
    rc(2) = 2.70
    status = 'hard guess                    '
  elseif (name == 'Sn') then
    rc(0) = 2.30
    rc(1) = 2.30
    rc(2) = 2.30
    status = 'hard guess                    '
  elseif (name == 'Sb') then
    rc(0) = 2.40
    rc(1) = 2.40
    rc(2) = 2.40
    status = 'hard guess                    '
  elseif (name == 'Te') then
    rc(0) = 2.60
    rc(1) = 2.60
    rc(2) = 2.60
    status = 'hard guess                    '
  elseif (name == 'I ' .or. name == ' I') then
    rc(0) = 2.30
    rc(1) = 2.30
    rc(2) = 2.30
    status = 'hard guess                    '
  elseif (name == 'Xe') then
    rc(0) = 3.20
    rc(1) = 3.20
    rc(2) = 3.20
    status = 'hard who cares                '
  elseif (name == 'Cs') then
    rc(0) = 4.00
    rc(1) = 4.50
    rc(2) = 4.00
    status = 'hard guess based on used      '
  elseif (name == 'Ba') then
    rc(0) = 3.70
    rc(1) = 4.20
    rc(2) = 3.70
    status = 'hard guess                    '
  elseif (name == 'La') then
    rc(0) = 3.70
    rc(1) = 3.90
    rc(2) = 3.70
    status = 'hard guess                    '
  elseif (name == 'Ce') then
    rc(0) = 3.30
    rc(1) = 3.60
    rc(2) = 2.80
    rc(3) = 2.90
    status = 'hard guess  based on used     '
  elseif (name == 'Pr') then
    rc(0) = 3.00
    rc(1) = 3.50
    rc(2) = 3.00
    rc(3) = 3.00
    status = 'hard guess                    '
  elseif (name == 'Nd') then
    rc(0) = 3.00
    rc(1) = 3.50
    rc(2) = 3.00
    rc(3) = 3.00
    status = 'hard guess                    '
  elseif (name == 'Pm') then
    rc(0) = 3.00
    rc(1) = 3.50
    rc(2) = 3.00
    rc(3) = 3.00
    status = 'hard guess                    '
  elseif (name == 'Sm') then
    rc(0) = 3.00
    rc(1) = 3.50
    rc(2) = 3.00
    rc(3) = 3.00
    status = 'hard guess                    '
  elseif (name == 'Eu') then
    rc(0) = 3.00
    rc(1) = 3.50
    rc(2) = 3.00
    rc(3) = 3.00
    status = 'hard guess                    '
  elseif (name == 'Gd') then
    rc(0) = 3.00
    rc(1) = 3.50
    rc(2) = 3.00
    rc(3) = 3.00
    status = 'hard guess                    '
  elseif (name == 'Tb') then
    rc(0) = 2.90
    rc(1) = 3.40
    rc(2) = 2.90
    rc(3) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Dy') then
    rc(0) = 2.90
    rc(1) = 3.40
    rc(2) = 2.90
    rc(3) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Ho') then
    rc(0) = 2.90
    rc(1) = 3.40
    rc(2) = 2.90
    rc(3) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Er') then
    rc(0) = 2.90
    rc(1) = 3.30
    rc(2) = 2.90
    rc(3) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Tm') then
    rc(0) = 2.90
    rc(1) = 3.30
    rc(2) = 2.90
    rc(3) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Yb') then
    rc(0) = 2.90
    rc(1) = 3.30
    rc(2) = 2.90
    rc(3) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Lu') then
    rc(0) = 2.90
    rc(1) = 3.20
    rc(2) = 2.90
    rc(3) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Hf') then
    rc(0) = 2.90
    rc(1) = 3.20
    rc(2) = 2.90
    status = 'hard guess (try 4f in valence)'
  elseif (name == 'Ta') then
    rc(0) = 2.50
    rc(1) = 2.70
    rc(2) = 2.50
    status = 'hard guess                    '
  elseif (name == 'W ' .or. name == ' W') then
    rc(0) = 2.40
    rc(1) = 2.60
    rc(2) = 2.40
    status = 'hard guess                    '
  elseif (name == 'Re') then
    rc(0) = 2.40
    rc(1) = 2.60
    rc(2) = 2.40
    status = 'hard guess                    '
  elseif (name == 'Os') then
    rc(0) = 2.35
    rc(1) = 2.55
    rc(2) = 2.35
    status = 'hard guess                    '
  elseif (name == 'Ir') then
    rc(0) = 2.40
    rc(1) = 2.60
    rc(2) = 2.40
    status = 'hard guess                    '
  elseif (name == 'Pt') then
    rc(0) = 2.40
    rc(1) = 2.60
    rc(2) = 2.40
    status = 'hard guess                    '
  elseif (name == 'Au') then
    rc(0) = 2.40
    rc(1) = 2.60
    rc(2) = 2.40
    status = 'hard guess                    '
  elseif (name == 'Hg') then
    rc(0) = 2.50
    rc(1) = 2.50
    rc(2) = 2.50
    status = 'hard guess                    '
  elseif (name == 'Tl') then
    rc(0) = 2.80
    rc(1) = 2.80
    rc(2) = 3.10
    status = 'hard guess                    '
  elseif (name == 'Pb') then
    rc(0) = 3.00
    rc(1) = 3.00
    rc(2) = 3.00
    status = 'hard guess based on used      '
  elseif (name == 'Bi') then
    rc(0) = 2.50
    rc(1) = 2.50
    rc(2) = 2.80
    status = 'hard guess                    '
  elseif (name == 'Po') then
    rc(0) = 2.70
    rc(1) = 2.70
    rc(2) = 3.00
    status = 'hard guess                    '
  elseif (name == 'At') then
    rc(0) = 2.70
    rc(1) = 2.70
    rc(2) = 2.70
    status = 'hard guess                    '
  elseif (name == 'Rn') then
    rc(0) = 3.40
    rc(1) = 3.40
    rc(2) = 3.40
    status = 'hard who cares                '
  elseif (name == 'Fr') then
    rc(0) = 5.00
    rc(1) = 5.50
    status = 'hard who cares                '
  elseif (name == 'Ra') then
    rc(0) = 4.20
    rc(1) = 4.20
    rc(2) = 4.20
    status = 'hard guess                    '
  elseif (name == 'Ac') then
    rc(0) = 3.20
    rc(1) = 4.00
    rc(2) = 2.80
    rc(3) = 2.80
    status = 'hard guess                    '
  elseif (name == 'Th') then
    rc(0) = 3.20
    rc(1) = 4.00
    rc(2) = 2.80
    rc(3) = 2.80
    status = 'hard guess                    '
  elseif (name == 'Pa') then
    rc(0) = 3.20
    rc(1) = 4.00
    rc(2) = 2.80
    rc(3) = 2.80
    status = 'hard guess                    '
  elseif (name == 'U ' .or. name == ' U') then
    rc(0) = 3.10
    rc(1) = 3.90
    rc(2) = 2.80
    rc(3) = 2.80
    status = 'hard guess                    '
  elseif (name == 'Np') then
    rc(0) = 3.00
    rc(1) = 3.90
    rc(2) = 2.80
    rc(3) = 2.80
    status = 'hard guess                    '
  elseif (name == 'Pu') then
    rc(0) = 3.00
    rc(1) = 3.90
    rc(2) = 2.80
    rc(3) = 2.80
    status = 'hard guess                    '
  elseif (name == 'Am') then
    rc(0) = 2.90
    rc(1) = 3.80
    rc(2) = 2.90
    rc(3) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Cm') then
    rc(0) = 2.90
    rc(1) = 3.80
    rc(2) = 2.90
    rc(3) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Bk') then
    rc(0) = 2.90
    rc(1) = 3.80
    rc(2) = 2.90
    rc(3) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Cf') then
    rc(0) = 2.90
    rc(1) = 3.80
    rc(2) = 2.90
    rc(3) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Es') then
    rc(0) = 3.00
    rc(1) = 3.80
    rc(2) = 2.90
    rc(3) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Fm') then
    rc(0) = 2.90
    rc(1) = 3.80
    rc(2) = 2.90
    rc(3) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Md') then
    status = 'hard guess                    '
    rc(0) = 2.90
    rc(1) = 3.80
    rc(2) = 2.90
    rc(3) = 2.90
  elseif (name == 'No') then
    rc(0) = 2.90
    rc(1) = 3.80
    rc(2) = 2.90
    rc(3) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Lr') then
    rc(0) = 2.90
    rc(1) = 3.80
    rc(2) = 2.90
    rc(3) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Rf') then
    rc(0) = 2.90
    rc(1) = 3.80
    rc(2) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Db') then
    rc(0) = 2.90
    rc(1) = 3.80
    rc(2) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Sg') then
    rc(0) = 2.90
    rc(1) = 3.80
    rc(2) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Bh') then
    rc(0) = 2.90
    rc(1) = 3.80
    rc(2) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Hs') then
    rc(0) = 2.90
    rc(1) = 3.80
    rc(2) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Mt') then
    rc(0) = 2.90
    rc(1) = 3.80
    rc(2) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Ds') then
    rc(0) = 2.90
    rc(1) = 3.80
    rc(2) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Rg') then
    rc(0) = 2.90
    rc(1) = 3.80
    rc(2) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Cn') then
    rc(0) = 2.90
    rc(1) = 3.80
    rc(2) = 2.90
    status = 'hard guess                    '
  elseif (name == 'Nh') then
    rc(0) = 3.00
    rc(1) = 3.80
    rc(2) = 3.30
    status = 'hard guess                    '
  elseif (name == 'Fl') then
    rc(0) = 3.00
    rc(1) = 3.80
    rc(2) = 3.00
    status = 'hard guess                    '
  elseif (name == 'Mc') then
    rc(0) = 3.00
    rc(1) = 3.80
    rc(2) = 3.00
    status = 'hard guess                    '
  elseif (name == 'Lv') then
    rc(0) = 3.00
    rc(1) = 3.80
    rc(2) = 3.00
    status = 'hard guess                    '
  elseif (name == 'Ts') then
    rc(0) = 3.00
    rc(1) = 3.80
    rc(2) = 3.00
    status = 'hard guess                    '
  elseif (name == 'Og') then
    rc(0) = 3.00
    rc(1) = 3.80
    rc(2) = 3.00
    status = 'hard guess                    '
  else
    write(6,*)
    write(6,*)
    write(6,'("  Element ",a2," unknown")') name
    write(6,'("  Using rc = 0.0")')
    write(6,*)
    status = 'unknown element               '
  endif

  return

end subroutine atom_p_tbl_psd_tm2
