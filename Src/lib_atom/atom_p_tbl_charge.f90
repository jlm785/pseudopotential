!>  Function determines the nuclear charge of an element.
!>
!>  \author       Jose Luis Martins
!>  \version      6.013
!>  \date         22 June 2021
!>  \copyright    GNU Public License v2

       subroutine atom_p_tbl_charge(name,iz)

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

!      cleanup and new interface, July 2019. JLM
!      jlm  version 6.013

       implicit none

!      input

       character(len=2), intent(in)      ::  name                        !<  chemical symbol of the atom

!      output

       integer, intent(out)              ::  iz                          !<  atomic number

       if (name == 'H ' .or. name == ' H') then
         iz = 1
       elseif (name == 'He') then
         iz = 2
       elseif (name == 'Li') then
         iz = 3
       elseif (name == 'Be') then
         iz = 4
       elseif (name == 'B ' .or. name == ' B') then
         iz = 5
       elseif (name == 'C ' .or. name == ' C') then
         iz = 6
       elseif (name == 'N ' .or. name == ' N') then
         iz = 7
       elseif (name == 'O ' .or. name == ' O') then
         iz = 8
       elseif (name == 'F ' .or. name == ' F') then
         iz = 9
       elseif (name == 'Ne') then
         iz = 10
       elseif (name == 'Na') then
         iz = 11
       elseif (name == 'Mg') then
         iz = 12
       elseif (name == 'Al') then
         iz = 13
       elseif (name == 'Si') then
         iz = 14
       elseif (name == 'P ' .or. name == ' P') then
         iz = 15
       elseif (name == 'S ' .or. name == ' S') then
         iz = 16
       elseif (name == 'Cl') then
         iz = 17
       elseif (name == 'Ar') then
         iz = 18
       elseif (name == 'K ' .or. name == ' K') then
         iz = 19
       elseif (name == 'Ca') then
         iz = 20
       elseif (name == 'Sc') then
         iz = 21
       elseif (name == 'Ti') then
         iz = 22
       elseif (name == 'V ' .or. name == ' V') then
         iz = 23
       elseif (name == 'Cr') then
         iz = 24
       elseif (name == 'Mn') then
         iz = 25
       elseif (name == 'Fe') then
         iz = 26
       elseif (name == 'Co') then
         iz = 27
       elseif (name == 'Ni') then
         iz = 28
       elseif (name == 'Cu') then
         iz = 29
       elseif (name == 'Zn') then
         iz = 30
       elseif (name == 'Ga') then
         iz = 31
       elseif (name == 'Ge') then
         iz = 32
       elseif (name == 'As') then
         iz = 33
       elseif (name == 'Se') then
         iz = 34
       elseif (name == 'Br') then
         iz = 35
       elseif (name == 'Kr') then
         iz = 36
       elseif (name == 'Rb') then
         iz = 37
       elseif (name == 'Sr') then
         iz = 38
       elseif (name == 'Y ' .or. name == ' Y') then
         iz = 39
       elseif (name == 'Zr') then
         iz = 40
       elseif (name == 'Nb') then
         iz = 41
       elseif (name == 'Mo') then
         iz = 42
       elseif (name == 'Tc') then
         iz = 43
       elseif (name == 'Ru') then
         iz = 44
       elseif (name == 'Rh') then
         iz = 45
       elseif (name == 'Pd') then
         iz = 46
       elseif (name == 'Ag') then
         iz = 47
       elseif (name == 'Cd') then
         iz = 48
       elseif (name == 'In') then
         iz = 49
       elseif (name == 'Sn') then
         iz = 50
       elseif (name == 'Sb') then
         iz = 51
       elseif (name == 'Te') then
         iz = 52
       elseif (name == 'I ' .or. name == ' I') then
         iz = 53
       elseif (name == 'Xe') then
         iz = 54
       elseif (name == 'Cs') then
         iz = 55
       elseif (name == 'Ba') then
         iz = 56
       elseif (name == 'La') then
         iz = 57
       elseif (name == 'Ce') then
         iz = 58
       elseif (name == 'Pr') then
         iz = 59
       elseif (name == 'Nd') then
         iz = 60
       elseif (name == 'Pm') then
         iz = 61
       elseif (name == 'Sm') then
         iz = 62
       elseif (name == 'Eu') then
         iz = 63
       elseif (name == 'Gd') then
         iz = 64
       elseif (name == 'Tb') then
         iz = 65
       elseif (name == 'Dy') then
         iz = 66
       elseif (name == 'Ho') then
         iz = 67
       elseif (name == 'Er') then
         iz = 68
       elseif (name == 'Tm') then
         iz = 69
       elseif (name == 'Yb') then
         iz = 70
       elseif (name == 'Lu') then
         iz = 71
       elseif (name == 'Hf') then
         iz = 72
       elseif (name == 'Ta') then
         iz = 73
       elseif (name == 'W ' .or. name == ' W') then
         iz = 74
       elseif (name == 'Re') then
         iz = 75
       elseif (name == 'Os') then
         iz = 76
       elseif (name == 'Ir') then
         iz = 77
       elseif (name == 'Pt') then
         iz = 78
       elseif (name == 'Au') then
         iz = 79
       elseif (name == 'Hg') then
         iz = 80
       elseif (name == 'Tl') then
         iz = 81
       elseif (name == 'Pb') then
         iz = 82
       elseif (name == 'Bi') then
         iz = 83
       elseif (name == 'Po') then
         iz = 84
       elseif (name == 'At') then
         iz = 85
       elseif (name == 'Rn') then
         iz = 86
       elseif (name == 'Fr') then
         iz = 87
       elseif (name == 'Ra') then
         iz = 88
       elseif (name == 'Ac') then
         iz = 89
       elseif (name == 'Th') then
         iz = 90
       elseif (name == 'Pa') then
         iz = 91
       elseif (name == ' U' .or. name == 'U ') then
         iz = 92
       elseif (name == 'Np') then
         iz = 93
       elseif (name == 'Pu') then
         iz = 94
       elseif (name == 'Am') then
         iz = 95
       elseif (name == 'Cm') then
         iz = 96
       elseif (name == 'Bk') then
         iz = 97
       elseif (name == 'Cf') then
         iz = 98
       elseif (name == 'Es') then
         iz = 99
       elseif (name == 'Fm') then
         iz = 100
       elseif (name == 'Md') then
         iz = 101
       elseif (name == 'No') then
         iz = 102
       elseif (name == 'Lr') then
         iz = 103
       elseif (name == 'Rf') then
         iz = 104
       elseif (name == 'Db') then
         iz = 105
       elseif (name == 'Sg') then
         iz = 106
       elseif (name == 'Bh') then
         iz = 107
       elseif (name == 'Hs') then
         iz = 108
       elseif (name == 'Mt') then
         iz = 109
       elseif (name == 'Ds') then
         iz = 110
       elseif (name == 'Rg') then
         iz = 111
       elseif (name == 'Cn') then
         iz = 112
       elseif (name == 'Nh') then
         iz = 113
       elseif (name == 'Fl') then
         iz = 114
       elseif (name == 'Mc') then
         iz = 115
       elseif (name == 'Lv') then
         iz = 116
       elseif (name == 'Ts') then
         iz = 117
       elseif (name == 'Og') then
         iz = 118
       else
         write(6,*)
         write(6,*)
         write(6,'("  Element ",a2," unknown")') name
         write(6,'("  Using charge 200")')
         write(6,*)
         iz = 200
       endif

       return

       end subroutine atom_p_tbl_charge
