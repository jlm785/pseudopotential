!>  computes and prints the total energy from the
!>  electron charge density.
!>
!>  \author       Sverre Froyen,  Norm Troullier, Jose Luis Martins
!>  \version      6.0.3
!>  \date         1980s, 22 June 2021, 12 September 2021.
!>  \copyright    GNU Public License v2

  subroutine atom_atm_etotal(itype, nameat, norb,                        &
       no, lo, iso, zo, etot, ev, ek, ep,                                &
       iowrite, mxdorb)


!>  tot(i)    i=1,10 contains various contributions to the total
!>            energy.
!>            (1)   sum of eigenvalues ev
!>            (2)   sum of orbital kinetic energies ek
!>            (3)   el-ion interaction from sum of orbital
!>                  potential energies ep
!>            (4)   electrostatic el-el interaction  (from velect)
!>            (5)   vxc (exchange-correlation) correction to sum
!>                  of eigenvalues                   (from velect)
!>            (6)   3 * vc - 4 * ec
!>                  correction term for virial theorem
!>                  when correlation is included     (from velect)
!>            (7)   exchange and correlation energy  (from velect)
!>            (8)   kinetic energy from eigenvalues  (1,3,4,5)
!>            (9)   potential energy
!>            (10)  total energy


! converted to fortran 90, March 1st 2018
! cleanup and new interface, July 2019. JLM
! jlm  version 6.01
! zsh removed 6 August 2021. JLM
! so->iso. 12 September 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdorb                           !<  dimension of the number of orbitals

  integer, intent(in)               ::  itype                            !<  type of calculation (-1 signals end of calculation)

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbot of the atom

  integer, intent(in)               ::  norb                             !<  number of orbitals

  integer, intent(in)               ::  no(mxdorb)                       !<  principal quantum number n
  integer, intent(in)               ::  lo(mxdorb)                       !<  angular quantum number l
  integer, intent(in)               ::  iso(mxdorb)                      !<  2*spin or 2*(j-l)
  real(REAL64), intent(in)          ::  zo(mxdorb)                       !<  orbital occupation

  real(REAL64), intent(in)          ::  ev(mxdorb)                       !<  orbital energy
  real(REAL64), intent(in)          ::  ek(mxdorb)                       !<  orbital kinetic energy
  real(REAL64), intent(in)          ::  ep(mxdorb)                       !<  orbital potential energy

! input and output

  real(REAL64), intent(inout)       ::  etot(10)                         !<  total energy components

! local variables

  real(REAL64)   ::  vsum

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  character(len=1), parameter  ::  IL(7) = (/'s','p','d','f','g','h','i'/)

! counters

  integer     ::  i


! sum up eigenvalues ev, kinetic energies ek, and
! el-ion interaction ep

  etot(1) = ZERO
  etot(2) = ZERO
  etot(3) = ZERO
  do i=1,norb
    etot(1) = etot(1) + zo(i)*ev(i)
    etot(2) = etot(2) + zo(i)*ek(i)
    etot(3) = etot(3) + zo(i)*ep(i)
  enddo

! kinetic energy

  etot(8) = etot(1) - etot(3) - 2*etot(4) - etot(5)

! potential energy

  etot(9) = etot(3) + etot(4) + etot(7)

! total energy

  etot(10) = etot(1) - etot(4) - etot(5) + etot(7)

! printout

  write(iowrite,'(//,1x,a2," output data for orbitals",/,1x,             &
     &     28("-"),//," nl    s      occ",9x,"eigenvalue",4x,            &
     &     "kinetic energy",6x,"pot energy",/)') nameat

  do i=1,norb
    if(lo(i) < 7) then
       write(iowrite,'(1x,i1,a1,f6.1,f10.4,3f17.8)')                     &
        no(i), IL(lo(i)+1), iso(i)*0.5, zo(i), ev(i), ek(i), ep(i)
    else
       write(iowrite,'(1x,i1,a1,f6.1,f10.4,3f17.8)')                     &
              no(i), 'L', iso(i)*0.5, zo(i), ev(i), ek(i), ep(i)
    endif
  enddo

  write(iowrite,'(//," total energies",/,1x,14("-"),/,                   &
     &    /," sum of eigenvalues        =",f18.8,                        &
     &    /," kinetic energy from ek    =",f18.8,                        &
     &    /," el-ion interaction energy =",f18.8,                        &
     &    /," el-el  interaction energy =",f18.8,                        &
     &    /," vxc    correction         =",f18.8,                        &
     &    /," virial correction         =",f18.8,                        &
     &    /," exchange + corr energy    =",f18.8,                        &
     &    /," kinetic energy from ev    =",f18.8,                        &
     &    /," potential energy          =",f18.8,/,1x,45("-"),           &
     &    /," total energy              =",f18.8)') (etot(i),i=1,10)

  if (itype < 4) then

!   virial theorem

    vsum = 2*etot(8) + etot(9) + etot(6)

    write(iowrite,'(//," virial theorem(nonrelativistic)",               &
     &   /,1x,14("-"),/,                                                 &
     &   /," kinetic energy  *  2      =",f18.8,                         &
     &   /," potential energy          =",f18.8,                         &
     &   /," virial correction         =",f18.8,/,1x,45("-"),            &
     &   /," virial sum                =",f18.8,/)')                     &
         2*etot(8),etot(9),etot(6),vsum

  endif

  return

end subroutine atom_atm_etotal

