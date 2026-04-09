!>  Finds the core correction recommendations
!>
!>  \author       Jose Luis Martins
!>  \version      6.1.0
!>  \date         26 March 2026.
!>  \copyright    GNU Public License v2

subroutine atom_psd_change_cfac(ifcore, cfac, rcfac, cfactau)


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  ifcore                           !<  0: no partial core correction; 1 : xc partial core correction; 2 : not recommended.

! input and output

  real(REAL64), intent(inout)       ::  cfac                             !<  criteria for pseudo-core charge
  real(REAL64), intent(inout)       ::  rcfac                            !<  pseudo-core radius
  real(REAL64), intent(inout)       ::  cfactau                          !<  core correction factor for kinetic energy density

! local variables

  character(len=1)         ::  yesno

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  real(REAL64), parameter    ::  EPS = 1.0E-4_REAL64


  if(ifcore >= 0) then

    if (rcfac <= ZERO) then
      write(6,*)
      write(6,'("  The current value of cfac is ",f10.3)') cfac
      write(6,*) '  The xc partial core correction will be set'
      write(6,*) '  at the radius where the core charge density'
      write(6,*) '  is cfac times the valence cahrge density'
      write(6,*) '  (there is a correction for ionized reference atoms)'
      write(6,*)
      write(6,*) '  Do you want to change cfac? (y/n)'
      write(6,*)

      read(5,*) yesno

      if(yesno == 'y' .or. yesno == 'Y') then
        write(6,*)
        write(6,*) '  Enter new value of cfac'
        write(6,*)

        read(5,*) cfac

        if(cfac < 0.1 .or. cfac > 10.0) then
          write(6,*)
          write(6,*) '  The value looks strange, enter it again'
          write(6,*)

          read(5,*) cfac

          if(cfac < ZERO) then
            write(6,*)
            write(6,*) '  cfac cannot be negative, setting it to 1'
            write(6,*)
            cfac = ONE
          endif
        endif
      endif
    else
      write(6,*)
      write(6,'("  The current value of rcfac is ",f10.3)') rcfac
      write(6,*) '  The xc partial core correction will be set at that valur'
      write(6,*)
      write(6,*) '  Do you want to change cfac? (y/n)'
      write(6,*)

      read(5,*) yesno

      if(yesno == 'y' .or. yesno == 'Y') then
        write(6,*)
        write(6,*) '  Enter new value of rcfac'
        write(6,*)

        read(5,*) rcfac

        if(rcfac < 0.1 .or. rcfac > 10.0) then
          write(6,*)
          write(6,*) '  The value looks strange,  check it and enter it again'
          write(6,*)

          read(5,*) rcfac

          if(rcfac <= ZERO) then
            write(6,*)
            write(6,*) '  rcfac cannot be negative, not changing it'
            write(6,*)
          endif
        endif
      endif
    endif

  endif


  write(6,*)
  write(6,'("  The current value of cfactau is ",f10.3)') cfactau
  write(6,*) '  The kinetic energy density partial core correction will be set'
  write(6,*) '  at the radius where the core charge density'
  write(6,*) '  is cfac times the valence cahrge density'
  write(6,*) '  (there is a correction for ionized reference atoms)'
  write(6,*)
  write(6,*) '  Do you want to change cfactau? (y/n)'
  write(6,*)

  read(5,*) yesno

  if(yesno == 'y' .or. yesno == 'Y') then
    write(6,*)
    write(6,*) '  Enter new value of cfactau'
    write(6,*)

    read(5,*) cfactau

    if(cfactau < 0.1 .or. cfactau > 10.0) then
      write(6,*)
      write(6,*) '  The value looks strange, check it and enter it again'
      write(6,*)

      read(5,*) cfactau

      if(cfactau <= ZERO) then
        write(6,*)
        write(6,*) '  cfactau cannot be negative, setting it to 1'
        write(6,*)
        cfactau = ONE
      endif
    endif
  endif


  return

end subroutine atom_psd_change_cfac
