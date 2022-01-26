!>  Finds the core correction recommendations
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.3
!>  \date         23 September 2021.
!>  \copyright    GNU Public License v2

subroutine atom_psd_change_ifcore(ifcore_atm, ifcore_dat, ifcore_tab,    &
            ifcore)


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  ifcore_atm                       !<  0: no partial core correction; 1 : xc partial core correction; 2 : not recommended. From filein
  integer, intent(in)               ::  ifcore_dat                       !<  0: no partial core correction; 1 : xc partial core correction; 2 : not recommended. From fileae
  integer, intent(in)               ::  ifcore_tab                       !<  0: no partial core correction; 1 : xc partial core correction; 2 : not recommended. From p_tbl_psd_pcc

! output

  integer, intent(out)              ::  ifcore                           !<  0: no partial core correction; 1 : xc partial core correction; 2 : not recommended.

! local variables

  character(len=1)         ::  yesno


  if(ifcore_atm == ifcore_dat .and. ifcore_tab == ifcore_dat) then

    ifcore = ifcore_atm
    if(ifcore == 0) then
      write(6,*)
      write(6,*) '  With the current value of ifcore'
      write(6,*) '  the xc partial core correction will NOT be used'
      write(6,*)
    elseif(ifcore == 1) then
      write(6,*)
      write(6,*) '  With the current value of ifcore'
      write(6,*) '  the xc partial core correction will be used'
      write(6,*)
    elseif(ifcore == 2) then
      write(6,*)
      write(6,*) '  With the current value of ifcore the partial'
      write(6,*) '  core correction including Hartree'
      write(6,*) '  which is not recommended will be used.'
      write(6,*)
    endif
    write(6,*)
    write(6,*) ' Do you want to change that? (y/n)'
    write(6,*)

    read(5,*) yesno

    if(yesno == 'y' .or. yesno == 'Y') then
      write(6,*)
      write(6,*) '  Do you want to use xc partial core correction? (y/n)'
      write(6,*)

      read(5,*) yesno

      if(yesno == 'y' .or. yesno == 'Y') then
        ifcore = 1
      else
        ifcore = 0
      endif

    endif

  else

    write(6,*)
    if(ifcore_atm == 0) then
      write(6,*) '  From filein (atom.dat) the xc-pcc is not recommended.'
    elseif(ifcore_atm == 1) then
      write(6,*) '  From filein (atom.dat) the xc-pcc is recommended.'
    elseif(ifcore_atm == 2) then
      write(6,*) '  From filein (atom.dat) the hxc-pcc is recommended.'
      write(6,*) '  but you shouldn not use it.'
    endif
    if(ifcore_dat == 0) then
      write(6,*) '  From fileae (datafile.dat) the xc-pcc is not recommended.'
    elseif(ifcore_dat == 1) then
      write(6,*) '  From fileae (datafile.dat) the xc-pcc is recommended.'
    elseif(ifcore_dat == 2) then
      write(6,*) '  From fileae (datafile.dat) the hxc-pcc is recommended.'
      write(6,*) '  but you should not use it.'
    endif
    if(ifcore_tab == 0) then
      write(6,*) '  From default tables the xc-pcc is not recommended.'
    elseif(ifcore_dat == 1) then
      write(6,*) '  From default tables the xc-pcc is recommended.'
    endif

    write(6,*)
    write(6,*) '  Do you want to use xc partial core correction? (y/n)'
    write(6,*)

    read(5,*) yesno

    if(yesno == 'y' .or. yesno == 'Y') then
      ifcore = 1
    else
      ifcore = 0
    endif

  endif

  return

end subroutine atom_psd_change_ifcore
