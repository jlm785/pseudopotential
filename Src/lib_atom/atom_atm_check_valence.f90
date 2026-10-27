!>  Checks the valence orbitals read from atom.dat:
!>  0 <= l < n, no orbital repeated, and no valence orbital
!>  that is also a core orbital.  Stops on error.
!>  Repeated angular momenta with different n are allowed.
!>
!>  \author       claude
!>  \version      6.1.0
!>  \date         27 October 2026.
!>  \copyright    GNU Public License v2

subroutine atom_atm_check_valence(ncore, nocore, locore, nval, ni, li,   &
      cname)

! written 27 October 2026. claude

  implicit none

! input

  integer, intent(in)               ::  ncore                            !<  number of core orbitals to check against
  integer, intent(in)               ::  nocore(*)                        !<  principal quantum number of core orbitals
  integer, intent(in)               ::  locore(*)                        !<  angular momentum of core orbitals

  integer, intent(in)               ::  nval                             !<  number of valence orbitals
  integer, intent(in)               ::  ni(nval)                         !<  principal quantum number of valence orbitals
  integer, intent(in)               ::  li(nval)                         !<  angular momentum of valence orbitals

  character(len=*), intent(in)      ::  cname                            !<  name of the calling subroutine

! counters

  integer     ::  i, j


  do i = 1,nval

    if(li(i) < 0 .or. ni(i) <= li(i)) then
      write(6,*)
      write(6,'("   Stopped in ",a,": valence orbital ",i3,              &
          &     " has n = ",i3," l = ",i3)') cname, i, ni(i), li(i)
      write(6,*) '  0 <= l < n is required'
      write(6,*)

      STOP

    endif

    do j = 1,i-1
      if(ni(j) == ni(i) .and. li(j) == li(i)) then
        write(6,*)
        write(6,'("   Stopped in ",a,": valence orbitals ",i3," and ",i3, &
            &     " are the same, n = ",i3," l = ",i3)') cname, j, i, ni(i), li(i)
        write(6,*)

        STOP

      endif
    enddo

    do j = 1,ncore
      if(nocore(j) == ni(i) .and. locore(j) == li(i)) then
        write(6,*)
        write(6,'("   Stopped in ",a,": valence orbital ",i3,              &
            &     " n = ",i3," l = ",i3," is also a core orbital")') cname, i, ni(i), li(i)
        write(6,*)

        STOP

      endif
    enddo

  enddo

  return

end subroutine atom_atm_check_valence
