!>  Asks the user for changes in core radii
!>
!>  \author       Jose Luis Martins
!>  \version      6.013
!>  \date         30 August June 2021
!>  \copyright    GNU Public License v2

subroutine atom_psd_change_radii(lpmax, rc, rc_tab, rz, rx)

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter                ::  lc = 4                           !  maximum valence angular momentum

! input

  real(REAL64), intent(in)          ::  rc_tab(0:lc)                     !<  tabulated core radii
  real(REAL64), intent(in)          ::  rz(0:lc), rx(0:lc)               !<  avearge zero and extrma for l

! modified

  integer, intent(inout)            ::  lpmax                            !<  maximum value of l
  real(REAL64), intent(inout)       ::  rc(0:lc)                         !<  core radius r_c(l)

! local variables

  character(len=1)                  ::  yesno
  integer                           ::  l_in
  real(REAL64)                      ::  rc_in
  integer                           ::  lpmax_old
  integer                           ::  lpc

! counter

  integer        ::  lp


  lpmax_old = lpmax

  write(6,*)
  write(6,*)
  write(6,*) '  The current value of core radii are:'
  write(6,*)
  do lp = 0,lpmax
    write(6,'("   l = ",i2,"  rc = ",f12.3)') lp, rc(lp)
  enddo
  write(6,*)
  write(6,*) '  Do you want ot modify them? (y/n)'
  write(6,*)
  read(5,*) yesno

  if(yesno == 'y' .or. yesno == 'Y') then

    write(6,*)
    write(6,*)  '  l =     current    tabulated   zero of psi  max of psi'
    write(6,*)
    do lp = 0,lpmax
      write(6,'(i5,4(2x,f10.3))') lp, rc(lp), rc_tab(lp), rz(lp), rx(lp)
    enddo

    write(6,*)
    write(6,*) '  Enter new value of maximum angular momentum 0,...,',lc
    write(6,*) '  Current value is: ', lpmax
    write(6,*)
    read(5,*) l_in
    if(l_in >= 0 .and. l_in <= lc) then
      lpmax = l_in
      do lp = 0,lpmax
        write(6,*)
        write(6,*) '  Enter new core radius for angular momentum l = ',lp
        if(lp > lpmax_old) then
          write(6,*) '  Validity of rc not checked (new l)'
          lpc = lpmax_old
        else
          write(6,'("   Current value is: ",f12.3)') rc(lp)
          lpc = lp
        endif
        write(6,*)
        read(5,*) rc_in
        if(rc_in > rz(lp)) then
          if(rc_in < 1.2*rz(lpc) .or. rc_in > 4*max(rz(lpc),rx(lpc))) then
            write(6,*)
            write(6,*) '  Value looks suspicious.'
            write(6,*) '  Do you want to keep it? (y/n)'
            write(6,*)
            read(5,*) yesno
            if(yesno == 'y' .or. yesno == 'Y') then
              rc(lp) = rc_in
            else
              write(6,*)
              write(6,*) '  Enter new value (not checked)'
              write(6,*)
              read(5,*) rc(lp)
            endif
          else
            rc(lp) = rc_in
          endif
        else
          write(6,*)
          write(6,*) '  Wrong value proceeding with original value'
          write(6,*)
        endif
      enddo
    else
      write(6,*)
      write(6,*) '  Wrong answer proceeding with the value above'
      write(6,*)
    endif



  elseif(yesno == 'n' .or. yesno == 'N') then
  else
    write(6,*)
    write(6,*) '  Wrong answer proceeding with the above values'
    write(6,*)
  endif


  return

end subroutine atom_psd_change_radii
