!>  Modifies the parameters for the non-local KB pseudopotential.
!>
!>  \author       J.L.Martins
!>  \version      6.1.0
!>  \date         14 January 2026.
!>  \copyright    GNU Public License v2

subroutine atom_kb_setup_pseudo(llocal, lmax_pot, nql, delql)

! split from atom_kb_setup. 14 January 2026. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer,intent(in)                ::  lmax_pot                         !<  maximum angular momentum for potential

! input and output

  integer, intent(inout)            ::  llocal                           !<  angular momentum for local potential (negative: maximum of l-dependent)
  integer, intent(inout)            ::  nql                              !<  number of points for the Fourier grid of local potential and densities
  real(REAL64), intent(inout)       ::  delql                            !<  spacing of points in Fourier grid

! local variables

  character(len=1)        ::  yesno
  integer                 ::  n_in
  real(REAL64)            ::  r_in

! constants

  real(REAL64), parameter           ::  ZERO = 0.0_REAL64


  write(6,*)
  write(6,*) '  The current parameters for the Kleiman-Bylander'
  write(6,*) '  non-local pseudopotential transformation are:'
  write(6,*)
  if(llocal < 0) then
    write(6,*) '  The maximum of all potentials will be used as local'
  else
    write(6,*) '  The local potential will be  l = ',llocal
  endif
  write(6,*)
  write(6,'("   The Fourier grid will range from q = 0.0000 to",         &
       &  f10.4, "  with steps of ",f10.4)') nql*delql, delql
  write(6,'("   Enough for a cutoff of ",f10.3," Ha.")')                 &
             (nql*delql)*(nql*delql)/8
  write(6,*)


  write(6,*)
  write(6,*) '  Do you want to change any of those parameters? (y/n)'
  write(6,*)

  read(5,*) yesno

  if(yesno == 'Y' .or. yesno == 'y') then

    write(6,*)
    write(6,*) '  Do you want to change the local potential? (y/n)'
    write(6,*)

    read(5,*) yesno

    if(yesno == 'Y' .or. yesno == 'y') then

      write(6,*)
      write(6,*) '  Enter the angular momentum for local potential'
      write(6,*) '  between 0 and ',lmax_pot,'.  Enter negative for maximum'
      write(6,*)

      read(5,*) n_in

      if(n_in > lmax_pot) then
        llocal = lmax_pot
        write(6,*)
        write(6,*) '  Wrong value. Will use l = ',llocal
        write(6,*)
      else
        llocal = n_in
      endif

    endif

    write(6,*)
    write(6,*) '  Do you want to change the Fourier grid? (y/n)'
    write(6,*)

    read(5,*) yesno

    if(yesno == 'Y' .or. yesno == 'y') then

      write(6,*)
      write(6,*) '  Enter the number of steps and step value (nql, delql)'
      write(6,*)

      read(5,*) n_in, r_in

      if(n_in < 0) then
        write(6,*)
        write(6,*) '  Wrong value. Will use nql = ',nql
        write(6,*)
      else
        if(n_in < 100 .or. n_in > 50000) write(6,*) '  Strange value of nql...'
        nql = n_in
      endif

      if(r_in < ZERO) then
        write(6,*)
        write(6,*) '  Wrong value. Will use delql = ',delql
        write(6,*)
      else
        if(r_in < 0.0001 .or. r_in > 0.5) write(6,*) '  Strange value of delql...'
        delql = r_in
      endif
    endif

  else

    write(6,*)
    write(6,*) '  Continuing with the above parameters'
    write(6,*)

  endif


  return

end subroutine atom_kb_setup_pseudo
