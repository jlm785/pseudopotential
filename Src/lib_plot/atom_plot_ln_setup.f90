!>  Sets up the files and radii for the plots a comparison of the
!>  logarithmic derivative ( d psi / dr ) / psi
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.1.4
!>  \date         30 August 2021
!>  \copyright    GNU Public License v2

subroutine atom_plot_ln_setup(datafile, pseudo, pseudokb, lkb, rpoint)


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! output

  character(len=*), intent(out)     ::  datafile                         !<  file name for all electron results
  character(len=*), intent(out)     ::  pseudo                           !<  file name for semi-local pseudopotential results
  character(len=*), intent(out)     ::  pseudokb                         !<  file name for non-local KB pseudopotential results

  logical, intent(out)              ::  lkb                              !<  If KB results are included in the plots
  real(REAL64), intent(out)         ::  rpoint                           !<  radius where log derivatives are calculated

! local variables

  character(len=2)                  ::  kbtype

  write(6,*)
  lkb = .FALSE.
  write(6,*)' Do you want to test also the  Kleinman-Bylander,'
  write(6,*)' form of the pseudopotential?'
  write(6,*)' Enter kb for yes'
  read(5,'(a2)') kbtype
  if(kbtype == 'kb' .or. kbtype == 'KB') lkb = .TRUE.


  write(6,*)
  write(6,*)' Enter data file with all-electron information (max 15 chr.):'
  write(6,*)'(default value = datafile.dat)'
!  write(6,*)'(atom.exe < atom.dat !!!WITH MIN OPTION!!!'
!  write(6,*)'          > datafile.dat)'

  read(5,'(a15)') datafile
  if (datafile == ' ') datafile = 'datafile.dat   '
  write(6,*)


  write(6,*)
  write(6,*)' Enter semi-local pseudopotential file (max 15 chr.):'
  write(6,*)'(default value = pseudo.dat01)'

  read(5,'(a15)') pseudo
  if (pseudo == ' ') THEN
    pseudo = 'pseudo.dat     '
!    pseudo='pseudo.dat     '
  endif
  write(6,*)

  if(lkb) then
    write(6,*)
    write(6,*)' Enter pseudopotential file (max 15 chr.):'
    write(6,*)'(default value = pseudokb.dat)'

    read(5,'(a15)') pseudokb
    if (pseudokb == ' ') THEN
      pseudokb = 'pseudokb.dat   '
    endif
    write(6,*)
  endif

  write(6,*)' Enter distance (in a.u.) at which the'
  write(6,*)' logaritmic derivatives are calculated'
  read(5,*) rpoint
  write(6,*)
  write(6,*)

  return

end subroutine atom_plot_ln_setup
