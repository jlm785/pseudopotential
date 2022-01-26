!>  Reads from a file (plot_kb.dat) the data for subsequent plotting
!>  and figures out the required array dimensions
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.2
!>  \date         11 August 2021
!>  \copyright    GNU Public License v2

subroutine atom_plot_kb_input_size(ioplot, nmax, numw, numb, nump, nftp)


! adapted from all-electron+pseudopotential plot code.



  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  ioplot                           !<  tape to read input

! output

  integer, intent(out)              ::  nmax                             !<  dimension of number of points
  integer, intent(out)              ::  numw                             !<  dimension of number of wave-functions
  integer, intent(out)              ::  numb                             !<  dimension of number of basis-functions
  integer, intent(out)              ::  nump                             !<  dimension of number of projectors
  integer, intent(out)              ::  nftp                             !<  dimension of number of Fourier transform points

! local variables

  real(REAL64)          ::  absc, ord

  character(len=3)      ::  marker

  integer     ::  ndata, number, npoints, ioerr

  integer     ::  iw, ib, ip, iq

  nmax = 0

  numw = 0
  numb = 0
  nump = 0

  nftp = 0

  iw = 0
  ib = 0
  ip = 0
  iq = 0

! Read data columns max iterations should be more than enough

  do ndata = 1,1000


!   (it is not yet possible to read data fixed formatted, because there
!    is not a unique format of the data in the file, unit3)

    do npoints = 1,10000
       number = npoints
       read(ioplot,*,iostat = ioerr) absc, ord

       if(ioerr /= 0) exit

    enddo

    if(ioerr == 0) then
            write(6,*)
            write(6,*)' set new max dimension in atom_plot_input_size', npoints
            write(6,*)

            STOP

    elseif(ioerr < 0) then

!     arrived at the end of file

      exit

    else

!     end of date set
!     Backup one record and get marker

      backspace (unit=ioplot)
      read(ioplot,'(8x,a3)') marker

!     local potential

      if(marker(1:3) == 'loc') then

         if(number > nmax) nmax = number

!     wave-function

      elseif(marker(1:3) == 'wvf') then

         iw = iw+1
         if(number > nmax) nmax = number

!     fourier markers - pseudopotentials

      elseif(marker(1:3) == 'bas') then

         ib = ib+1
         if(number > nmax) nmax = number

      elseif(marker(1:1) == 'p' .and. marker(3:3) == 'r') then

        ip = ip+1
        if(number > nmax) nmax = number

      elseif(marker(1:1) == 'p' .and. marker(3:3) == 'q') then

        iq = iq+1
        if(number > nmax) nmax = number

      elseif(marker(1:3) == 'eki') then

        if(number > nmax) nmax = number

     else

        write(6,*)
        write(6,*)
        write(6,*) ' Skipping data in atom_plot_kb_input_size.  Unknown marker'
        write(6,*) ' check plot file (default plot_kb.dat)'
        write(6,*) ' offending marker is :',marker
        write(6,*)

      endif

    endif

  enddo

  numw = iw
  numb = ib
  nump = ip
  nftp = iq

  return

end subroutine atom_plot_kb_input_size
