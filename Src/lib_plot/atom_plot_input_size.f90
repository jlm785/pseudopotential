!>  Reads from a file (plot.dat) the data for subsequent plotting
!>  and figures out the required array dimensions
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.4
!>  \date         11 August 2021, 9 October 2021.
!>  \copyright    GNU Public License v2

subroutine atom_plot_input_size(ioplot, nmax, numw, nump, nftp)

! adapted from old code
! originally written by Peter Schuster, April 1993
! interative ploting added by Manuel Maria Alemany, January 2000
! removed H formats 22 January 2008. JLM
! converted to fortran 90, 6 August 2021
! kinetic, 9 October 2021. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  ioplot                           !<  tape to read input

! output

  integer, intent(out)              ::  nmax                             !<  dimension of number of points
  integer, intent(out)              ::  numw                             !<  dimension of number of wave-functions
  integer, intent(out)              ::  nump                             !<  dimension of number of pseudo-wave-functions
  integer, intent(out)              ::  nftp                             !<  dimension of number of Fourier transform points

! local variables

  real(REAL64)          ::  absc, ord

  character(len=3)      ::  marker

  integer     ::  ndata, number, npoints, ioerr

  integer     ::  iw, ip, iv, ifp, ifw

  nmax = 0
  numw = 0
  nump = 0
  nftp = 0

  iw = 0
  ip = 0
  iv = 0
  ifp = 0
  ifw = 0


! Read data columns

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

!     look if 'ion-change marker'

      if(marker == 'zio') then
         read(ioplot,*) absc


!     true wave-function markers

      elseif(marker(1:1) == 'w' .and. marker(3:3) == 't')then

         iw = iw+1
         if(number > nmax) nmax = number

!     pseudo wave-function markers

      elseif(marker(1:1) == 'w' .and. marker(3:3) == 'p')then

         ip = ip+1
         if(number > nmax) nmax = number

!     potential markers

      elseif(marker(1:2) == 'vn')then

         iv = iv+1
         if(number > nmax) nmax = number

!     fourier markers - pseudopotentials

      elseif(marker(1:2) == 'fn')then

         ifp = ifp+1
         if(number > nftp) nftp = number

!     fourier markers - pseudo wave-functions

      elseif(marker(1:2) == 'fw')then

         ifw = ifw+1
         if(number > nftp) nftp = number

      else

        write(6,*)
        write(6,*)
        write(6,*) ' Stopped in atom_plot_input_size.  Error in markers'
        write(6,*) ' check plot file (default plot.dat)'
        write(6,*) ' offending marker is :',marker
        write(6,*)


         STOP

      endif

    endif

  enddo

  numw = max(iw,ip)
  nump = max(iv,ifp,ifw)

  return

end subroutine atom_plot_input_size
