!>  Reads from a file (plot.dat) the data for subsequent plotting
!>  of wave-functions and potentials
!>
!>  \author       Peter Schuster, Manuel Maria Alemany, Jose Luis Martins
!>  \version      6.0.5
!>  \date         April 1993, January 2000, 6 August 2021, 21 october 2021.
!>  \copyright    GNU Public License v2

subroutine atom_plot_input (ioplot, nmax, numw, nump, nftp,              &
     & wx, wy, px, py, vx, vy, fx, fy, fwx, fwy,                         &
     & numbw, numbp, numbv, numbf, numbfw,                               &
     & iw, ip, iv, ifp, ifw, labelw,                                     &
     & labelp, labelv, labelf, labelfw, zion,                            &
     & iowrite)

! adapted from old code
! originally written by Peter Schuster, April 1993
! interative ploting added by Manuel Maria Alemany, January 2000
! removed H formats 22 January 2008. JLM
! converted to fortran 90, 6 August 2021
! nbig, 18 october 2021. JLM
! printing. New interface. 21 October 2021. JLM



  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  ioplot                           !<  tape to read input

  integer, intent(in)               ::  nmax                             !<  dimension of number of points
  integer, intent(in)               ::  numw                             !<  dimension of number of wave-functions
  integer, intent(in)               ::  nump                             !<  dimension of number of pseudo-wave-functions
  integer, intent(in)               ::  nftp                             !<  dimension of number of Fourier transform points

! output

  integer, intent(out)              ::  iw                               !<  number of true wave-functions
  character(len=12), intent(out)    ::  labelw(numw)                     !<  label of true wave-functions
  integer, intent(out)              ::  numbw(numw)                      !<  index of true wave-functions
  real(REAL64), intent(out)         ::  wx(nmax,numw), wy(nmax,numw)     !<  true wave-functions

  integer, intent(out)              ::  ip                               !<  number of pseudo-wave-functions
  character(len=12), intent(out)    ::  labelp(numw)                     !<  label of pseudo-wave-functions
  integer, intent(out)              ::  numbp(numw)                      !<  index of pseudo-wave-functions
  real(REAL64), intent(out)         ::  px(nmax,numw), py(nmax,numw)     !<  pseudo wave-functions

  integer, intent(out)              ::  iv                               !<  number of pseudo-potentials
  character(len=12), intent(out)    ::  labelv(nump)                     !<  label of pseudo-potentials
  integer, intent(out)              ::  numbv(nump)                      !<  index of pseudopotentials
  real(REAL64), intent(out)         ::  vx(nmax,nump), vy(nmax,nump)     !<  pseudo-potentials

  integer, intent(out)              ::  ifp                              !<  number of pseudo-potentials transforms
  character(len=12), intent(out)    ::  labelf(nump)                     !<  label of pseudo-potentials transforms
  integer , intent(out)             ::  numbf(nump)                      !<  index of pseudopotentials transforms
  real(REAL64), intent(out)         ::  fx(nftp,nump), fy(nftp,nump)     !<  pseudo-potentials transforms

  integer , intent(out)             ::  ifw                              !<  number of pseudo-wave-function transforms
  character(len=12), intent(out)    ::  labelfw(nump)                     !<  label of pseudo-wave-function transforms
  integer , intent(out)             ::  numbfw(nump)                     !<  index of pseudo-wave-function transforms
  real(REAL64), intent(out)         ::  fwx(nftp,nump), fwy(nftp,nump)   !<  pseudo-wave-function transforms

  real(REAL64), intent(out)         ::  zion                             !<  ionic charge

! local variables

  real(REAL64), allocatable          ::  absc(:), ord(:)

  character(len=3)      ::  marker

  integer     ::  ndata, number, npoints, ioerr, lp1
  integer     ::  nbig

! constants

  character(len=1), parameter   ::  IL(7) = (/'s','p','d','f','g','h','i'/)

! counters

  integer     ::  i


  nbig = max(nmax,nftp)
  allocate(absc(nbig),ord(nbig))

  iw = 0
  ip = 0
  iv = 0
  ifp = 0
  ifw = 0


! Read data columns

  do ndata = 1,1000

!   (it is not yet possible to read data fixed formatted, because there
!    is not a unique format of the data in the file, unit3)

    do npoints = 1,nbig
       number = npoints
       read(ioplot,*,iostat = ioerr) absc(npoints), ord(npoints)

       if(ioerr /= 0) exit

    enddo

    if(ioerr == 0) then

            write(6,*)
            write(6,*)' NOT ALL NUMBERS READ! SET NEW DIMENSION FOR NMAX!'
            write(6,*)

            STOP

    elseif(ioerr < 0) then

!     arrived at the end of file
!     information for program-user

      write(iowrite,*)
      write(iowrite,'(" numbers of rows in columns of true wave-functions:   ", 10i6)')   &
           (numbw(i),i=1,iw)
      write(iowrite,*)
      write(iowrite,'(" numbers of rows in columns of pseudo wave-functions: ", 10i6)')   &
           (numbp(i),i=1,ip)
      write(iowrite,*)
      write(iowrite,'(" numbers of rows in columns of pseudopotentials:      ", 10i6)')   &
           (numbv(i),i=1,iv)
      write(iowrite,*)
      write(iowrite,'(" numbers of rows in columns of FT - ppot:             ", 10i6)')   &
           (numbf(i),i=1,ifp)
      write(iowrite,*)
      write(iowrite,'(" numbers of rows in columns of FT - wave-functions:   ", 10i6)')   &
           (numbfw(i),i=1,ifw)
      write(iowrite,*)

      exit

    else

!     end of date set
!     Backup one record and get marker

      backspace (unit=ioplot)
      read(ioplot,'(8x,a3)') marker

!     look if 'ion-change marker'

      if(marker == 'zio') then
         read(ioplot,'(2x,f5.2)') zion

!     true wave-function markers

      elseif(marker(1:1) == 'w' .and. marker(3:3) == 't')then

         iw = iw+1
         labelw(iw) = 'true wf '//marker(2:2)
         numbw(iw) = number-1
         do i = 1,numbw(iw)
           wx(i,iw) = absc(i)
           wy(i,iw) = ord(i)
         enddo
!         call atom_plot_field(nmax,numw,numbw,absc,ord,wx,wy,iw)

!     pseudo wave-function markers

      elseif(marker(1:1) == 'w' .and. marker(3:3) == 'p')then

         ip = ip+1
         labelp(ip) = 'pseudo wf '//marker(2:2)
         numbp(ip) = number-1
         do i = 1,numbp(ip)
           px(i,ip) = absc(i)
           py(i,ip) = ord(i)
         enddo
!         call atom_plot_field(nmax,numw,numbp,absc,ord,px,py,ip)

!     potential markers

      elseif(marker(1:2) == 'vn')then
         iv = iv+1
         labelv(iv) = marker(3:3)//' - nonlocal'
         numbv(iv) = number-1
         do i = 1,numbv(iv)
           vx(i,iv) = absc(i)
           vy(i,iv) = ord(i)
         enddo
!         call atom_plot_field(nmax,nump,numbv,absc,ord,vx,vy,iv)

!     fourier markers - pseudopotentials

      elseif(marker(1:2) == 'fn')then

         ifp = ifp+1
         read(marker(3:3),'(i1)') lp1
         labelf(ifp) = IL(lp1)//' - nonlocal'
         numbf(ifp) = number-1
         do i = 1,numbf(ifp)
           fx(i,ifp) = absc(i)
           fy(i,ifp) = ord(i)
         enddo
!         call atom_plot_field(nftp,nump,numbf,absc,ord,fx,fy,ifp)

!     fourier markers - pseudo wave-functions

      elseif(marker(1:2) == 'fw')then
         ifw = ifw+1
         read(marker(3:3),'(i1)') lp1
         labelfw(ifw) = 'FT-wf '//IL(lp1+1)
         numbfw(ifw) = number-1
         do i = 1,numbfw(ifw)
           fwx(i,ifw) = absc(i)
           fwy(i,ifw) = ord(i)
         enddo
!         call atom_plot_field(nmax,nump,numbfw,absc,ord,fwx,fwy,ifw)

      else

        write(6,*)
        write(6,*)
        write(6,*) ' Stopped in atom_plot_input.  Error in markers'
        write(6,*) ' check plot file (default plot.dat)'
        write(6,*) ' offending marker is :',marker
        write(6,*)

         STOP

      endif

    endif

  enddo

  deallocate(absc,ord)

  return

end subroutine atom_plot_input
