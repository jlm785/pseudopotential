!>  Reads from a file (plot.dat) the data for subsequent plotting
!>  of wave-functions and potentials
!>
!>  \author       Peter Schuster, Manuel Maria Alemany, Jose Luis Martins
!>  \version      6.0.5
!>  \date         4 September 2021, 21 October 2021.
!>  \copyright    GNU Public License v2

subroutine atom_plot_kb_input (iowrite, ioplot,                          &
     nmax, numw, numb, nump, nftp,                                       &
     iw, labelw, numbw, lo_w, wx, wy,                                    &
     ib, labelb, numbb, lo_b, bx, by,                                    &
     ip, labelp, numbp, lo_p, px, py,                                    &
     iq, labelq, numbq, lo_q, qx, qy,                                    &
     labelek, numbek, qek, ek,                                           &
     labelv, numbv, vr, vk)

! adapted from all-electron+pseudopotential plot code.
! kinetic, 9 October 2021. JLM
! printing. New interface. 21 October 2021. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  ioplot                           !<  tape to read input

  integer, intent(in)               ::  nmax                             !<  dimension of number of points
  integer, intent(in)               ::  numw                             !<  dimension of number of wave-functions
  integer, intent(in)               ::  numb                             !<  dimension of number of pseudo-wave-functions
  integer, intent(in)               ::  nump                             !<  dimension of number of projectors
  integer, intent(in)               ::  nftp                             !<  dimension of number of Fourier transform points

! output

  integer, intent(out)              ::  iw                               !<  number of wave-functions
  character(len=12), intent(out)    ::  labelw(numw)                     !<  label of wave-functions
  integer, intent(out)              ::  numbw(numw)                      !<  index of wave-functions
  integer, intent(out)              ::  lo_w(numw)                       !<  angular momentum of wave-function
  real(REAL64), intent(out)         ::  wx(nmax,numw), wy(nmax,numw)     !<  wave-functions

  integer, intent(out)              ::  ib                               !<  number of basis-functions
  character(len=12), intent(out)    ::  labelb(numb)                     !<  label of basis-functions
  integer, intent(out)              ::  numbb(numb)                      !<  index of basis-functions
  integer, intent(out)              ::  lo_b(numb)                       !<  angular momentum of wave-function
  real(REAL64), intent(out)         ::  bx(nmax,numb), by(nmax,numb)     !<  basis-functions

  integer, intent(out)              ::  ip                               !<  number of projectors
  character(len=12), intent(out)    ::  labelp(nump)                     !<  label of projectors
  integer, intent(out)              ::  numbp(nump)                      !<  index of projectors
  integer, intent(out)              ::  lo_p(nump)                       !<  angular momentum of projectors
  real(REAL64), intent(out)         ::  px(nmax,nump), py(nmax,nump)     !<  KB projector (real space)

  integer, intent(out)              ::  iq                               !<  number of projectors
  character(len=12), intent(out)    ::  labelq(nftp)                     !<  label of projectors
  integer, intent(out)              ::  numbq(nftp)                      !<  index of projectors
  integer, intent(out)              ::  lo_q(nftp)                       !<  angular momentum of projectors
  real(REAL64), intent(out)         ::  qx(nmax,nftp), qy(nmax,nftp)     !<  KB projectors (reciprocal space)

  character(len=12), intent(out)    ::  labelek                          !<  label of kinetic energy
  integer, intent(out)              ::  numbek                           !<  number of q points for kinetic energy
  real(REAL64), intent(out)         ::  qek(nmax), ek(nmax)              !<  kinetic energy integral

  character(len=12), intent(out)    ::  labelv                           !<  label of local potential
  integer, intent(out)              ::  numbv                            !<  number of q points for local potential
  real(REAL64), intent(out)         ::  vr(nmax), vk(nmax)               !<  local potential

! local variables

  real(REAL64), allocatable          ::  absc(:), ord(:)

  character(len=9)      ::  marker
  character(len=1)      ::  c1

  integer     ::  ndata, number, npoints, ioerr, lp1

! constants

  character(len=1), parameter   ::  IL(7) = (/'s','p','d','f','g','h','i'/)

! counters

  integer     ::  i


  allocate(absc(nmax),ord(nmax))

  iw = 0
  ib = 0
  ip = 0
  iq = 0

! Read data columns

  do ndata = 1,1000


!   (it is not yet possible to read data fixed formatted, because there
!    is not a unique format of the data in the file, unit3)

    do npoints = 1,nmax
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
      write(iowrite,'(" numbers of rows in columns of wave-functions:   ", 10i6)')   &
           (numbw(i),i=1,iw)
      write(iowrite,*)
      write(iowrite,'(" numbers of rows in columns of basis-functions: ", 10i6)')    &
           (numbb(i),i=1,ib)
      write(iowrite,*)

      exit

    else

!     end of date set
!     Backup one record and get marker

      backspace (unit=ioplot)
      read(ioplot,'(8x,a9)') marker

!     look if 'ion-change marker'

!     local potential

      if(marker(1:3) == 'loc') then

         labelv = 'Local potential'
         numbv = number-1
         do i = 1,numbv
           vr(i) = absc(i)
           vk(i) = ord(i)
         enddo

!     wave-function markers

      elseif(marker(1:3) == 'wvf') then

         iw = iw+1
         labelw(iw) = 'pseudo wf '//marker(4:4)
         numbw(iw) = number-1
         do i = 1,numbw(iw)
           wx(i,iw) = absc(i)
           wy(i,iw) = ord(i)
         enddo
         do  lp1 = 1,7
           if(marker(4:4) == IL(lp1)) lo_w(iw) = lp1-1
         enddo
         do lp1 = 1,10
           write(c1,'(i1)') lp1-1
           if(marker(4:4) == c1) lo_w(iw) = lp1-1
         enddo

!     basis-function markers

      elseif(marker(1:3) == 'bas') then

         ib = ib+1
         labelb(ib) = 'basis wf '//marker(4:4)
         numbb(ib) = number-1
         do i = 1,numbb(ib)
           bx(i,ib) = absc(i)
           by(i,ib) = ord(i)
         enddo
         do  lp1 = 1,7
           if(marker(4:4) == IL(lp1)) lo_b(ib) = lp1-1
         enddo
         do lp1 = 1,10
           write(c1,'(i1)') lp1-1
           if(marker(4:4) == c1) lo_b(ib) = lp1-1
         enddo

      elseif(marker(1:1) == 'p' .and. marker(3:3) == 'r') then

         ip = ip+1
         labelp(ip) = 'projector '//marker(2:2)
         numbp(ip) = number-1
         do i = 1,numbp(ip)
           px(i,ip) = absc(i)
           py(i,ip) = ord(i)
         enddo
         do  lp1 = 1,7
           if(marker(2:2) == IL(lp1)) lo_p(ip) = lp1-1
         enddo
         do lp1 = 1,10
           write(c1,'(i1)') lp1-1
           if(marker(2:2) == c1) lo_p(ip) = lp1-1
         enddo

      elseif(marker(1:1) == 'p' .and. marker(3:3) == 'q') then

         iq = iq+1
         labelq(iq) = 'FT proj. '//marker(2:2)
         numbq(iq) = number-1
         do i = 1,numbq(iq)
           qx(i,iq) = absc(i)
           qy(i,iq) = ord(i)
         enddo

         do  lp1 = 1,7
           if(marker(2:2) == IL(lp1)) lo_q(iq) = lp1-1
         enddo
         do lp1 = 1,10
           write(c1,'(i1)') lp1-1
           if(marker(2:2) == c1) lo_q(iq) = lp1-1
         enddo

      elseif(marker(1:3) == 'eki') then

         labelek = 'Integral E_k'
         numbek = number-1
         do i = 1,numbek
           qek(i) = absc(i)
           ek(i) = ord(i)
         enddo

      else

         write(iowrite,*)
         write(iowrite,*)
         write(iowrite,*) ' Skipping data in atom_plot_kb_input.  Unknown marker'
         write(iowrite,*) ' check plot file (default plot_kb.dat)'
         write(iowrite,*) ' offending marker is :',marker
         write(iowrite,*)

      endif

    endif

  enddo

  return

end subroutine atom_plot_kb_input
