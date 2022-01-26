!>  Prints one set of xy data to itplot and identifies it with a marker
!>
!>  \author       Sverre Froyen,  Norm Troullier, Jose Luis Martins
!>  \version      6.0.4
!>  \date         1980s and 1990s, 30 June 2021, 9 October 2021.
!>  \copyright    GNU Public License v2

subroutine atom_plot_one(nfirst, nlast, r, v, rstep, l, ist,             &
            prefix, suffix, lspd,                                        &
            itplot, mxdnr)

! ***********************************************************
! *                                                         *
! *    This is a plotting routine; the user should          *
! *  adjust for their own needs.  Prints                    *
! *  out the potential to the current plot.dat              *
! *  file (unit=itplot) for later ploting.  A marker        *
! *  is placed at the end of each group of data.            *
! *                                                         *
! ***********************************************************

! number of decimals. 9 October 2021. JLM
  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension  of the number of radial points
  integer, intent(in)               ::  itplot                           !<  default plot file

  integer, intent(in)               ::  nfirst, nlast                    !<  first and last points to plot
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid
  real(REAL64), intent(in)          ::  v(mxdnr)                         !<  potential or wave-function or...

  real(REAL64), intent(in)          ::  rstep                            !<  approximate step for plot
  integer, intent(in)               ::  l                                !<  angular momentum

  integer, intent(in)               ::  ist                              !<  sign/scaling
  logical, intent(in)               ::  lspd                             !<  l = 0,1,2,... -> s,p,d,,...

  character(len=4)                  ::  prefix                           !<  prefix of the marker
  character(len=4)                  ::  suffix                           !<  suffix of the marker


! local variables

  real(REAL64)           ::  step
  character(len=17)      ::  marker
  character(len=2)       ::  ang

! constants

  real(REAL64), parameter       ::  ZERO = 0.0_REAL64
  character(len=1), parameter  ::  IL(7) = (/'s','p','d','f','g','h','i'/)

! counter

  integer             ::  j


   if(ist == 0) return

!  Step size should be ~0.05 but it is adjustable as seen fit to give
!  a reasonalble plot.  Plot grid is coarser than computational grid.
!  Beware that for potential v(1) may be infinity, hence nfirst...

  step = ZERO
  do j = nfirst,nlast
    if (r(j) >= step) then
      write(itplot,'(1x,f18.6,3x,f20.10)') r(j), v(j)*ist
      step = r(j) + rstep
    endif
  enddo

  if(l < 0) then
    ang = "  "
  elseif(l < 7) then
    if(lspd) then
      write(ang,'(a1," ")') IL(l+1)
    else
      write(ang,'(i1," ")') l
    endif
  elseif(l < 10) then
    write(ang,'(i1," ")') l
  elseif(l < 100) then
    write(ang,'(i2)') l
  else
    write(ang,'("lx")')
  endif

  marker = 'marker '//adjustl(trim(prefix))//trim(ang)//adjustl(trim(suffix))

  write(itplot,'(1x,a17)') marker

  return

end subroutine atom_plot_one
