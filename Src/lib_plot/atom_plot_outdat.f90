!>  Writes the files with data for subsequent plotting
!>  of wave-functions and potentials with gnuplot
!>
!>  \author       Peter Schuster, Manuel Maria Alemany, Jose Luis Martins
!>  \version      6.013
!>  \date         April 1993, January 2000, 6 August 2021, 19 October 2021.
!>  \copyright    GNU Public License v2

subroutine atom_plot_outdat (iost, nmax, numw, nump, nftp,               &
      wx, wy, py, vx, vy, fx, fy, fwx, fwy,                              &
      numbw, numbv, numbf, numbfw,                                       &
      iw, ip, iv, ifp, ifw,                                              &
      labelw, labelp, labelv, labelf, labelfw, yminwp, ymaxwp, yminv)


! adapted from old code
! originally written by Peter Schuster, April 1993
! interative ploting added by Manuel Maria Alemany, January 2000
! removed H formats 22 January 2008. JLM
! converted to fortran 90, 6 August 2021
! corrected gnuplot comment statements. 19 October 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  iost                             !<  files use tapes iost+1 to iost+4

  integer, intent(in)               ::  nmax                             !<  dimension of number of points
  integer, intent(in)               ::  numw                             !<  dimension of number of true wave-functions
  integer, intent(in)               ::  nump                             !<  dimension of number of pseudo-wave-functions
  integer, intent(in)               ::  nftp                             !<  dimension of number of Fourier transform points

  integer, intent(inout)            ::  iw                               !<  number of true wave-functions
  character(len=12), intent(in)     ::  labelw(numw)                     !<  label of true wave-functions
  integer, intent(in)               ::  numbw(numw)                      !<  index of true wave-functions
  real(REAL64), intent(in)          ::  wx(nmax,numw), wy(nmax,numw)     !<  true wave-functions

  integer, intent(in)               ::  ip                               !<  number of pseudo-wave-functions
  character(len=12), intent(in)     ::  labelp(numw)                     !<  label of pseudo-wave-functions
!  integer, intent(in)               ::  numbp(numw)                      !<  index of pseudo-wave-functions
!  real(REAL64), intent(in)          ::  px(nmax,numw), py(nmax,numw)     !<  pseudo wave-functions
  real(REAL64), intent(in)          ::  py(nmax,numw)                    !<  pseudo wave-functions

  integer, intent(in)               ::  iv                               !<  number of pseudo-potentials
  character(len=12), intent(in)     ::  labelv(nump)                     !<  label of pseudo-potentials
  integer, intent(in)               ::  numbv(nump)                      !<  index of pseudo-potentials
  real(REAL64), intent(in)          ::  vx(nmax,nump), vy(nmax,nump)     !<  pseudo-potentials

  integer, intent(in)               ::  ifp                              !<  number of pseudo-potentials transforms
  character(len=12), intent(in)     ::  labelf(nump)                     !<  label of pseudo-potentials transforms
  integer , intent(in)              ::  numbf(nump)                      !<  index of pseudopotentials transforms
  real(REAL64), intent(in)          ::  fx(nftp,nump), fy(nftp,nump)     !<  pseudo-potentials transforms

  integer , intent(in)              ::  ifw                              !<  number of pseudo-wave-function transforms
  character(len=12), intent(in)     ::  labelfw(nump)                    !<  label of pseudo-wave-function transforms
  integer , intent(in)              ::  numbfw(nump)                     !<  index of pseudo-wave-function transforms
  real(REAL64), intent(in)          ::  fwx(nftp,nump), fwy(nftp,nump)   !<  pseudo-wave-function transforms

! output

  real(REAL64), intent(out)         ::  yminwp, ymaxwp, yminv            !<  min and max

! local variables

  integer     ::  iowf, iofw, iopp, iofp

! parameters

  real(REAL64), parameter                ::  ZERO = 0.0_REAL64

! counters

  integer     ::  i, j, k, kc


! opening gnuplot files to write to

!  n+1: 'wfct.gp'.....ae and pseudo wave-functions
!  n+2: 'fwfct.gp'....ae and pseudo wave-functions + FT's(pwf's)
!  n+3: 'ppot.gp'.....pseudo potentials
!  n+4: 'fppot.gp'....pseudo potentials + FT's(pp's)

  iowf = iost+1
  iofw = iost+2
  iopp = iost+3
  iofp = iost+4

  open(unit = iowf, file = 'wfct.gp', status = 'unknown', form = 'formatted')
  open(unit = iofw, file = 'fwft.gp', status = 'unknown', form = 'formatted')
  open(unit = iopp, file = 'ppot.gp', status = 'unknown', form = 'formatted')
  open(unit = iofp, file = 'fppt.gp', status = 'unknown', form = 'formatted')

!     write true wave-functions and pseudo wafe-functions to 'wfct.gp'
!     output format:
!                    x; wst;wsp;(wst;wsp);...;wgp;wgt;(wgt;wgp)

  if(iw /= 0 .or. ip /= 0) then
    if(iw < ip) iw = ip
    write(iowf,'("#     r (a.u.)",3x,20(1x,a12))') (labelw(i),labelp(i),i=1,iw)
    write(iowf,'("#")')
    yminwp = ZERO
    ymaxwp = ZERO
    do i = 1,numbw(iw)
       write(iowf,'(20(3x,f10.6))') wx(i,1),(wy(i,kc),py(i,kc),kc=1,iw)
       do k = 1,iw
          if(wy(i,k) < yminwp) yminwp = wy(i,k)
          if(wy(i,k) > ymaxwp) ymaxwp = wy(i,k)
          if(py(i,k) < yminwp) yminwp = py(i,k)
          if(py(i,k) > ymaxwp) ymaxwp = py(i,k)
       enddo
    enddo

  endif

! write FT(pseudo wave-function)                     to 'fwft.gp'
!  output format:   x; fw0;(fw0);...;fw4;(fw4)


  if(ifw == 0) then
     write(iofw,*)'# NO NUMBERS'
  else
     write(iofw,'("#   q (1/a.u.)",3x,10(1x,a12))') (labelfw(i),i=1,ifw)
     write(iofw,'("#")')
     do i = 1,numbfw(ifw)
        write(iofw,'(20(3x,f10.6))') fwx(i,1),(fwy(i,kc),kc=1,ifw)
     enddo
  endif


! write ionic pseudopotentials                        to 'ppot.gp'
! output format:
!                x; vn0;(vn0);...;vn4;(vn4)


!  continue
  if(iv == 0) then
     write(iopp,*)'# NO NUMBERS'
  else
     yminv = 0.0
     write(iopp,'("#      r (a.u.)",1x,10(3x,a12))') (labelv(j),j=1,iv)
     write(iopp,'("#")')
     do i = 1,numbv(iv)
        write(iopp,'(11(4x,f11.6))') vx(i,1),(vy(i,kc),kc=1,iv)
       do k = 1,iv
           if(vy(i,k) < yminv) yminv = vy(i,k)
        enddo
     enddo
  endif


! write FT(ionic pseudopotentials)                   to 'fppot.gp'
! output format:
!                x; fn1;(fn1);...;fn5;(fn5)

  if(ifp == 0)then
     write(iofp,*)'# NO NUMBERS'
  else
     write(iofp,'("#    q (1/a.u.)",1x,10(2x,a12))') (labelf(i),i=1,ifp)
     write(iofp,'("#")')
     do i = 1,numbf(ifp)
        write(iofp,'(20(4x,f10.6))') fx(i,1),(fy(i,kc),kc=1,ifp)
     enddo
  endif


  close(iowf)
  close(iofw)
  close(iopp)
  close(iofp)

  return

end subroutine atom_plot_outdat
