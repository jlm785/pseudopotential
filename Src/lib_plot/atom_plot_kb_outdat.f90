!>  Writes the files with data for subsequent plotting
!>  of wave-functions and potentials with gnuplot
!>
!>  \author       Peter Schuster, Manuel Maria Alemany, Jose Luis Martins
!>  \version      6.0.5
!>  \date         4 September 2021, 9 October 2021.
!>  \copyright    GNU Public License v2

subroutine atom_plot_kb_outdat (iost, nmax, numw, numb, nump, nftp,      &
     iw, labelw, numbw, lo_w, wx, wy,                                    &
     ib, labelb, numbb, lo_b, bx, by,                                    &
     ip, labelp, numbp,       px, py,                                    &
     iq, labelq, numbq,       qx, qy,                                    &
     labelek, numbek, qek, ek,                                           &
     labelv, numbv, vr, vk,                                              &
     yminwp, ymaxwp, yminp, ymaxp, yminq, ymaxq,                         &
     xmaxek, yminek, ymaxek, yminv, ymaxv)

! adapted from all-electron+pseudopotential plot code.
! kinetic, 9 October 2021. JLM
! nump can be zero. 18 October 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  iost                             !<  files use tapes iost+1 to iost+4

  integer, intent(in)               ::  nmax                             !<  dimension of number of points
  integer, intent(in)               ::  numw                             !<  dimension of number of wave-functions
  integer, intent(in)               ::  numb                             !<  dimension of number of basis-functions
  integer, intent(in)               ::  nump                             !<  dimension of number of projectors
  integer, intent(in)               ::  nftp                             !<  dimension of number of Fourier transform points

  integer, intent(inout)            ::  iw                               !<  number of wave-functions
  character(len=12), intent(in)     ::  labelw(numw)                     !<  label of wave-functions
  integer, intent(in)               ::  numbw(numw)                      !<  index of wave-functions
  integer, intent(in)               ::  lo_w(numw)                       !<  angular momentum of wave-function
  real(REAL64), intent(in)          ::  wx(nmax,numw), wy(nmax,numw)     !<  wave-functions

  integer, intent(in)               ::  ib                               !<  number of basis-functions
  character(len=12), intent(in)     ::  labelb(numb)                     !<  label of basis-functions
  integer, intent(in)               ::  numbb(numb)                      !<  index of basis-functions
  integer, intent(in)               ::  lo_b(numb)                       !<  angular momentum of wave-function
  real(REAL64), intent(in)          ::  bx(nmax,numb), by(nmax,numb)     !<  basis-functions

  integer, intent(in)               ::  ip                               !<  number of projectors
  character(len=12), intent(in)     ::  labelp(nump)                     !<  label of projectors
  integer, intent(in)               ::  numbp(nump)                      !<  index of projectors
  real(REAL64), intent(in)          ::  px(nmax,nump), py(nmax,nump)     !<  KB projectors (real space)

  integer, intent(in)               ::  iq                               !<  number of projectors
  character(len=12), intent(in)     ::  labelq(nftp)                     !<  label of projectors
  integer, intent(in)               ::  numbq(nftp)                      !<  index of projectors
  real(REAL64), intent(in)          ::  qx(nmax,nftp), qy(nmax,nftp)     !<  KB projectors (reciprocal space)

  character(len=12), intent(in)     ::  labelek                          !<  label of kinetic energy
  integer, intent(in)               ::  numbek                           !<  number of q points for kinetic energy
  real(REAL64), intent(in)          ::  qek(nmax), ek(nmax)              !<  kinetic energy integral

  character(len=12), intent(in)     ::  labelv                           !<  label of local potential
  integer, intent(in)               ::  numbv                            !<  number of q points for local potential
  real(REAL64), intent(in)          ::  vr(nmax), vk(nmax)               !<  local potential

! output

  real(REAL64), intent(out)         ::  yminwp, ymaxwp                   !<  min and max
  real(REAL64), intent(out)         ::  yminp, ymaxp                     !<  min and max
  real(REAL64), intent(out)         ::  yminq, ymaxq                     !<  min and max
  real(REAL64), intent(out)         ::  xmaxek, yminek, ymaxek           !<  min and max
  real(REAL64), intent(out)         ::  yminv, ymaxv                     !<  min and max

! local variables

  integer, allocatable              ::  ipoint_b(:)

  integer     ::  iob, lmax
  integer     ::  nbas, ipoint_w

  character(len=13)     ::  filename

! parameters

  real(REAL64), parameter                ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64
  character(len=1), parameter   ::  IL(7) = (/'s','p','d','f','g','h','i'/)

! counters

  integer     ::  i, kc, l


! opening gnuplot files to write to


  iob = iost+1

  lmax = 0
  do i = 1,numw
    if(lmax < lo_w(i)) lmax = lo_w(i)
  enddo
  do i = 1,numb
    if(lmax < lo_b(i)) lmax = lo_b(i)
  enddo

  yminwp = ZERO
  ymaxwp = ZERO

  yminp = ZERO
  ymaxp = ZERO

  yminq = ZERO
  ymaxq = ZERO

  yminv = ZERO
  ymaxv = ZERO

! loop over angular momentum of basis functions

  do l = 0,lmax

    write(filename,'("bas",i1,".gp")') l
    open(unit = iob, file = filename, status = 'unknown', form = 'formatted')
    ipoint_w = 0
    do i = 1,iw
      if(lo_w(i) == l) ipoint_w = i
    enddo

    nbas = 0
    do i = 1,ib
      if(lo_b(i) == l) nbas = nbas + 1
    enddo
    allocate(ipoint_b(nbas))
    nbas = 0
    do i = 1,ib
      if(lo_b(i) == l) then
        nbas = nbas + 1
        ipoint_b(nbas) = i
      endif
    enddo
    if(ipoint_w == 0) then

      write(iob,'("#     r (a.u.)",3x,20(1x,a12))') (labelb(ipoint_b(i)),i=1,nbas)
      write(iob,'("#")')
      do i = 1,numbb(ipoint_b(1))
        write(iob,'(20(3x,f10.6))') bx(i,1),(by(i,ipoint_b(kc)),kc=1,nbas)
        do kc = 1,nbas
          if(by(i,kc) < yminwp) yminwp = by(i,kc)
          if(by(i,kc) > ymaxwp) ymaxwp = by(i,kc)
        enddo
      enddo

    else

      write(iob,'("#     r (a.u.)",3x,20(1x,a12))') labelw(ipoint_w),    &
                 (labelb(ipoint_b(i)),i=1,nbas)
      write(iob,'("#")')
      do i = 1,numbw(ipoint_w)
        write(iob,'(20(3x,f10.6))') wx(i,1),wy(i,ipoint_w),              &
                 (by(i,ipoint_b(kc)),kc=1,nbas)
        if(wy(i,1) < yminwp) yminwp = wy(i,1)
        if(wy(i,1) > ymaxwp) ymaxwp = wy(i,1)
        do kc = 1,nbas
          if(by(i,ipoint_b(kc)) < yminwp) yminwp = by(i,ipoint_b(kc))
          if(by(i,ipoint_b(kc)) > ymaxwp) ymaxwp = by(i,ipoint_b(kc))
        enddo
      enddo

    endif

    close(iob)
    deallocate(ipoint_b)

  enddo

! projectors (real space)

  if(nump > 0) then

    open(unit = iob, file = 'kbop.gp', status = 'unknown', form = 'formatted')

    write(iob,'("#     r (a.u.)",3x,20(1x,a12))') (labelp(i),i=1,ip)
    write(iob,'("#")')
    do i = 1,numbp(1)
      write(iob,'(20(3x,f10.6))') px(i,1),(py(i,kc),kc=1,ip)
      do kc = 1,ip
        if(py(i,kc) < yminp) yminp = py(i,kc)
        if(py(i,kc) > ymaxp) ymaxp = py(i,kc)
      enddo
    enddo

    close(iob)

  endif

! projectors (reciprocal space)

  if(nump > 0) then

    open(unit = iob, file = 'kbtr.gp', status = 'unknown', form = 'formatted')

    write(iob,'("#     q (1 / a.u.)",3x,20(1x,a12))') (labelq(i),i=1,iq)
    write(iob,'("#")')
    do i = 1,numbq(1)
      write(iob,'(20(3x,f10.6))') qx(i,1),(qy(i,kc),kc=1,iq)
      do kc = 1,iq
        if(qy(i,kc) < yminq) yminq = qy(i,kc)
        if(qy(i,kc) > ymaxq) ymaxq = qy(i,kc)
      enddo
    enddo

    close(iob)

  endif

! integral of kinetic energy (reciprocal space)

  open(unit = iob, file = 'kbek.gp', status = 'unknown', form = 'formatted')

  write(iob,'("#     q (1 / a.u.)",3x,20(1x,a12))') labelek
  write(iob,'("#")')
  do i = 1,numbek
    write(iob,'(f12.6,3x,f18.10)') qek(i),ek(i)
  enddo
  xmaxek = 5*(floor(qek(numbek)/5)+1)
  ymaxek = (10*ONE)**floor(log10(ek(1))+1)
  yminek = (10*ONE)**floor(log10(ek(numbek)))

  close(iob)

! Local potential

  open(unit = iob, file = 'kbvk.gp', status = 'unknown', form = 'formatted')

  write(iob,'("#     r (a.u.)",7x,20(1x,a12))') labelv
  write(iob,'("#")')
  do i = 1,numbv
    write(iob,'(f12.6,3x,f18.10)') vr(i), vk(i)
    if(yminv > vk(i)) yminv = vk(i)
    if(ymaxv < vk(i)) ymaxv = vk(i)
  enddo
  yminv = yminv - 0.2

  close(iob)

  return

end subroutine atom_plot_kb_outdat
