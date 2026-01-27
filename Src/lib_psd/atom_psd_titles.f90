!>  prepares the titles for the output files
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.1.0
!>  \date         1980s, 29 June 2021, 26 January 2026.
!>  \copyright    GNU Public License v2

subroutine atom_psd_titles(typ, vers, indv, no, zo,                      &
       ispp, ifcore, rc, cfac, zratio, ncore, norb,                      &
       iray, psdtitle, siestatitle, nicore, irel, npot,                  &
       mxdorb)

! extracted from pseud2 code. JLM
! indv. 15 September 2021. JLM
! siestatitle. 26 January 2026. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter                ::  lc = 4                           !  maximum valence angular momentum

! input

  integer, intent(in)               ::  mxdorb                           !<  maximum number of orbitals

  character(len=3), intent(in)      ::  typ                              !<  type of pseudopotential
  character(len=5), intent(in)      ::  vers                             !<  pseudopotential version

  integer, intent(in)               ::  indv(0:lc,-1:1)                  !<  index of orbital

  integer, intent(in)               ::  norb                             !<  number of orbitals
  integer, intent(in)               ::  ncore                            !<  number of orbitals treated as core
  integer, intent(in)               ::  no(mxdorb)                       !<  principal quantum number
  real(REAL64), intent(in)          ::  zo(mxdorb)                       !<  core radii

  character(len=1), intent(in)      ::  ispp                             !<  spin-polarization
  integer, intent(in)               ::  ifcore                           !<  core correction
  real(REAL64), intent(in)          ::  rc(0:lc)                         !<  core radii

  real(REAL64), intent(in)          ::  cfac                             !<  pseudocore stops where rho_core = cfac * rho_val
  real(REAL64), intent(in)          ::  zratio                           !<  ratio between valence charge present and valence charge of neutral atom


! output

  character(len=10), intent(out)    ::  iray(6)                          !<  first title
  character(len=10), intent(out)    ::  psdtitle(20)                     !<  second title
  character(len=70), intent(out)    ::  siestatitle                      !<  second title in new siesta format
  character(len=4), intent(out)     ::  nicore                           !<  type of core correction
  character(len=3), intent(out)     ::  irel                             !<  relativistic (or not)
  integer, intent(out)              ::  npot(-1:1)                       !<  number of potentials

! local variables

  integer            ::  ifull
  real(REAL64)       ::  zelt, zelu, zeld
  integer            ::  noi, ncp
  integer            ::  itext
  character(len=2)   ::  orb
  integer            ::  lmax

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  character(len=1), parameter  ::  IL(0:6) = (/'s','p','d','f','g','h','i'/)

! counters

  integer     ::  i, l


  ncp = ncore+1

  if(typ == 'tm2') then
    iray(1) = 'atom '//vers
    iray(2) = '          '
    call zedate(iray(2))
    iray(3) = '  Improved'
    iray(4) = ' Troullier'
    iray(5) = ' - Martins'
    iray(6) = ' potential'
  else
    do i = 1,6
      iray(i) = '          '
    enddo
  endif


! encode the title array.

  do i = 1,20
    psdtitle(i) = '          '
  enddo

  do i = 1,70
    siestatitle(i:i) = ' '
  enddo

  lmax = 0
  do l = 0,lc
    do i = ncp,norb
      if(ispp == ' ') then
        if(indv(l, 0) == i) then
          if(l > lmax) lmax =l
        endif
      else
        if(indv(l, 1) == i .or. indv(l,-1) == i) then
          if(l > lmax) lmax =l
        endif
      endif
    enddo
  enddo

  if(lmax > 4) lmax = 4                                                  !  in case lc is increased...

! paranoid check for future changes

  if(2*lc + 2 > 20) then
    write(6,*)
    write(6,*) '  Stopped in atom_psd_titles.  Check LC'
    write(6,*)

    STOP

  endif

  if(ispp == ' ') then

    do l = 0,lc

      noi = 0
      orb = '  '

      do i = ncp,norb
        if(indv(l, 0) == i) then
          noi = no(i)
          zelt = zo(i)
          write(orb,'(i1,a1)') noi,IL(l)
        endif
      enddo

      if(noi /= 0) then
        write(psdtitle(2*l+1),'(i1,a1,"(",f6.2,")")') noi,IL(l), zelt
        write(psdtitle(2*l+2),'(a1," rc=",f5.2)') ispp, rc(l)
      endif

!     second line for siesta

      if(l < lmax+1) then
        itext=l*17
        write(siestatitle(itext+1:),'(a2,f5.2,"  r=",f5.2,"/")')             &
             orb, zelt, rc(l)
      endif

    enddo

  elseif(ispp == 'r') then

    do l = 0,lc
      noi = 0
      zelt = ZERO
      orb = '  '

      do i = ncp,norb
        if(indv(l, 1) == i) then
          noi = no(i)
          zelt = zelt + zo(i)
          write(orb,'(i1,a1)') noi,IL(l)
        endif
        if(indv(l,-1) == i) then
          noi = no(i)
          zelt = zelt + zo(i)
          write(orb,'(i1,a1)') noi,IL(l)
        endif
      enddo

      if(noi /= 0) then
        write(psdtitle(2*l+1),'(i1,a1,"(",f6.2,")")') noi,IL(l), zelt
        write(psdtitle(2*l+2),'(a1," rc=",f5.2)') ispp, rc(l)
      endif

!     second line for siesta

      if(l < lmax+1) then
        itext=l*17
        write(siestatitle(itext+1:),'(a2,f5.2,"  r=",f5.2,"/")')             &
             orb, zelt, rc(l)
      endif

   enddo

  elseif(ispp == 's') then

    do l = 0,lc

      noi = 0
      orb = '  '
      zelu = ZERO
      zeld = ZERO

      do i = ncp,norb
        if(indv(l, 1) == i) then
          noi = no(i)
          zeld = zo(i)
          write(orb,'(i1,a1)') noi,IL(l)
        endif
        if(indv(l,-1) == i) then
          noi = no(i)
          zelu = zo(i)
          write(orb,'(i1,a1)') noi,IL(l)
        endif
      enddo

      if(noi /= 0) then
        write(psdtitle(2*l+1),'(i1,a1,"(",f6.2,")")') noi,IL(l), zelu
        write(psdtitle(2*l+2),'(f4.2,")",a1,f4.2)') zeld ,ispp, rc(l)
      endif

!     second line for siesta (is zdown zup reversed?)

      if(l < lmax+1) then
        itext=l*17
        write(siestatitle(itext+1:),'(a2,f4.2,1x,f4.2,1x,f4.2,"/")')     &
             orb, zeld, zelu, rc(l)
      endif

    enddo

  endif

! Determine the number of  potentials.

  do i = -1,1
    npot(i) = 0
    do l = 0,lc
      if (indv(l,i)  /=  0) npot(i) = npot(i) + 1
    enddo
  enddo

! Prepareother heading information

  ifull = 0
  if (cfac <= ZERO .or. zratio == ZERO) ifull = 1
  if (ifcore == 1) then
    if (ifull == 0) then
      nicore = 'pcec'
    else
      nicore = 'fcec'
    endif
  elseif (ifcore == 2) then
    if (ifull == 0) then
      nicore = 'pche'
    else
      nicore = 'fche'
    endif
  else
    nicore = 'nc  '
  endif
  if (ispp == 's') then
    irel='isp'
  elseif (ispp == 'r') then
    irel='rel'
  else
    irel = 'nrl'
  endif

  return

end subroutine atom_psd_titles
