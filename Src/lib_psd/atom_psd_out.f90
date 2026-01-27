!>  Writes an output file that is compatible with the parsec and siesta codes
!>
!>  \author       FDP (parsec), Jose Luis Martins
!>  \version      6.0.8
!>  \date         200s?, May 2018. 18 May 2022. JLM
!>  \copyright    GNU Public License v2

subroutine atom_psd_out(itpars, filename, nameat, cftype,                &
      icorr, irel, nicore, iray, psdtitle, siestatitle,                  &
      npot, nr, a, b, r, zion, indv, ifcore,                             &
      vpsd, cdc, cdpsd,                                                  &
      cfac, rcfac, ncp, norb, lo, rc, zo, rpsi_ps,                       &
      mxdnr, mxdorb, mxdlc, mxdnw)

! it is based on FDP modifications in parsec's fork of the code.

! modified 25 May 2012  difference charge density.
! modified 12 April 2012 zo(norb)
! converted to f90, 21 May 2018. JLM
! Modified rpsi_ps, dimensions, 9 July 2021
! rpsi_ps. 2 November 2021. JLM
! psdtitle, 18 May 2022. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter                ::  lc = 4                           !  maximum valence angular momentum

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points
  integer, intent(in)               ::  mxdorb                           !<  dimension of the number of orbitals
  integer, intent(in)               ::  mxdlc                            !<  maximum angular momentum scattering channel
  integer, intent(in)               ::  mxdnw                            !<  dimension of number of wave-functions same l

  integer, intent(in)               ::  itpars                           !<  tape number
  character(len=*), intent(in)      ::  filename                         !<  file name for writing pseudopotential

  character(len=2), intent(in)      ::  cftype                           !<  indicates output for parsec (PA), siesta (SI) or internal (any other)

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbot of the atom

  character(len=2), intent(in)      ::  icorr                            !<  correlation type
  character(len=3), intent(in)      ::  irel                             !<  relativistic/spin calculation isp/rel/nrl
  character(len=4), intent(in)      ::  nicore                           !<  codes type of partial core correction  nc/pcec/...

  character(len=10), intent(in)     ::  iray(6)                          !<  code version, date, type of pseudopotential
  character(len=10), intent(in)     ::  psdtitle(20)                     !<  configuration, code details
  character(len=70), intent(in)     ::  siestatitle                      !<  second title in new siesta format

  integer, intent(in)               ::  npot(-1:1)                       !<  number of pseudopotentials

  integer, intent(in)               ::  nr                               !<  number of the number of radial points
  real(REAL64), intent(in)          ::  a                                !<  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64), intent(in)          ::  b                                !<  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64), intent(in)          ::  r(nr)                            !<  radial grid points

  real(REAL64), intent(in)          ::  zion                             !<  ionic charge
  integer, intent(in)               ::  indv(0:lc,-1:1)                  !<  orbital for scattering channel l

  integer, intent(in)               ::  ifcore                           !<  0 no partial core correction, 1 partial core correction

  real(REAL64), intent(in)          ::  vpsd(mxdnr,0:lc,-1:1)            !<  r*ionic potential in Rydberg

  real(REAL64), intent(in)          ::  cdpsd(mxdnr,-1:1)                !<  valence charge density
  real(REAL64), intent(in)          ::  cdc(mxdnr)                       !<  core charge density

  real(REAL64), intent(in)          ::  cfac                             !<  core correction factor
  real(REAL64), intent(in)          ::  rcfac                            !<  core correction radius

  integer, intent(in)               ::  ncp                              !<  number of core orbitals + 1
  integer, intent(in)               ::  norb                             !<  number of orbitals

  integer, intent(in)               ::  lo(mxdorb)                       !<  angular quantum number l
  real(REAL64), intent(in)          ::  rc(0:lc)                         !<  orbital core radius
  real(REAL64), intent(in)          ::  zo(mxdorb)                       !<  orbital occupation

  real(REAL64), intent(in)      ::  rpsi_ps(mxdnr,mxdnw,0:mxdlc,-1:1)    !<  r*pseudo-wave-function

! local

  real(REAL64)              ::  zelt

  integer       ::  nrm
  integer       ::  id                         !  identifies which is the old "down" potential

! counters

  integer     ::  i, j, l

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64


! skip writing file

  if(itpars == 0) RETURN

! identify the type of calculation

  if(npot(1) /= 0) then
    if(npot(0) /=0) then
      write(6,*)
      write(6,*) '  Inconsistent input, not writing ',filename
      write(6,*) '  npot(0), npot(1) = ', npot(0), npot(1)
      write(6,*)

      RETURN

    else
      id = 1
    endif
  else
    if(npot(0) == 0) then
      write(6,*)
      write(6,*) '  Inconsistent input, not writing ',filename
      write(6,*) '  npot(0), npot(1) = ', npot(0), npot(1)
      write(6,*)

      RETURN

    else
      id = 0
    endif
  endif

  nrm = nr-1

  open(unit = itpars, file = filename, status = 'unknown', form = 'formatted')

  write(itpars,'(1x,a2,1x,a2,1x,a3,1x,a4)') nameat, icorr, irel, nicore
  write(itpars,'(1x,6a10)') (iray(j),j=1,6)
  if(cftype == 'SI') then
    write(itpars,'(1x,a70)') siestatitle
  elseif(cftype == 'PA') then
    write(itpars,'(1x,7a10)') (psdtitle(j),j=1,7)
  else
    write(itpars,'(1x,20a10)') (psdtitle(j),j=1,20)
  endif

  if(cftype == 'SI') then
    write(itpars,'(1x,2i3,i5,3g20.12)') npot(id), npot(-1), nrm, a, b, zion
  else
    if (ifcore == 0) then
      write(itpars,'(1x,2i3,i5,3g20.12,a20)') npot(id), npot(-1), nrm,     &
           a, b, zion,' nl nls nr a b zion '
    else
      write(itpars,'(1x,2i3,i5,5g20.12,a29)') npot(id), npot(-1), nrm,     &
           a, b, zion, cfac, rcfac,' nl nls nr a b zion cfac rcfac'
    endif
  endif

  write(itpars,'(" Radial grid follows")')
  write(itpars,'(4(g20.12))') (r(j),j=2,nr)

! Write the potentials

  if(irel == 'rel') then

    do l = 0,lc
      if (indv(l,1) /= 0) then
        write(itpars,'(" Average pseudopotential follows (l on next line)")')
        write(itpars,'(1x,i2)') l
        write(itpars,'(4(g20.12))') (((l+1)*vpsd(j,l, 1) + l*vpsd(j,l,-1)) / (2*l+1) ,j=2,nr)
      endif
    enddo

    do l = 1,lc
      if (indv(l,1) /= 0) then
        write(itpars,'(" Spin-orbit pseudopotential follows (l on next line)")')
        write(itpars,'(1x,i2)') l
        write(itpars,'(4(g20.12))') (2*(vpsd(j,l, 1) - vpsd(j,l,-1)) / (2*l+1) ,j=2,nr)
     endif
    enddo

  else

    do l = 0,lc
      if (indv(l,id) /= 0) then
        write(itpars,'(" Pseudopotential follows (l on next line)")')
        write(itpars,'(1x,i2)') l
        write(itpars,'(4(g20.12))') (vpsd(j,l,id),j=2,nr)
      endif
    enddo

    do l = 0,lc
      if (indv(l,-1) /= 0) then
       write(itpars,'(" Minor comp. of Pseudop. follows (l on next line)")')
       write(itpars,'(1x,i2)') l
       write(itpars,'(4(g20.12))') (vpsd(j,l,-1),j=2,nr)
      endif
    enddo

  endif

  write(itpars,'(" Core charge follows")')

  if (ifcore == 0) then
    write(itpars,'(4(g20.12))') (ZERO,i=2,nr)
  else
    write(itpars,'(4(g20.12))') (cdc(i),i=2,nr)
  end if

  write(itpars,'(" Valence charge follows")')
!jlm
       write(itpars,'(4(g20.12))') ((cdpsd(i, 1)+cdpsd(i,-1)),i=2,nr)
!jlm
  do l = 0,lc
!   calculate number of electrons with this angular momentum
    if (indv(l,id) /= 0 .or. indv(l,-1) /= 0) then

      zelt = ZERO
      do j = ncp,norb
        if (lo(j) == l) zelt = zelt + zo(j)
      enddo
      if (indv(l,id) /= 0) then
        write(itpars,'(" Pseudo-wave-function follows (l, zelect, rc)")')
        write(itpars,'(1x,i2,1x,g20.12,1x,g20.12)') l,zelt,rc(l)
        write(itpars,'(4(g20.12))') (rpsi_ps(j,1,l,id),j=2,nr)
      endif

    endif
  enddo

  close(unit = itpars)

  return

end subroutine atom_psd_out
