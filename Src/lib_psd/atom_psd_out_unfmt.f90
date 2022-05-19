!>  Writes the old-style pseudopotential unformatted file
!>
!>  \author       Sverre Froyen,  Norm Troullier, Jose Luis Martins
!>  \version      6.0.8
!>  \date         1980s and 1990s, 18 May 2022.
!>  \copyright    GNU Public License v2

subroutine atom_psd_out_unfmt(iopsd, filepsd, nameat,                    &
       icorr, irel, nicore, iray, psdtitle,                              &
       npot, nr, a, b, r, zion, indv, ifcore,                            &
       vpsd, cdc, cdpsd,                                                 &
       mxdnr)


! corrected bug reported by Alexander Dobin
! converted to f90, 21 May 2018. JLM
! psdtitle, 18 May 2022. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter                ::  lc = 4                           !  maximum valence angular momentum


! input

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points

  integer, intent(in)               ::  iopsd                            !<  default tape for pseudopotential in old format
  character(len=*), intent(in)      ::  filepsd                          !<  name of default tape for reading pseudopotential in old format

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbot of the atom

  character(len=2), intent(in)      ::  icorr                            !<  correlation type
  character(len=3), intent(in)      ::  irel                             !<  relativistic/spin calculation isp/rel/nrl
  character(len=4), intent(in)      ::  nicore                           !<  codes type of partial core correction  nc/pcec/...

  character(len=10), intent(in)     ::  iray(6)                          !<  code version, date, type of pseudopotential
  character(len=10), intent(in)     ::  psdtitle(20)                     !<  configuration, code details

  integer, intent(in)               ::  npot(-1:1)                       !<  number of  pseudopotentials

  integer, intent(in)               ::  nr                               !<  number of the number of radial points
  real(REAL64), intent(in)          ::  a                                !<  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64), intent(in)          ::  b                                !<  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points

  real(REAL64), intent(in)          ::  zion                             !<  ionic charge
  integer, intent(in)               ::  indv(0:lc,-1:1)                  !<  orbital for scattering channel l

  integer, intent(in)               ::  ifcore                           !<  0 no partial core correction, 1 partial core correction

  real(REAL64), intent(in)          ::  vpsd(mxdnr,0:lc,-1:1)            !<  r*ionic potential in Rydberg (down or total)

  real(REAL64), intent(in)          ::  cdpsd(mxdnr,-1:1)                !<  valence charge density
  real(REAL64), intent(in)          ::  cdc(mxdnr)                       !<  core charge density

! local variables

  integer       ::  nrm
  integer       ::  id                         !  identifies which is the old "down" potential

! counters

  integer     ::  i, j, l

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64


! skip writing file

  if(iopsd == 0) RETURN

! identify the type of calculation

  if(npot(1) /= 0) then
    if(npot(0) /=0) then
      write(6,*)
      write(6,*) '  Inconsistent input, not writing ',filepsd
      write(6,*) '  npot(0), npot(1) = ', npot(0), npot(1)
      write(6,*)

      RETURN

    else
      id = 1
    endif
  else
    if(npot(0) == 0) then
      write(6,*)
      write(6,*) '  Inconsistent input, not writing ',filepsd
      write(6,*) '  npot(0), npot(1) = ', npot(0), npot(1)
      write(6,*)

      RETURN

    else
      id = 0
    endif
  endif

  nrm = nr-1

! old file is iopsd

  open(unit = iopsd, file = filepsd, status = 'unknown', form = 'unformatted')

! keep compatibility with old codes

  write(iopsd) nameat, icorr, irel, nicore, (iray(i),i=1,6),             &
     (psdtitle(i),i=1,7), npot(id), npot(-1), nrm, a, b, zion,           &
     (psdtitle(i),i=8,20)
  write(iopsd) (r(i),i=2,nr)

! Write the potentials to the current pseudo.dat file (unit=1).

  if(irel == 'rel') then

    do l = 0,lc
      if (indv(l,1) /= 0) then
        write(iopsd) l, (((l+1)*vpsd(j,l, 1) + l*vpsd(j,l,-1)) / (2*l+1) ,j=2,nr)
      endif
    enddo

    do l = 1,lc
      if (indv(l,1) /= 0) then
        write(iopsd) l, (2*(vpsd(j,l, 1) - vpsd(j,l,-1)) / (2*l+1) ,j=2,nr)
      endif
    enddo

  else

    do l = 0,lc
      if (indv(l,id) /= 0) then
        write(iopsd) l, (vpsd(j,l,id),j=2,nr)
      endif
    enddo

    do l = 0,lc
      if (indv(l,-1) /= 0) then
        write(iopsd) l,(vpsd(j,l,-1),j=2,nr)
      endif
    enddo

  endif

! Write the charge densities to the current pseudo.dat file (unit=1).

  if (ifcore == 0) then
    write(iopsd) (ZERO,i=2,nr)
  else
    write(iopsd) (cdc(i),i=2,nr)
  endif

!jlm   bug reported by Alexander Dobin

! write(iopsd) (zratio*(cdd(i)+cdu(i)),i=2,nr)
  write(iopsd) ((cdpsd(i, 1)+cdpsd(i,-1)),i=2,nr)

!jlm

  close(unit = iopsd)

  return

end subroutine atom_psd_out_unfmt
