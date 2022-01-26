!>  Generates the DFT electronic screening potential from
!>  the electron charge density.  The ionic part is
!>  added in dsolv1/dsolv2.
!>
!>  \author       Sverre Froyen, Norm troullier, Carlos Balbas, Jose Luis Martins
!>  \version      6.0.3
!>  \date         1980's, 22 June 2021,12 September 2021.
!>  \copyright    GNU Public License v2

subroutine atom_atm_velect(iter, iconv, icorr, ispp, ifcore,             &
      nr, r, drdi, zel, cdv, cdc, vhxc, etot,                            &
      iowrite, mxdnr)

! This version is a  revision of the exchange-correlation part
! to introduce the interface ATOMXC which allows GGA calculations
! with Perdew'96 and others xc functionals.
! J.M. SOLER and L.C. BALBAS, January 1997.


! converted to f90, February 2018
! cleanup and new interface, July 2019. JLM
! cdv, 12 September 2021. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points

! input

  integer, intent(in)               ::  iter                             !<  iteration number
  integer, intent(in)               ::  iconv                            !<  convergence flag (if iconv = 1, calculates Hartree energy)
  character(len=2), intent(in)      ::  icorr                            !<  correlation type
  character(len=1), intent(in)      ::  ispp                             !<  spin polarization  ' ', 's', 'r'
  integer, intent(in)               ::  ifcore                           !<  0 no partial core correction, 1 partial xc, 2 partial

  integer, intent(in)               ::  nr                               !<  dimension of the number of radial points

  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i

  real(REAL64), intent(in)          ::  cdc(mxdnr)                       !<  4*pi*r**2 * core charge density (total)

! output

  real(REAL64), intent(out)         ::  vhxc(mxdnr,-1:1)                 !<  screening potential in Rydberg

  real(REAL64), intent(out)         ::  etot(10)                         !<  components of total energy


! input and output

  real(REAL64), intent(inout)       ::  zel                              !<  electron charge
  real(REAL64), intent(inout)       ::  r(mxdnr)                         !<  radial grid points

  real(REAL64), intent(inout)       ::  cdv(mxdnr,-1:1)                  !<  4*pi*r**2 * valence charge density

! local allocatable arrays

  real(real64), allocatable         ::  y(:),  yp(:),  ypp(:)            !  work arrays for spline interpolation
  real(real64), allocatable         ::  s1(:), s2(:), w(:)               !  work arrays

  real(REAL64), allocatable         ::  dens(:,:)                        !  charge density
  real(REAL64), allocatable         ::  vxc(:,:)                         !  exchange and correlation potential

! local

  real(REAL64)       ::  pi
  integer            ::  isx
  real(REAL64)       ::  a1, an
  real(REAL64)       ::  b1, bn
  integer            ::  nrm
  integer            ::  ierr
  real(REAL64)       ::  xlo
  integer            ::  ll
  real(REAL64)       ::  xnorm
  real(REAL64)       ::  xccor
  real(REAL64)       ::  ehart
  real(REAL64)       ::  ex, ec, dx, dc
  integer            ::  irel
  integer            ::  nspin

! counters

  integer     ::  i

! constants

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  PFIVE = 0.5_REAL64
  real(REAL64), parameter    ::  OPF = 1.5_REAL64
  real(REAL64), parameter    ::  PNN = 0.99_REAL64


  allocate(y(mxdnr),yp(mxdnr),ypp(mxdnr))
  allocate(s1(mxdnr),s2(mxdnr),w(6*mxdnr))

  allocate(dens(mxdnr,2))
  allocate(vxc(mxdnr,2))

  do i = 1,10
    etot(i) = ZERO
  enddo

  pi = 4*atan(ONE)

! fit cd/r by splines

  y(1) = ZERO
  do i = 2,nr
    y(i) = cdv(i,0) / r(i)
  enddo

  if (ifcore == 2) then
    do i = 2,nr
      y(i) = y(i) + cdc(i)/r(i)
    enddo
  endif

  isx = 0
  a1 = ZERO
  an = ZERO
  b1 = ZERO
  bn = ZERO
  nrm=nr

  call splift(r,y,yp,ypp,nrm,w,ierr,isx,a1,b1,an,bn)

  if(ierr /= 1) then
    write(iowrite,'(1x,"****** Error in splift ierr =",i2)') ierr

    stop

  endif

! compute the integrals of cd/r and cd from
! r(1)=0 to r(i)

  xlo = ZERO

  call spliq(r,y,yp,ypp,nrm,xlo,r,nrm,s2,ierr)

  if(ierr /= 1) then
    write(iowrite,'(1x,"****** Error in spliq ierr =",i2)') ierr

    stop

  endif

  do i=1,nr
    ypp(i) = r(i)*ypp(i) + 2*yp(i)
    yp(i)  = r(i)*yp(i)  + y(i)
    y(i)   = r(i)*y(i)
  enddo

  call spliq(r,y,yp,ypp,nrm,xlo,r,nrm,s1,ierr)

  if(ierr /= 1) then
    write(iowrite,'(1x,"****** Error in spliq ierr =",i2)') ierr

    stop

  endif

! check normalization

  xnorm = ZERO
  if (ifcore == 2 .and. iter == 0 ) zel = s1(nr)
  if (zel /= ZERO) xnorm = zel/s1(nr)
  if (iter  >  3 .and. abs(zel-s1(nr))  >  0.01) then
    if (zel < s1(nr) + ONE) then
      write(iowrite,'(/," warning *** charge density rescaled in velect",  &
            & /," iteration number",i4,3x,"scaling factor =",f6.3,/)')   &
             iter,xnorm
    else
      xnorm = pnn*xnorm
      write(iowrite,'(/," warning *** charge density partially ",        &
            & "rescaled in velect",/," iteration number",i4,3x,          &
            & "scaling factor =",f6.3,/)') iter,xnorm
    endif
  endif

! compute new hartree potential
! renormalize the charge density

  do i = 2,nr
    vhxc(i,-1) = 2 * xnorm*(s1(i)/r(i) + s2(nr) - s2(i))
    vhxc(i, 0) = vhxc(i,-1)
    vhxc(i, 1) = vhxc(i,-1)
    cdv(i,-1) = xnorm*cdv(i,-1)
    cdv(i, 0) = xnorm*cdv(i, 0)
    cdv(i, 1) = xnorm*cdv(i, 1)
  enddo

! compute hartree contribution to total energy

  ehart = ZERO

  if (iconv == 1) then
    ehart = ZERO
    ll = 4
    do i = 2,nr
      ehart = ehart + ll*cdv(i, 0)*vhxc(i,0)*drdi(i)
      ll = 6 - ll
    enddo
    ehart = ehart / 6
  endif

!Add exchange and correlation
!This part is totally new. It is written to use
!the package XC of J.M. Soler.
!J.M. Soler and L.C. Balbas. January 1997

!Compute dens(i,nspin) = density up, density down

  if (nr > mxdnr) then
     write(iowrite,'(" velect: mxdnr too small. nr, mxdnr = ",2i8)')     &
             nr, mxdnr

    stop

  endif

  do i = 2,nr
    if (ispp == 's') then
      dens(i,1) = cdv(i,-1) / (4*pi*r(i)**2)
      dens(i,2) = cdv(i, 1) / (4*pi*r(i)**2)
    else
      dens(i,1) = PFIVE*cdv(i, 0) / (4*pi*r(i)**2)
      dens(i,2) = dens(i,1)
    endif
    if (ifcore > 0) then
          dens(i,1) = dens(i,1) + PFIVE*cdc(i) / (4*pi*r(i)**2)
          dens(i,2) = dens(i,2) + PFIVE*cdc(i) / (4*pi*r(i)**2)
    endif
  enddo

! Extrapolate the density at r=0

  dens(1,1) = dens(2,1) - (dens(3,1)-dens(2,1))*r(2)/(r(3)-r(2))
  dens(1,2) = dens(2,2) - (dens(3,2)-dens(2,2))*r(2)/(r(3)-r(2))

! Define 'irel' and 'nspin' for the interface ATOMXC

  if (ispp == 'r') then
    irel = 1
  else
    irel = 0
  endif
  nspin = 2

  r(1) = ZERO

  if (icorr == 'ca') then
    call atom_xc('LDA', 'ca', irel, nr, mxdnr, r, nspin, dens,           &
                ex, ec, dx, dc, vxc, iowrite)
  elseif(icorr == 'pw') then
    call atom_xc('LDA', 'pw92', irel, nr, mxdnr, r, nspin, dens,         &
                ex, ec, dx, dc, vxc, iowrite)
  elseif(icorr == 'pb') then
    call atom_xc('GGA', 'pbe', irel, nr, mxdnr, r, nspin, dens,          &
                 ex, ec, dx, dc, vxc, iowrite)
  endif

! Add vxc to total potential and energies

  do i = 2,nr
    vhxc(i,-1) = vhxc(i,-1) + vxc(i,1)
    vhxc(i, 0) = vhxc(i, 0) + (vxc(i,1) + vxc(i,2)) / 2
    vhxc(i, 1) = vhxc(i, 1) + vxc(i,2)
  enddo

!lcb
  xccor = ZERO
  ll = 4
  do i = 2,nr
    xccor = xccor + ll * drdi(i) * (vxc(i,1)*cdv(i,-1) + vxc(i,2)*cdv(i, 1))
    ll = 6 - ll
  enddo

  etot(4) = ehart
  etot(5) = xccor / 3
  etot(6) = xccor - 4*(ex + ec)
  etot(7) = ex + ec
!lcb

! Obtain total potential at r(1) = 0

  do i = -1,1
    vhxc(1,i) = vhxc(2,i) - (vhxc(3,i)-vhxc(2,i))*r(2)/(r(3)-r(2))
  enddo

! *** lcb-jms modification end ********************

  deallocate(y,yp,ypp)
  deallocate(s1,s2,w)

  deallocate(dens)
  deallocate(vxc)

  return

end subroutine atom_atm_velect

