!>  Plots a comparison of the logarithmic derivative ( d psi / dr ) / psi
!>  at the radius rpoint between the all-electron potential,
!>  the semi-local form of the pseudopotential and the non-local
!>  Kleinman-Bylander pseudopotential
!>
!>  \author       Norm Troullier, Manuel Maria Alemany, Jose Luis Martins
!>  \version      6.0.5
!>  \date         1980s and 1990s, 21 August 2021, 25 September 2021.
!>  \copyright    GNU Public License v2

subroutine atom_plot_ln_sub(iowrite, iodata, datafile, iopsd, pseudo,    &
         iokb, pseudokb, ioplot, iotmp, lkb, rpoint, lint)

!  OLD Warnings:

!  This program will construct a d(ln[R])/dr vs Energy plot
!  for a semilocal or K & B pseudopotential.
!  General warnings:  This program is not bullet proof and the
!  user should carefully check their results.  Some of the
!  more probable errors are:
!  1)  The pseudopotential file has fewer potentials then the
!  datafile, ie. lmax of the pseudofile is less then lmax of the
!  datafile.
!  2)  The configuration of the datafile is not the same as the
!  pseudofile.  In this case either create a datafile with the
!  same configuration or use the atom program to modify the
!  charge density of the pseudofile('pm' option).
!  3)  The datafile and pseudfile are not the same method,
!  ie. relativistic, spinpolarized, nonspinpolarized.  Since this
!  routine was designed for nonspin test, a relativistic spin are
!  reduced to nonspin.  A spinpolarized should only be done with
!  a both up and down spin occupations equal(nonspin) for the
!  datafile.  This does not mean a pseudofile with spinpolarized
!  cannot be used, but that a spinpolarized pseduofile shopuld be
!  modified with the atom program ('pm' option) into a
!  nonspinpolarized configuration.

!  The program was written by Norm Troullier at the University of Minnesota
!  It was modified by Manuel Maria Alemany (mmga) at INESC
!  It is maintained by Jose Luis Martins
!  jose.l.martins@tecnico.ulisboa.pt
!  jlmartins@inesc-mn.pt

! printing. 21 October 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  iowrite                          !<  default output

  integer, intent(in)               ::  iodata                           !<  tape number for datafile
  character(len=*), intent(in)      ::  datafile                         !<  file name for all electron results
  integer, intent(in)               ::  iopsd                            !<  tape number for pseudo
  character(len=*), intent(in)      ::  pseudo                           !<  file name for semi-local pseudopotential results
  integer, intent(in)               ::  iokb                             !<  tape number for pseudokb
  character(len=*), intent(in)      ::  pseudokb                         !<  file name for non-local KB pseudopotential results

  integer, intent(in)               ::  ioplot                           !<  tape number for plots
  integer, intent(in)               ::  iotmp                            !<  tape number for temporary command files


  logical, intent(in)               ::  lkb                              !<  If KB results are included in the plots
  real(REAL64), intent(in)          ::  rpoint                           !<  radius where log derivatives are calculated

  logical, intent(in)               ::  lint                             !<  interactive run

! dimensions

  integer                           ::  mxdnr                            !  dimension of radial points
  integer                           ::  mxdorb                           !  dimension of the number of orbitals
  integer                           ::  mxdl                             !  dimension of maximum angular momentum
  integer                           ::  mxdsc                            !  dimension of maximum scattering channels
  integer                           ::  mxdpts                           !  dimension of number of histogram points

! variables

  integer                           ::  nr                               !  dimension of the number of radial points
  real(REAL64), allocatable         ::  r(:)                             !  radial grid points
  real(REAL64), allocatable         ::  drdi(:)                          !  d r(i) / d i
  real(REAL64), allocatable         ::  d2rodr(:)                        !  (d^2 r(i) / d i^2) / (d r / d i)

  integer                           ::  ncore                            !  number of orbitals treated as core

  integer                           ::  norb                             !  number of orbitals

  integer                           ::  lmax                             !  log devs shown up to lmax

  integer, allocatable              ::  no_ae(:)                         !  principal quantum number n
  integer, allocatable              ::  lo_ae(:)                         !  angular quantum number l
  integer, allocatable              ::  iso_ae(:)                        !  2*spin or 2*(j-l)
  real(REAL64), allocatable         ::  zo_ae(:)                         !  orbital occupation
  real(REAL64), allocatable         ::  ev(:)                            !  orbital energy

  integer                           ::  nsc                              !  number of scattering channels
  integer, allocatable              ::  lo_sc(:)                         !  angular quantum number l
  integer, allocatable              ::  iso_sc(:)                        !  2*spin or 2*(j-l)

  integer                           ::  norb_ps                          !  number of orbitals in pseudopotential
  integer, allocatable              ::  lo_ps(:)                         !  angular quantum number l
  integer, allocatable              ::  iso_ps(:)                        !  2*spin or 2*(j-l)

  real(REAL64)                      ::  znuc                             !  nuclear charge

  real(REAL64), allocatable         ::  vionic(:,:,:)                    !  r*ionic potential in Rydberg

  real(REAL64), allocatable         ::  vhxc(:,:)                        !  screening potential in Rydberg

  real(REAL64), allocatable         ::  cdc(:)                           !  4*pi*r**2 charge density of core
  real(REAL64), allocatable         ::  cdv(:,:)                         !  4*pi*r**2 charge density of valence.

  real(REAL64), allocatable         ::  vlocal(:)                        !  local potentia in KB pseudopotentials
  real(REAL64), allocatable         ::  vkbproj(:,:,:)                   !  kb-projector
  integer, allocatable              ::  inorm(:,:)                       !  sign of denominator of KB operator

  real(REAL64), allocatable         ::  ehist(:)                         !  energy in histogram

  real(REAL64), allocatable         ::  dlog_t(:,:)                      !  log derivatives all-electron
  real(REAL64), allocatable         ::  dlog_p(:,:)                      !  log derivatives semi-local pseudopotential
  real(REAL64), allocatable         ::  dlog_k(:,:)                     !  log derivatives non-local KB pseudopotential
  real(REAL64), allocatable         ::  dreg(:)

  real(REAL64), allocatable         ::  esing(:), dsingx(:,:)

  character(len=4), allocatable     ::  c_val(:)

  real(REAL64)                      ::  etot(10)

  character(len=1)                  ::  ispp
  character(len=2)                  ::  icorr, nameat
  character(len=3)                  ::  irel
  character(len=10)                 ::  iray(6), ititle(7)

  integer                           ::  npoint                           !  r(npoint) ~ rpoint

  integer                           ::  itype, ifcore                    !  core correction
  integer                           ::  lmax_v                           !  maximum angular momentum in potential

  real(REAL64)                      ::  zval                             !  valence charge
  real(REAL64)                      ::  zel                              !  electron charge
  real(REAL64)                      ::  evlow, evhi, step                !  energy range and step for plot
  integer                           ::  numtotal                         !  number of points in energy histogram
  integer                           ::  nval

  real(REAL64)                      ::  deltae

  real(REAL64)                      ::  a, b
  character(len=1)                  ::  yesno
  integer                           ::  lmax_in

! parameters

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, ONE = 1.0_REAL64

! counters

  integer     ::  i, j, l, numev


  write(6,*)
  write(6,*) '  Calculating log-derivatives and generating plots'
  write(6,*)


! find dimensions

  call atom_psd_dat_in_size(iodata, datafile, mxdl, mxdnr, mxdorb)
  lmax = mxdl

  if(lint) then
    write(6,*)
    write(6,*) '  Log derivative plots will be shown up to l = ',lmax
    write(6,*) '  Do you want to change that value? (y/n)'
    write(6,*)

    read(5,*) yesno

    if(yesno == 'y' .or. yesno == 'Y') then
      write(6,*)
      write(6,*) '  Enter new maximum angular momentum for plots'
      write(6,*)
      read(5,*) lmax_in
      if(lmax_in < 0 .or. lmax_in > 10) then
        write(6,*)
        write(6,*) '  Wrong value, keeping old value.'
        write(6,*)
      else
        lmax = lmax_in
      endif
    endif

  endif

  if(lint) then
    evlow = -3*ONE
    evhi = 2*ONE
  else
    evlow = -2*ONE
    evhi = 0.5*ONE
  endif
  step = 0.01

  mxdpts = nint((evhi - evlow)/step) +1

! allocate arrays

  allocate(r(mxdnr), drdi(mxdnr), d2rodr(mxdnr))

  allocate(no_ae(mxdorb), lo_ae(mxdorb), iso_ae(mxdorb), zo_ae(mxdorb))
  allocate(ev(mxdorb))

  allocate(lo_ps(mxdorb), iso_ps(mxdorb))

  allocate(vionic(mxdnr,0:mxdl,-1:1))
  allocate(vhxc(mxdnr,-1:1))

  allocate(cdc(mxdnr), cdv(mxdnr,-1:1))

  allocate(vlocal(mxdnr))
  allocate(vkbproj(mxdnr,0:mxdl,-1:1))
  allocate(inorm(0:mxdl,-1:1))

! open and read in data from file DATAFILE

  call atom_psd_dat_in(iodata, datafile, itype, icorr, ispp, lmax_v,     &
      nr, a, b, r, drdi, d2rodr,                                         &
      nameat, norb, ncore, no_ae, lo_ae, iso_ae, zo_ae, znuc, zel,       &
      cdc, vionic, vhxc, ev,                                             &
      mxdnr, mxdorb, mxdl)

  mxdsc = 2*(lmax + 1)
  if(ispp == ' ') then
    mxdsc = lmax + 1
  elseif(ispp == 'r') then
    mxdsc = 2*lmax+1
  endif
  nsc = mxdsc

  allocate(lo_sc(mxdsc),iso_sc(mxdsc))

  if(ispp == ' ') then
    do l = 0,lmax
      lo_sc(l+1) = l
      iso_sc(l+1) = 0
    enddo
  elseif(ispp == 's') then
     do l = 0,lmax
      lo_sc(2*l+1) = l
      iso_sc(2*l+1) =-1
      lo_sc(2*l+2) = l
      iso_sc(2*l+2) = 1
    enddo
  elseif(ispp == 'r') then
    lo_sc(1) = 0
    iso_sc(1) = 1
    if(lmax > 0) then
      do l = 1,lmax
        lo_sc(2*l) = l
        iso_sc(2*l) =-1
        lo_sc(2*l+1) = l
        iso_sc(2*l+1) = 1
      enddo
    endif
  endif

  allocate(ehist(mxdpts))

  allocate(dlog_t(mxdpts,mxdsc))
  allocate(dlog_p(mxdpts,mxdsc), dlog_k(mxdpts,mxdsc))

  allocate(dreg(mxdorb))

  allocate(c_val(mxdsc))


  do j = -1,1
  do i = 1,mxdnr
    cdv(i,j) = ZERO
  enddo
  enddo

  do i = 1,mxdsc
  do numev = 1,mxdpts
    dlog_t(numev,i) = ZERO
  enddo
  enddo

  do i = 1,mxdsc
    if(i < 10) then
      write(c_val(i),'(a3,i1)') 'log', i
    elseif(i < 100) then
      write(c_val(i),'(a2,i2)') 'ln', i
    else

      STOP

    endif
  enddo

  write(iowrite,*)
  write(iowrite,*) '   Plotting log derivatives for :',nameat
  write(iowrite,*)

  zval = ZERO
  do i = ncore+1,norb
    zval = zval + zo_ae(i)
  enddo
!ccccc
! mmga

  write(iowrite,*) '... CALCULATING LOGARITMIC DERIVATIVES ... '
  write(iowrite,*)

! mmga
! mmga You can modify the range of energies changing
! mmga the values of EVLOW and EVHI in this file,
! mmga

!jlm


  numtotal = nint((evhi - evlow)/step) +1
  if(numtotal > mxdpts) numtotal = mxdpts

! mmga
!ccccc

!  Find the point closest to RPOINT

  do i = 1,nr
    if (rpoint < r(i)) then
      npoint = i

      exit

    endif
  enddo

  deltae = (evhi - evlow)/(numtotal -1)
  do numev = 1,numtotal
    ehist(numev) = evlow + (numev-1)*deltae
  enddo

  call atom_plot_ln_nokb(ispp, npoint, numtotal, ehist, dlog_t,          &
      nr, r, drdi, d2rodr,                                               &
      nsc, lo_sc, iso_sc, znuc,                                          &
      lmax_v, vionic, vhxc,                                              &
      iowrite, mxdnr, mxdl, mxdsc, mxdpts)

! The ev(i) and dreg(i) are not passed to the plotting subroutines.
! They are here in case someone wants to add it to the plots.

  allocate(esing(1), dsingx(1,mxdsc))

  do i = ncore+1,norb

    esing(1) = ev(i)
    call atom_plot_ln_nokb(ispp, npoint, 1, esing, dsingx,               &
        nr, r, drdi, d2rodr,                                             &
        1, lo_ae(i:i), iso_ae(i:i), znuc,                                &
        lmax_v, vionic, vhxc,                                            &
        iowrite, mxdnr, mxdl, mxdsc, 1)

    dreg(i) = dsingx(1,1)

  enddo

  write(iowrite,*)
  write(iowrite,*) '  Log derivatives at the eigenvalue energies'
  write(iowrite,*)
  write(iowrite,*)  '  l      s           E         d(ln[R])/dr (E)'

  do i = ncore+1,norb
    write(iowrite,'(i4,f8.2,3x,f12.4,3x,f12.4)') lo_ae(i), 0.5*iso_ae(i:i), ev(i), dreg(i)
  enddo
  write(iowrite,*)

  deallocate(esing, dsingx)


! open and read semi-local data from file pseudo

  call atom_plot_ln_in_psd(iopsd, pseudo, .FALSE.,                       &
      nameat, icorr, irel, ifcore, iray, ititle,                         &
      norb_ps, lo_ps, iso_ps, nr, r, drdi, d2rodr, znuc,                 &
      vionic, cdc, cdv,  vlocal, vkbproj, inorm,                         &
      mxdorb, mxdl, mxdnr)


  nval = norb_ps

  ZVAL = ZNUC
  ncore = 0

  call atom_atm_velect(1, 0, icorr, ' ', ifcore,                         &
      nr, r, drdi, zval, cdv, cdc, vhxc, etot,                           &
      iowrite, mxdnr)

  call atom_plot_ln_nokb(' ', npoint, numtotal, ehist, dlog_p,            &
        nr, r, drdi, d2rodr,                                              &
        nsc, lo_sc, iso_sc, znuc,                                         &
        lmax_v, vionic, vhxc,                                             &
        iowrite, mxdnr, mxdl, mxdsc, mxdpts)


  if(lkb) then

!   open and read non-local KB data from file pseudokb

    call atom_plot_ln_in_psd(iokb, pseudokb, .TRUE.,                     &
        nameat, icorr, irel, ifcore, iray, ititle,                       &
        norb_ps, lo_ps, iso_ps, nr, r, drdi, d2rodr, znuc,               &
        vionic, cdc, cdv,  vlocal, vkbproj, inorm,                       &
        mxdorb, mxdl, mxdnr)

    nval = norb
    itype = ifcore + 1

    ZVAL = ZNUC
    ncore = 0

!   repeated, for consistency

    call atom_atm_velect(1, 0, icorr, ' ', ifcore,                       &
        nr, r, drdi, zval, cdv, cdc, vhxc, etot,                         &
        iowrite, mxdnr)

    call atom_plot_ln_kb(npoint, numtotal, ehist, dlog_k,                &
        nr, r, drdi, d2rodr,                                             &
        nsc, lo_sc, iso_sc,                                              &
        lmax_v, vkbproj, vionic, vhxc, vlocal, inorm,                    &
        iowrite, mxdnr, mxdl, mxdsc, mxdpts)

  endif

  if(lint) then

    call atom_plot_ln_outdis(ioplot, iotmp, lkb,                         &
        numtotal, evlow, evhi, rpoint,                                   &
        c_val, nsc, lo_sc, iso_sc, ehist, dlog_t, dlog_p, dlog_k,        &
        iowrite, mxdsc, mxdpts)

  else

    call atom_plot_ln_outpdf(ioplot, iotmp, lkb,                         &
        numtotal, evlow, evhi, rpoint,                                   &
        c_val, nsc, lo_sc, iso_sc, ehist, dlog_t, dlog_p, dlog_k,        &
        iowrite, mxdsc, mxdpts)

  endif


! deallocate arrays

  deallocate(r, drdi, d2rodr)

  deallocate(no_ae, lo_ae, iso_ae, zo_ae)
  deallocate(lo_sc, iso_sc)
  deallocate(lo_ps, iso_ps)
  deallocate(ev)

  deallocate(vionic)
  deallocate(vhxc)

  deallocate(cdc, cdv)

  deallocate(vlocal)
  deallocate(vkbproj)
  deallocate(inorm)

  deallocate(ehist)

  deallocate(dlog_t, dlog_p, dlog_k)
  deallocate(dreg)

  deallocate(c_val)

  return

end subroutine atom_plot_ln_sub


