!>  tests the eigenvalues and recodes vpsd in average and difference
!>
!>  \author       Norm Troullier, Jose Luis Martins
!>  \version      6.0.5
!>  \date         1980s and 1990s, 30 June 2021, 15 Septenber 2021.
!>  \copyright    GNU Public License v2


subroutine atom_psd_test_simple(ifcore, icorr, ispp, rc,                 &
    nr, r, drdi, d2rodr,                                                 &
    nameat, norb, ncore, lo, iso, zo, znuc, zel,                         &
    cdpsd, cdc, vpsd,                                                    &
    ev, evi, indv, zratio, zion,                                         &
    iowrite, mxdnr, mxdorb)

! Writen 13 April 2018 from the old codes.  The old
! pseudopotential generation subroutines were broken in several subroutines
! and converted to f90.
! Modified (cleaning) 8 July 2021. JLM
! so->iso, vpsd. 15 Septenber 2021. JLM
! no virial en etotal. 18 October 2021. JLM

!mmga  modifications from early Sverre code by Manuel Maria Gonzalez Alemany
!njtj  modifications from early Sverre code by Norm Troullier

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter                ::  lc = 4                           !  maximum valence angular momentum

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points
  integer, intent(in)               ::  mxdorb                           !<  dimension of the number of orbitals

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  ifcore                           !<  0: no partial core correction; 1 : xc partial core correction; 2 : not recommended.
  character(len=2), intent(in)      ::  icorr                            !<  correlation type
  character(len=1), intent(in)      ::  ispp                             !<  spin polarization  ' ', 's', 'r'

  real(REAL64), intent(in)          ::  rc(0:lc)                         !<  core radius r_c(l)

  integer, intent(in)               ::  nr                               !<  number of radial grid points
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i
  real(REAL64), intent(in)          ::  d2rodr(mxdnr)                    !<  (d^2 r(i) / d i^2) /  (d r(i) / d i)


  character(len=2), intent(in)      ::  nameat                           !<  chemical symbot of the atom

  integer, intent(in)               ::  norb                             !<  number of orbitals
  integer, intent(in)               ::  ncore                            !<  number of orbitals treated as core
  integer, intent(in)               ::  lo(mxdorb)                       !<  angular quantum number l
  integer, intent(in)               ::  iso(mxdorb)                      !<  2*spin or 2*(j-l)
  real(REAL64), intent(in)          ::  zo(mxdorb)                       !<  orbital occupation

  real(REAL64), intent(in)          ::  znuc                             !<  nuclear charge
  real(REAL64), intent(in)          ::  zel                              !<  electron charge

  real(REAL64), intent(in)          ::  evi(norb)                        !<  fixed orbital energy

! output

  real(REAL64), intent(out)         ::  zratio
  real(REAL64), intent(out)         ::  zion

! input and output

  real(REAL64), intent(inout)       ::  cdpsd(mxdnr,-1:1)                !<  valence charge density (down or total)
  real(REAL64), intent(inout)       ::  cdc(nr)                          !<  core charge density

  real(REAL64), intent(inout)       ::  vpsd(mxdnr,0:lc,-1:1)            !<  pseudo-potential in Rydberg (down or total or average in output if relativistic)

  real(REAL64), intent(inout)       ::  ev(norb)                         !<  orbital energy

  integer, intent(in)               ::  indv(0:lc,-1:1)                  !<  orbital for scattering channel l


! local allocatable arrays

  real(real64), allocatable         ::  vhxc(:,:)
  real(REAL64), allocatable         ::  vout(:,:)                        !<  r*effective potential in Rydberg (up or not used)

  integer, allocatable              ::  nops(:)
  real(REAL64), allocatable         ::  ek(:), ep(:)

! local variables

  integer            ::  ncp

  real(REAL64)       ::  etot(10)                                        !  components of total energy

  real(REAL64)       ::  zval, zval2
  character(len=1)   ::  blank

  real(REAL64)       ::  zot

  real(REAL64)       ::  vpsdm, rmind
  real(REAL64)       ::  vps
  
  integer            ::  nlc(0:lc,-1:1)

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64
  character(len=1), parameter  ::  IL(0:6) = (/'s','p','d','f','g','h','i'/)

! counters

  integer     ::  i, j, l



! Reset the n quantum numbers to include all valence orbitals.
! Compute the ratio between the valence charge present and the
! valence charge of a neutral atom.
! Transfer pseudo valence charge to charge array

  allocate(nops(norb))
  allocate(ek(norb),ep(norb))

  ncp = ncore + 1
  zval = ZERO
  zratio = ZERO
  do i = 1,norb
    nops(i) = 0
  enddo
  do j = -1,1
  do l = 0,lc
    nlc(l,j) = 0
  enddo
  enddo
  do i = ncp,norb
    nlc(lo(i),iso(i)) = nlc(lo(i),iso(i)) + 1
    nops(i) = lo(i) + nlc(lo(i),iso(i))
    zval = zval + zo(i)
  enddo
  zion = zval + znuc - zel
  if (zval  /=  ZERO) zratio = zion/zval


! Convert spin-polarized potentials back to nonspin-polarized
! by occupation weight(zo).  Assumes core polarization is
! zero, ie. polarization is only a valence effect.

! ASSUMES A WELL ORDERED ORBITALS WITH ALTERNATING SPINS  WARNING

  if (ispp == 's' ) then
    do i = ncp,norb,2
      if(lo(i) /= lo(i+1)) then
        write(6,*) '  stopped in psd_test'
        write(6,*) '  non-alternating spin-orpitals'

        STOP

      endif
      zot = zo(i) + zo(i+1)
      if (zot /=  ZERO) then
        do j = 2,nr
          vpsd(j,lo(i), 0) = (vpsd(j,lo(i),iso(i))*zo(i) +               &
                           vpsd(j,lo(i+1),iso(lo(i)+1))*zo(i+1)) / zot
        enddo
      else
        do j = 2,nr
          vpsd(j,lo(i), 0) = (vpsd(j,lo(i),iso(i)) +               &
                           vpsd(j,lo(i+1),iso(lo(i)+1))) / 2
        enddo
      endif
    enddo
  endif

  allocate(vout(mxdnr,-1:1))

  zval2 = zval
  blank = ' '
  call atom_atm_velect(0, 1, icorr, blank, ifcore,                       &
      nr, r, drdi, zval, cdpsd, cdc, vout, etot,                         &
      iowrite, mxdnr)

  allocate(vhxc(mxdnr,mxdorb))

  do j = 1,norb
    do i = 1,nr
      vhxc(i,j) = vout(i,iso(j))
    enddo
  enddo


! Test the pseudopotential self consistency.  Spin-polarized
! is tested as non-spin-polarized(since up/down potentials are
! now the same)

  call atom_atm_dsolv2(0, 1 ,blank, ifcore,                              &
      nr, r, drdi, d2rodr,                                               &
      norb-ncore, 0, nops(ncp), lo(ncp), iso(ncp), zo(ncp), znuc,        &
      cdpsd, cdc, vpsd, vhxc(:,ncp),                                     &
      ev(ncp), ek(ncp), ep(ncp), evi(ncp),                               &
      iowrite, mxdnr, mxdorb-ncore, lc)

! Printout the pseudo eigenvalues after cutoff.

  write(iowrite,*)
  write(iowrite,*)
  write(iowrite,'(" test of eigenvalues")')
  write(iowrite,*)
  write(iowrite,'(" rcut =",8(2x,a1,f7.2))') (IL(lo(i)),rc(lo(i)),i=ncp,norb)
  write(iowrite,'(" eval =",8(2x,f8.5))') (ev(i),i=ncp,norb)

! Printout the data for potentials.

  write(iowrite,*)
  write(iowrite,*)
  write(iowrite,*)
  write(iowrite,'(" l    vps(0)    vpsmin      at r")')
  write(iowrite,*)

  do l = 0,lc
    do i = ncp,norb
      if(ispp == ' ') then

        if(indv(l, 0) == i) then
          vpsdm = ZERO
          do j = 2,nr
            if (r(j) >= .00001) then
              vps = vpsd(j,l, 0)/r(j)
              if (vps < vpsdm) then
                vpsdm = vps
                rmind = r(j)
              endif
            endif
          enddo
          write(iowrite,'(1x,a1,3f10.3)') IL(l),vpsd(2,l, 0)/r(2),vpsdm,rmind
        endif

      else

        if(indv(l, 1) == i) then
          vpsdm = ZERO
          do j = 2,nr
            if (r(j) >= .00001) then
              vps = vpsd(j,l, 1)/r(j)
              if (vps < vpsdm) then
                vpsdm = vps
                rmind = r(j)
              endif
            endif
          enddo
          write(iowrite,'(1x,a1,3f10.3)') IL(l),vpsd(2,l, 1)/r(2),vpsdm,rmind
        endif

        if(indv(l,-1) == i) then
          vpsdm = ZERO
          do j = 2,nr
            if (r(j) >= .00001) then
              vps = vpsd(j,l,-1)/r(j)
              if (vps < vpsdm) then
                vpsdm = vps
                rmind = r(j)
              endif
            endif
          enddo
          write(iowrite,'(1x,a1,3f10.3)') IL(l),vpsd(2,l,-1)/r(2),vpsdm,rmind
        endif

      endif
    enddo
  enddo

! Print out the energies from etotal.

  call atom_atm_etotal(4, nameat, norb-ncore,                            &
      nops(ncp), lo(ncp), iso(ncp), zo(ncp),                             &
      etot, ev(ncp), ek(ncp), ep(ncp),                                   &
      iowrite, norb)

  write(iowrite,*)
  write(iowrite,*)

  deallocate(nops)
  deallocate(vout)
  deallocate(vhxc)

  deallocate(ek,ep)

  return

end subroutine atom_psd_test_simple
