!>  Finds largest zero and maximum of all electron radial wave-functions
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.7
!>  \date         28 August 2021, 26 January 2022.
!>  \copyright    GNU Public License v2

subroutine atom_psd_max_zero(rz, rx, mxdlc, iowrite, ioae, fileae)

! new structure indx. 2 November 2021. JLM
! nzero = 0. 26 january 2022. JLM


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

  integer, parameter                ::  mxdnw = 2                        !  hard coded dimension of number of wave-functions same l

! input

  integer, intent(in)               ::  mxdlc                            !<  maximum angular momentum scattering channel

  integer, intent(in)               ::  iowrite                          !<  default output tape

  integer, intent(in)               ::  ioae                             !<  default tape for writing
  character(len=*), intent(in)      ::  fileae                           !<  name of default tape for all-electron results

! output

  real(REAL64), intent(out)         ::  rz(0:mxdlc)                      !<  average zero for l
  real(REAL64), intent(out)         ::  rx(0:mxdlc)                      !<  average extrema for l

! dimensions

  integer                           ::  mxdnr                            !  dimension of the number of radial points
  integer                           ::  mxdorb                           !  dimension of the number of orbitals
  integer                           ::  mxdl                             !  dimension of maximum angular momentum


! variables

  integer                           ::  nr                               !  number of radial points

  integer                           ::  norb                             !  number of orbitals
  integer                           ::  lmax                             !  maximum angular momentum+1

  integer                           ::  itype                            !  unused
  character(len=2)                  ::  icorr                            !  correlation type
  character(len=1)                  ::  ispp                             !  spin polarization  ' ', 's', 'r'

  real(REAL64)                      ::  a                                !  r(i) = a*(exp(b*(i-1))-1)
  real(REAL64)                      ::  b                                !  r(i) = a*(exp(b*(i-1))-1)

  character(len=2)                  ::  nameat                           !  chemical symbot of the atom

  integer                           ::  ncore                            !  number of orbitals treated as core

  real(REAL64)                      ::  znuc                             !  nuclear charge
  real(REAL64)                      ::  zel                              !  electron charge



  real(REAL64)                      ::  rzero
  real(REAL64)                      ::  rextr
  integer                           ::  nextr, nzero
  integer                           ::  nn

  integer                           ::  lcmax

! allocatable arrays

  real(REAL64), allocatable         ::  r(:)                             !  radial grid points
  real(REAL64), allocatable         ::  drdi(:)                          !  d r(i) / d i
  real(REAL64), allocatable         ::  d2rodr(:)                        !  (d^2 r(i) / d i^2) /  (d r(i) / d i)

  integer, allocatable              ::  no(:)                            !  principal quantum number n
  integer, allocatable              ::  lo(:)                            !  angular quantum number l
  integer, allocatable              ::  iso(:)                           !  2*spin or 2*(j-l)
  real(REAL64), allocatable         ::  zo(:)                            !  orbital occupation

  real(REAL64), allocatable         ::  cdc(:)                           !  core charge density (total)

  real(REAL64), allocatable         ::  vionic(:,:,:)                    !  r*ionic potential in Rydberg

  real(REAL64), allocatable         ::  vhxc(:,:)                        !  effective potential in Rydberg

  real(REAL64), allocatable         ::  ev(:)                            !  orbital energy

  real(REAL64), allocatable         ::  v(:)                             !  effective potential

  real(REAL64), allocatable         ::  ar(:)                            !  r*psi or major component
  real(REAL64), allocatable         ::  br(:)                            !  d ar / d r or minor component

  real(REAL64), allocatable         ::  rz_all(:)
  real(REAL64), allocatable         ::  rx_all(:)
  real(REAL64), allocatable         ::  ax_all(:)
  real(REAL64), allocatable         ::  bx_all(:)

  integer, allocatable              ::  indv(:,:)                        !  main orbital for scattering channel l (legacy)
  integer, allocatable              ::  indx(:,:,:)                      !  orbitals for scattering channel l
  integer, allocatable              ::  nindx(:,:)                       !  number of orbitals for scattering channel l

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64

! counter

  integer     ::  is, i, j, l, n, nxz


! get array sizes

  call atom_psd_dat_in_size(ioae, fileae, mxdl, mxdnr, mxdorb)

! allocates arrays

  allocate(r(mxdnr))
  allocate(drdi(mxdnr))
  allocate(d2rodr(mxdnr))

  allocate(no(mxdorb))
  allocate(lo(mxdorb))
  allocate(iso(mxdorb))
  allocate(zo(mxdorb))

  allocate(cdc(mxdnr))
  allocate(vionic(mxdnr,0:mxdl,-1:1))
  allocate(vhxc(mxdnr,-1:1))

  allocate(ev(mxdorb))

  allocate(v(mxdnr))

  allocate(ar(mxdnr))
  allocate(br(mxdnr))

! reads atomic data

  call atom_psd_dat_in(ioae, fileae, itype, icorr, ispp, lmax,           &
      nr, a, b, r, drdi, d2rodr,                                         &
      nameat, norb, ncore, no, lo, iso, zo, znuc, zel,                   &
      cdc, vionic, vhxc, ev,                                             &
      mxdnr, mxdorb, mxdl)

! fills indv, indx

  allocate(indv(0:mxdlc,-1:1))
  allocate(indx(mxdnw,0:mxdlc,-1:1))
  allocate(nindx(0:mxdlc,-1:1))

  call atom_psd_ind(indv, indx, nindx, lcmax, norb, ncore, lo, iso,      &
        mxdorb, mxdlc, mxdnw)

  nn = 3
  do i = ncore+1,norb
    nn = max(no(i)-lo(i) + 2,nn)
  enddo

  allocate(rz_all(nn))
  allocate(rx_all(nn))
  allocate(ax_all(nn))
  allocate(bx_all(nn))

! begin loop over valence orbitals

  do l = 0,lcmax
    rz(l) = ZERO
    rx(l) = ZERO
  enddo

  do l = 0,lcmax
    nxz = 0
    do is = -1,1

      if(nindx(l,is) > 0) then

!       only lowest n orbital is important

        n = 1
        do j = 1,nindx(l,is)
          if(no(indx(j,l,is)) < no(indx(n,l,is))) n = j
        enddo

        i = indx(n,l,is)
        call atom_psd_ae_orb(ispp, nr, r, drdi, d2rodr, v, ar, br,       &
            no(i), lo(i), iso(i), znuc, ev(i), .FALSE.,                  &
            vionic(:,lo(i),iso(i)), vhxc(:,iso(i)),                      &
            iowrite, mxdnr)

        call atom_atm_zero_extr(ispp, ar, br, nr, r,                     &
            no(i), lo(i), iso(i) ,vionic, vhxc, ev,                      &
            nzero, nextr, rz_all, rx_all, ax_all, bx_all,                &
            mxdnr)

        if(nzero == 0) then
          rzero = ZERO
        else
          rzero = rz_all(nzero)
        endif
        rextr = rx_all(nextr)

        rz(l) = rz(l) + rzero
        rx(l) = rx(l) + rextr
        nxz = nxz + 1

      endif

    enddo

    if(nxz /= 0) then
      rz(l) = rz(l) / nxz
      rx(l) = rx(l) / nxz
    endif

 enddo

! deallocates arrays

  deallocate(r)
  deallocate(drdi)
  deallocate(d2rodr)

  deallocate(no)
  deallocate(lo)
  deallocate(iso)
  deallocate(zo)

  deallocate(cdc)
  deallocate(vionic)
  deallocate(vhxc)

  deallocate(ev)

  deallocate(v)

  deallocate(ar)
  deallocate(br)

  deallocate(indv)
  deallocate(indx)
  deallocate(nindx)

  deallocate(rz_all)
  deallocate(rx_all)
  deallocate(ax_all)
  deallocate(bx_all)

  return

end subroutine atom_psd_max_zero
