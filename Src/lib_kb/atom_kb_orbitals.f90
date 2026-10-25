!>  Constructs the atomic orbitals used as basis functions.
!>
!>  A basis function with r_bas > 0 is a confined orbital
!>  (SIESTA soft-confinement recipe) with cutoff r_bas.
!>
!>  A basis function with r_bas = 0 is the "relaxed" energy derivative phi
!>  of the previous (confined) basis function psi, in the spirit of LAPW
!>  linearization, but short ranged:
!>  phi solves (H - E) phi = psi - w psi / <psi|w|psi>,  phi(r_c) = 0,
!>  <psi|phi> = 0, where w(r) vanishes for r < r_w.  It is obtained as
!>  phi = (d psi / d alpha) / <psi|w|psi> for the confined eigenfunction
!>  of H + alpha w.  For r < r_w, phi = d psi / d E + c psi.
!>
!>  \author       claude+J.L.Martins
!>  \version      6.1.0
!>  \date         24 October 2026.
!>  \copyright    GNU Public License v2

subroutine atom_kb_orbitals(nr, r, drdi, d2rodr, vionic, vscreen, ev,    &
         npot, lo,                                                       &
         n_bsets, lmax_bas, n_bas, r_bas, nz_bas,                        &
         rpsi_b, lo_b,                                                   &
         iowrite, mxdnr, mxdl, mxdbas, mxdset)

! split from atom_kb_sub. 24 October 2026. claude+JLM
! highest available l potential if l has none. 24 October 2026. JLM+claude
! relaxed energy derivative phi as basis function (r_bas = 0). 24 October 2026. claude+JLM


! The relaxed energy derivative phi follows the spirit of the
! linearization with u and du/dE of LAPW/LMTO:
!   O. K. Andersen, Phys. Rev. B 12, 3060 (1975).
!   D. D. Koelling and G. O. Arbman, J. Phys. F 5, 2041 (1975).
!
! The exact d psi / d E of a confined orbital cannot vanish at r_c:
! integrating <psi|(H-E)|dpsi> = <psi|psi> by parts with psi(r_c) = 0
! gives dpsi(r_c) = +/- 1 / psi'(r_c), large for soft confinement where
! psi'(r_c) is small.  phi relaxes the source term outside r_w, i.e. it
! is the first-order response of psi to the perturbation alpha w(r),
! obtained from an inhomogeneous radial equation as in:
!   R. M. Sternheimer, Phys. Rev. 96, 951 (1954).
!   A. Dalgarno and J. T. Lewis, Proc. R. Soc. London A 233, 70 (1955).
!
! Related ideas for a second basis function per channel:
!   G. Lippert, J. Hutter, P. Ballone, and M. Parrinello,
!   J. Phys. Chem. 100, 6231 (1996), "Response function basis sets",
!   where the second zeta is the derivative of the first with respect
!   to occupation (as quoted by Junquera et al. below).
!   E. Artacho, D. Sanchez-Portal, P. Ordejon, A. Garcia, and J. M. Soler,
!   Phys. Stat. Sol. (b) 215, 809 (1999), split-valence second zeta and
!   perturbative polarization orbitals of SIESTA.
!
! Soft confinement potential:
!   J. Junquera, O. Paz, D. Sanchez-Portal, and E. Artacho,
!   Phys. Rev. B 64, 235111 (2001).

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension of the number of radial points
  integer, intent(in)               ::  mxdl                             !<  dimension for angular momentum
  integer, intent(in)               ::  mxdbas                           !<  dimension for the number of atomic functions
  integer, intent(in)               ::  mxdset                           !<  dimension for number of atomic basis sets

  integer, intent(in)               ::  iowrite                          !<  default output tape

  integer, intent(in)               ::  nr                               !<  number of points in the radial grid
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points   r(i) = a*(exp(b*i)-1), i=1,...,nr
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i
  real(REAL64), intent(in)          ::  d2rodr(mxdnr)                    !<  (d^2 r / d i^2) / (d r / d i)

  real(REAL64), intent(in)          ::  vionic(mxdnr,0:mxdl,-1:1)        !<  r*pseudopotential. -1:  j=l-1/2.  0:  average.  1:  j=l+1/2
  real(REAL64), intent(in)          ::  vscreen(mxdnr)                   !<  screening potential
  real(REAL64), intent(in)          ::  ev(0:mxdl,-1:1)                  !<  eigenvalues (l,2j-2l).
  integer, intent(in)               ::  npot(-1:1)                       !<  number of orbitals (s,p,d,...). -1:  j=l-1/2.  0:  average.  1:  j=l+1/2
  integer, intent(in)               ::  lo(mxdl+1,-1:1)                  !<  angular momentum of orbital

  integer, intent(in)               ::  n_bsets                          !<  number of atomic basis sets
  integer,intent(in)                ::  lmax_bas(mxdset)                 !<  maximum angular momentum in basis
  integer,intent(in)                ::  n_bas(0:mxdl,mxdset)             !<  basis functions for angular momentum l
  real(REAL64), intent(in)          ::  r_bas(4,0:mxdl,mxdset)           !<  cutoff for the basis (up to quadruple zeta), 0: energy derivative of previous
  integer, intent(in)               ::  nz_bas(4,0:mxdl,mxdset)          !<  number of non-trivial zeroes in basis function

! output

  real(REAL64), intent(out)         ::  rpsi_b(mxdnr,mxdbas,mxdset)      !<  basis wavefunctions (r(i),n).
  integer, intent(out)              ::  lo_b(mxdbas,mxdset)              !<  angular momentum of basis function

! local allocatable arrays

  real(REAL64), allocatable         ::  drpsidr_b(:)                     !  drpsi_b / dr
  real(REAL64), allocatable         ::  veff_b(:)                        !  effective potential for basis
  real(REAL64), allocatable         ::  ev_b(:)                          !  energy of the confined orbital or <phi|H|phi>
  integer, allocatable              ::  nrc_b(:)                         !  maximum radius of basis function

  real(REAL64), allocatable         ::  rpsi_p(:), rpsi_m(:)             !  wavefunctions at alpha = +/- h
  real(REAL64), allocatable         ::  dpsi_h(:)                        !  finite difference with step HALPHA
  real(REAL64), allocatable         ::  wrel(:)                          !  weight function w(r) of relaxed derivative
  real(REAL64), allocatable         ::  veff_w(:)                        !  veff_b + alpha w

! local variables

  real(REAL64)                      ::  evl
  integer                           ::  nn
  real(REAL64)                      ::  revi
  real(REAL64)                      ::  vsiesta
  integer                           ::  iflag
  integer                           ::  lpot                             !  angular momentum of the potential used for l

  integer                           ::  nrp, nrm                         !  wave-functions calculated up to r(nrp), r(nrm)
  real(REAL64)                      ::  h                                !  finite difference step
  real(REAL64)                      ::  rw                               !  w(r) = 0 for r < rw
  real(REAL64)                      ::  xw                               !  <psi|w|psi>
  real(REAL64)                      ::  evw                              !  eigenvalue of H + alpha w
  real(REAL64)                      ::  errphi                           !  estimate of relative error in phi
  real(REAL64)                      ::  phinorm                          !  norm of phi
  real(REAL64)                      ::  phide                            !  <phi|H-E|phi> / <phi|phi>
  real(REAL64)                      ::  phimax                           !  max |phi|
  integer                           ::  nodes
  character(len=3)                  ::  kind_b                           !  psi (confined) or phi (energy derivative)
  integer                           ::  ll

! constants

  real(REAL64), parameter           ::  ZERO = 0.0_REAL64
  real(REAL64), parameter           ::  SMALL = 1.0E-6_REAL64
  real(REAL64), parameter           ::  TOL = 1.0E-10_REAL64
  real(REAL64), parameter           ::  EPS = 1.0E-9_REAL64

  integer, parameter                ::  IWTYPE = 2                       !  type of w(r):  1  ((r-rw)/(rc-rw))**NWPOW
                                                                         !                 2  exp(-(rc-rw)/(r-rw))  (smooth onset)
  integer, parameter                ::  NWPOW = 4                        !  power for IWTYPE = 1
  real(REAL64), parameter           ::  RWFAC = 0.5_REAL64               !  rw = RWFAC * r_c
  real(REAL64), parameter           ::  HALPHA = 1.0E-2_REAL64           !  step (Ry) in alpha for d psi / d alpha

  integer, parameter                ::  IOORB = 12                       !  tape for plotting the basis functions
  character(len=*), parameter       ::  FILEORB = 'basis_orbitals.dat'   !  file for plotting the basis functions

! counters

  integer                           ::  i, j, k, n, l, nb


  allocate(drpsidr_b(mxdnr))
  allocate(veff_b(mxdnr))
  allocate(ev_b(mxdbas))
  allocate(nrc_b(mxdbas))

  allocate(rpsi_p(mxdnr),rpsi_m(mxdnr))
  allocate(dpsi_h(mxdnr))
  allocate(wrel(mxdnr))
  allocate(veff_w(mxdnr))

  rpsi_b = ZERO
  lo_b = 0

  open(unit=IOORB, file=FILEORB, form='FORMATTED', status='UNKNOWN')

  write(iowrite,*)
  write(iowrite,'("   Atomic basis functions")')
  write(iowrite,*)
  write(iowrite,'("   psi: confined orbital, E is its eigenvalue (Ry)")')
  write(iowrite,'("   phi: relaxed energy derivative of the previous psi, ",  &
       &   "E = <phi|H|phi> (Ry)")')
  if(IWTYPE == 1) then
    write(iowrite,'("        w(r) = ((r-r_w)/(r_c-r_w))**",i2,                &
       &   ",  r_w = ",f6.3," r_c")') NWPOW, RWFAC
  else
    write(iowrite,'("        w(r) = exp(-(r_c-r_w)/(r-r_w)),  r_w = ",f6.3,   &
       &   " r_c")') RWFAC
  endif
  write(iowrite,*)
  write(iowrite,'("   set  l  zeta  kind  nodes     r_c          E         err")')

  do nb = 1, n_bsets

    n = 0

    do l = 0,lmax_bas(nb)
    do i = 1,n_bas(l,nb)
      n = n+1
      lo_b(n,nb) = l

      errphi = ZERO

      if(r_bas(i,l,nb) > EPS) then

!       confined orbital

        kind_b = 'psi'

        evl =  ev(l,0)
        nn = l+1 + nz_bas(i,l,nb)
        revi = r_bas(i,l,nb)

!       uses the highest available l if there is no potential for l

        lpot = -1
        do k = 1,npot(0)
          if(lo(k,0) == l) lpot = l
        enddo
        if(lpot < 0) then
          do k = 1,npot(0)
            if(lo(k,0) > lpot) lpot = lo(k,0)
          enddo
          if(i == 1) write(iowrite,'("   No potential for l = ",i2,       &
              &   " in basis set ",i3,", using the l = ",i2," potential")') l, nb, lpot
        endif

!       SIESTA recipe (only used up to revi)

        do j = 2,nr
          veff_b(j) = (vionic(j,lpot,0)  + (l*(l+1))/r(j)) / r(j)  + vscreen(j)
          if(r(j) > 0.9*revi .and. r(j) < revi) then
            vsiesta = 40*( exp(-0.1*revi/(r(j)-0.9*revi)) / (revi- r(j) + SMALL) )
            veff_b(j) = veff_b(j) + vsiesta
          endif
!         if(r(j) > 0.8*revi .and. r(j) < revi) then
!           vsiesta = 2.0*(r(j) - 0.8*revi)**3
!           veff_b(j) = veff_b(j) + vsiesta
!         endif
        enddo

        call atom_atm_difnrl_zr(nr, r, drdi, d2rodr,                     &
            veff_b, rpsi_b(:,n,nb), drpsidr_b,                           &
            nn, l, evl, iflag, revi, nrc_b(n), TOL,                      &
            iowrite, mxdnr)

        if(iflag /= 0) then
          write(6,*)
          write(6,*) '  Stopped in atom_kb_orbitals:'
          write(6,*) '  Could not get basis functions, iflag = ',iflag
          write(6,*)

          STOP

        endif

        ev_b(n) = evl

      else

!       relaxed energy derivative of the previous function.
!       veff_b, evl, nn, revi are those of the previous function.

        kind_b = 'phi'

        if(i == 1) then
          write(6,*)
          write(6,*) '  Stopped in atom_kb_orbitals:'
          write(6,*) '  Energy derivative requested for the first ',      &
                     'basis function of l = ',l,' in set ',nb
          write(6,*)

          STOP

        endif

        if(r_bas(i-1,l,nb) <= EPS) then
          write(6,*)
          write(6,*) '  Stopped in atom_kb_orbitals:'
          write(6,*) '  Energy derivative of an energy derivative ',      &
                     'requested for l = ',l,' zeta = ',i,' in set ',nb
          write(6,*)

          STOP

        endif

        rw = RWFAC*revi
        do j = 1,nr
          wrel(j) = ZERO
          if(r(j) > rw .and. r(j) < revi) then
            if(IWTYPE == 1) then
              wrel(j) = ((r(j)-rw)/(revi-rw))**NWPOW
            elseif((revi-rw) < 100*(r(j)-rw)) then
              wrel(j) = exp(-(revi-rw)/(r(j)-rw))
            endif
          endif
        enddo

        xw = ZERO
        ll = 4
        do j = 2,nrc_b(n-1)
          xw = xw + ll*rpsi_b(j,n-1,nb)*wrel(j)*rpsi_b(j,n-1,nb)*drdi(j)
          ll = 6 - ll
        enddo
        xw = xw / 3

!       phi = (d psi / d alpha) / <psi|w|psi>.  Central differences,
!       Richardson extrapolation of steps HALPHA and HALPHA/2

        nrp = 0
        nrm = 0
        do k = 1,2
          h = HALPHA / k

          do j = 1,nr
            veff_w(j) = veff_b(j) + h*wrel(j)
          enddo
          evw = evl
          call atom_atm_difnrl_zr(nr, r, drdi, d2rodr,                   &
              veff_w, rpsi_p, drpsidr_b,                                 &
              nn, l, evw, iflag, revi, nrp, TOL,                         &
              iowrite, mxdnr)
          if(iflag /= 0 .or. nrp /= nrc_b(n-1)) exit

          do j = 1,nr
            veff_w(j) = veff_b(j) - h*wrel(j)
          enddo
          evw = evl
          call atom_atm_difnrl_zr(nr, r, drdi, d2rodr,                   &
              veff_w, rpsi_m, drpsidr_b,                                 &
              nn, l, evw, iflag, revi, nrm, TOL,                         &
              iowrite, mxdnr)
          if(iflag /= 0 .or. nrm /= nrc_b(n-1)) exit

          if(k == 1) then
            do j = 1,nr
              dpsi_h(j) = (rpsi_p(j) - rpsi_m(j)) / (2*h*xw)
            enddo
          else
            do j = 1,nr
              rpsi_b(j,n,nb) = (rpsi_p(j) - rpsi_m(j)) / (2*h*xw)
              errphi = max(errphi, abs(rpsi_b(j,n,nb) - dpsi_h(j)))
              rpsi_b(j,n,nb) = (4*rpsi_b(j,n,nb) - dpsi_h(j)) / 3
            enddo
          endif
        enddo

        if(iflag /= 0 .or. nrp /= nrc_b(n-1) .or. nrm /= nrc_b(n-1)) then
          write(6,*)
          write(6,*) '  Stopped in atom_kb_orbitals:'
          write(6,*) '  Could not get relaxed energy derivative, iflag = ',iflag
          write(6,*) '  nrevi = ', nrc_b(n-1), nrp, nrm
          write(6,*)

          STOP

        endif

        nrc_b(n) = nrc_b(n-1)

!       <phi|H-E|phi> / <phi|phi> = - <phi|w|psi> / (<psi|w|psi> <phi|phi>)

        phinorm = ZERO
        phide = ZERO
        ll = 4
        do j = 2,nrc_b(n)
          phinorm = phinorm + ll*rpsi_b(j,n,nb)*rpsi_b(j,n,nb)*drdi(j)
          phide = phide + ll*rpsi_b(j,n,nb)*wrel(j)*rpsi_b(j,n-1,nb)*drdi(j)
          ll = 6 - ll
        enddo
        phide = -(phide / xw) / phinorm
        phinorm = sqrt(phinorm / 3)

        ev_b(n) = evl + phide

!       normalizes, positive near the origin like psi

        phimax = ZERO
        do j = 2,nrc_b(n)
          phimax = max(phimax, abs(rpsi_b(j,n,nb)))
        enddo
        do j = 2,nrc_b(n)
          if(abs(rpsi_b(j,n,nb)) > phimax/100) exit
        enddo
        if(rpsi_b(j,n,nb) < ZERO) phinorm = -phinorm

        do j = 1,nr
          rpsi_b(j,n,nb) = rpsi_b(j,n,nb) / phinorm
        enddo
        errphi = errphi / phimax

      endif

!     counts nodes ignoring the numerical noise in the tail near r_c

      phimax = ZERO
      do j = 2,nrc_b(n)
        phimax = max(phimax, abs(rpsi_b(j,n,nb)))
      enddo
      nodes = 0
      k = 0
      do j = 2,nrc_b(n)
        if(abs(rpsi_b(j,n,nb)) > phimax/10000) then
          if(k > 0) then
            if(rpsi_b(j,n,nb)*rpsi_b(k,n,nb) < ZERO) nodes = nodes+1
          endif
          k = j
        endif
      enddo

      if(kind_b == 'psi') then
        write(iowrite,'(3x,i3,i3,i5,3x,a3,i6,f11.4,f13.6)')                  &
            nb, l, i, kind_b, nodes, r(nrc_b(n)), ev_b(n)
      else
        write(iowrite,'(3x,i3,i3,i5,3x,a3,i6,f11.4,f13.6,es11.2)')           &
            nb, l, i, kind_b, nodes, r(nrc_b(n)), ev_b(n), errphi
      endif

!     writes the functions for plotting (gnuplot index blocks)

      write(IOORB,'("# set ",i3,"  l ",i3,"  zeta ",i3,"  kind ",a3,       &
          &   "  nodes ",i3,"  r_c ",f10.4,"  E (Ry) ",f14.6)')               &
          nb, l, i, kind_b, nodes, r(nrc_b(n)), ev_b(n)
      write(IOORB,'("#        r            rpsi")')
      do j = 1,nrc_b(n)
        write(IOORB,'(f12.6,es16.7)') r(j), rpsi_b(j,n,nb)
      enddo
      write(IOORB,*)
      write(IOORB,*)

    enddo
    enddo

  enddo

  write(iowrite,*)
  write(iowrite,'("   Basis functions written to ",a)') FILEORB
  write(iowrite,*)

  close(unit=IOORB)

  deallocate(drpsidr_b)
  deallocate(veff_b)
  deallocate(ev_b)
  deallocate(nrc_b)

  deallocate(rpsi_p,rpsi_m)
  deallocate(dpsi_h)
  deallocate(wrel)
  deallocate(veff_w)

  return

end subroutine atom_kb_orbitals
