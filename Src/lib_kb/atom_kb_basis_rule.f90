!>  Finds the parameters that describe the default atomic basis sets
!>  based on the SIESTA radii.  One set of rules for all elements.
!>
!>
!>  \author       claude+J.L.Martins
!>  \version      6.1.0
!>  \date         3 November 2026.
!>  \copyright    GNU Public License v2

subroutine atom_kb_basis_rule(n_bsets, tbasis, lmax_pot, l_occ,          &
         r_siesta, lmax_bas, n_bas, r_bas, nz_bas,                       &
         iowrite, mxdl, mxdset)

! written 3 November 2026 from the rules in atom_kb_basis_cutoff. claude+JLM


! Angular momenta:  atomic orbitals for l = 0,...,la with
! la = max(highest occupied l, lmax_pot - 1), and, if the type of basis
! ends with P, polarization orbitals for l = la + 1 (if la + 1 <= LPOLMAX).
! (H: s + p polarization.  Li, Be, Na: s + p pol.  Mg, K, Ca,...: s, p + d pol.
! p-block: s, p + d pol.  Transition metals: s, p, d + f pol.
! Lanthanides and actinides with f electrons: s, p, d, f.)
!
! Cutoff radii:
!   atomic orbitals       r_c(l) = min( max( F_l r_siesta(l), RMIN ), RMAX )
!                         F_l = FAC_S for l = 0, FAC_L for l > 0
!   polarization          r_c = min( max( FAC_POL max(r_siesta(0),r_siesta(1)), RMIN ), RMAX )
!
! SZ: one function for each l.  DZ: the second function is the relaxed
! energy derivative of the first (r_bas = 0).
!
! The constants were chosen comparing the atomic orbital (DZP) and
! plane-wave total energies with cpw for the elemental solids of
! H to U (with the confinement of atom_kb_orbitals), November 2026.

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdl                             !<  array dimension for basis angular momentum
  integer, intent(in)               ::  mxdset                           !<  dimension for number of atomic basis sets

  integer, intent(in)               ::  iowrite                          !<  default output tape

  integer, intent(in)               ::  n_bsets                          !<  number of atomic basis sets
  character(len=3), intent(in)      ::  tbasis(mxdset)                   !<  type of basis (SZ, DZ, SZP, DZP)

  integer, intent(in)               ::  lmax_pot                         !<  maximum angular momentum for potential
  integer, intent(in)               ::  l_occ                            !<  highest angular momentum of a bound occupied orbital
  real(REAL64), intent(in)          ::  r_siesta(0:mxdl)                 !<  cutoff using the SIESTA recipe

! output

  integer, intent(out)              ::  lmax_bas(mxdset)                 !<  maximum angular momentum in basis
  integer, intent(out)              ::  n_bas(0:mxdl,mxdset)             !<  basis functions for angular momentum l
  real(REAL64), intent(out)         ::  r_bas(4,0:mxdl,mxdset)           !<  cutoff for the basis (up to quadruple zeta), 0: energy derivative of previous
  integer, intent(out)              ::  nz_bas(4,0:mxdl,mxdset)          !<  number of non-trivial zeroes in basis function

! local variables

  integer                           ::  la                               !  highest l of the atomic (not polarization) orbitals
  integer                           ::  nzeta                            !  1 (SZ) or 2 (DZ)
  real(REAL64)                      ::  rc_l(0:mxdl)                     !  cutoff of the atomic orbitals
  real(REAL64)                      ::  rc_pol                           !  cutoff of the polarization orbitals
  character(len=3)                  ::  ch_temp

! constants

  real(REAL64), parameter           ::  ZERO = 0.0_REAL64

  real(REAL64), parameter           ::  FAC_S = 1.45_REAL64              !  r_c / r_siesta for l = 0
  real(REAL64), parameter           ::  FAC_L = 1.15_REAL64              !  r_c / r_siesta for l > 0
  real(REAL64), parameter           ::  FAC_POL = 1.0_REAL64             !  r_c / max(r_siesta(0),r_siesta(1)) for polarization
  real(REAL64), parameter           ::  RMIN = 6.0_REAL64                !  smallest cutoff (bohr)
  real(REAL64), parameter           ::  RMAX = 9.0_REAL64                !  largest cutoff (bohr)
  integer, parameter                ::  LPOLMAX = 3                      !  no polarization orbitals above f

! counters

  integer                           ::  i, l


  la = max(l_occ, lmax_pot-1, 0)
  if(la > mxdl) then
    write(6,*)
    write(6,*) '  Stopped in atom_kb_basis_rule:'
    write(6,*) '  la = ',la,' mxdl = ',mxdl

    STOP

  endif

  do l = 0,la
    if(l == 0) then
      rc_l(l) = min( max( FAC_S*r_siesta(l), RMIN ), RMAX )
    else
      rc_l(l) = min( max( FAC_L*r_siesta(l), RMIN ), RMAX )
    endif
  enddo

  if(lmax_pot > 0) then
    rc_pol = min( max( FAC_POL*max(r_siesta(0),r_siesta(1)), RMIN ), RMAX )
  else
    rc_pol = min( max( FAC_POL*r_siesta(0), RMIN ), RMAX )
  endif

  do i = 1,n_bsets

    ch_temp = tbasis(i)
    call chrcap(ch_temp,3)

    nzeta = 1
    if(ch_temp(1:2) == 'DZ') nzeta = 2

    n_bas(:,i) = 0
    r_bas(:,:,i) = ZERO
    nz_bas(:,:,i) = 0

    lmax_bas(i) = la
    do l = 0,la
      n_bas(l,i) = nzeta
      r_bas(1,l,i) = rc_l(l)
    enddo

    if(ch_temp(3:3) == 'P') then
      if(la+1 <= LPOLMAX .and. la+1 <= mxdl) then
        lmax_bas(i) = la+1
        n_bas(la+1,i) = nzeta
        r_bas(1,la+1,i) = rc_pol
      else
        write(iowrite,'(/,"   No polarization orbitals (l = ",i2,          &
            &   ") in basis set ",i3)') la+1, i
      endif
    endif

!   second zeta is the relaxed energy derivative of the first (r_bas = 0)

  enddo

  write(iowrite,'(/,"   Basis cutoff radii (bohr)",/)')
  write(iowrite,'("    l    r_siesta      r_c")')
  do l = 0,la
    write(iowrite,'(i5,2f11.4)') l, r_siesta(l), rc_l(l)
  enddo
  if(la+1 <= LPOLMAX) then
    write(iowrite,'(i5,11x,f11.4,"   polarization")') la+1, rc_pol
  endif

  return

end subroutine atom_kb_basis_rule
