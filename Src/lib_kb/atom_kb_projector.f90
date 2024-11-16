!>  Calculates the KB projectors from the semi-local potentials
!>  See L.Kleinman and D.M.Bylander, Phys.Rev.Lett. 48, 1425 (1982).
!>
!>  \author       Norm troullier, J.L.Martins
!>  \version      6.0.9
!>  \date         early 90s, May 2012, 18 September 2021.
!>  \copyright    GNU Public License v2

subroutine atom_kb_projector(npot, lo, nr, r, drdi, irel, llocal,        &
      vlocal, vionic, psi,  inorm, vkbproj,                              &
      iowrite, mxdl, mxdnr)


! adapted from the old program jlm 22/5/2012
!
! The projector is |vkb Ylm>(inorm)<Ylm vkb| with inorm =-1,0,1
! bug in inorm(l=0,j=1) corrected 25/5/2014
! mxdl, mxdnr, iowrite, 20 September 2021. JLM
! Initialize output to zero (not needed). 16 November 2024. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum components
  integer, intent(in)               ::  mxdnr                            !<  dimension of radial grid points.

  integer, intent(in)               ::  iowrite                          !<  default tape for writing

  integer, intent(in)               ::  npot(-1:1)                       !<  number of orbitals to be calculated
  integer, intent(in)               ::  lo(mxdl+1,-1:1)                  !<  angular momentum of orbital
  integer, intent(in)               ::  nr                               !<  number of points in the radial grid
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial grid points   r(i) = a*(exp(b*i)-1), i=1,...,nr
  real(REAL64), intent(in)          ::  drdi(mxdnr)                      !<  d r(i) / d i
  character(len=3), intent(in)      ::  irel                             !<  flag for relativistic (r) or spin-polarized (s) original calculation
  real(REAL64), intent(in)          ::  vionic(mxdnr,0:mxdl,-1:1)        !<  r*pseudopotential
  integer, intent(in)               ::  llocal                           !<  choice of local potential                                                      !
  real(REAL64), intent(in)          ::  vlocal(mxdnr)                    !<  r*local pseudopotential
  real(REAL64), intent(in)          ::  psi(mxdnr,0:mxdl,-1:1)           !<  wavefunctions (r(i),l,2j-2l).

! output

  integer, intent(out)              ::  inorm(0:mxdl,-1:1)               !<  sign of denominator of KB operator
  real(REAL64), intent(out)         ::  vkbproj(mxdnr,0:mxdl,-1:1)       !<  kb-projector

! allocatable arrays

  real(REAL64), allocatable         ::  w(:)

! local variables

  integer                    ::  l, jmax
  real(REAL64)               ::  anorm, a2norm, ratio

! parameters

  real(REAL64), parameter    ::  ONE = 1.0_REAL64, ZERO = 0.0_REAL64
  real(REAL64), parameter    ::  EPS = 1.0E-8_REAL64

! counters

  integer                    ::  i, j, k


  write(iowrite,*)
  write(iowrite,'("   Constructing the projectors")')
  write(iowrite,*)

  allocate(w(nr))

  inorm(:,:) = 0
  vkbproj(:,:,:) = ZERO

  jmax=0
  if(irel == 'rel') jmax = 1
  do j = -jmax,jmax
    do i = 1,npot(j)

      l = lo(i,j)
      do k = 2,30
        vkbproj(k,l,j) = vionic(k,l,j)/r(k) - vlocal(k)/r(k)
      enddo

      do k = 31,nr
        vkbproj(k,l,j) = (vionic(k,l,j)-vlocal(k)) / r(k)
      enddo
      vkbproj(1,l,j) = vkbproj(2,l,j)

!     calculate the normalizing coefficient for nonlocal parts,
!     uses bode's rule for integration.

      do k = 1,nr
        w(k) = psi(k,l,j)*psi(k,l,j)*vkbproj(k,l,j)
      enddo

      call atom_atm_integ(anorm,w,drdi,nr)

      do k = 1,nr
        w(k) = abs(w(k))
      enddo

      call atom_atm_integ(a2norm,w,drdi,nr)


!jlm
      if (abs(anorm) < sqrt(EPS)*a2norm .and. llocal /= l .and. a2norm > EPS) then
        write(iowrite,*)
        write(iowrite,'(" WARNING in atom_kb_projector: small anorm")')
        write(iowrite,'(" anorm, a2norm = ",2e15.7)') anorm, a2norm
        write(iowrite,*)
      endif
      if (llocal == l .and. (j == 0 .or. l == 0)) then
!jlm
        inorm(l,j) = 0
        anorm = ONE
      else
        if (anorm < ZERO) then
          inorm(l,j) = -1
          anorm = sqrt(-anorm)
        elseif(anorm > ZERO) then
          inorm(l,j) = 1
          anorm = sqrt(anorm)
        else
          if(a2norm > EPS) then
            write(iowrite,*)
            write(iowrite,'("     EXPECT DISASTER   ")')
            write(iowrite,*)
          endif
          inorm(l,j) = 0
          anorm = ONE
        endif
      endif


      if (inorm(l,j) /= 0) then

        if ( a2norm /= ZERO) then
          ratio = anorm*anorm/a2norm
          if(j == 0) then
            write(iowrite,'(1x,"  ratio for l= ",i3,9x," is ",f12.4)')   &
                l,inorm(l,j)*ratio
          else
            write(iowrite,'(1x,"  ratio for l= ",i3," j= ",i3,"/2 is ",f12.4)')   &
                 l,2*l+j,inorm(l,j)*ratio
          endif
        endif
      endif

!     calculate projector, note that is stored as psi(r)deltav(r).

      if (inorm(l,j) == 0) then
        do k = 1,nr
          vkbproj(k,l,j) = ZERO
        enddo
      else
        do k = 2,nr
          vkbproj(k,l,j) = (psi(k,l,j)/r(k)) * vkbproj(k,l,j) / anorm
        enddo
        if(l == 0) then
          vkbproj(1,l,j) = vkbproj(2,l,j)
        else
          vkbproj(1,l,j) = ZERO
        endif
      endif

    enddo             !     end loop over orbitals

  enddo               !     end loop over spin-orbit

  write(iowrite,*)

  deallocate(w)

  return

end subroutine atom_kb_projector
