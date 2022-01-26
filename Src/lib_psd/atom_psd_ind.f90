!>  Calculates indv and indx that indicates which orbital should be used
!>  to calculate the pseudopotentials for a given angular momentum, spin or j.
!>
!>  \author       Sverre Froyen,  Norm Troullier, Jose Luis Martins
!>  \version      6.0.6
!>  \date         1980s and 1990s, 8 July 2021, 2 November 2021.
!>  \copyright    GNU Public License v2

subroutine atom_psd_ind(indv, indx, nindx, lcmax, norb, ncore, lo, iso,  &
        mxdorb, mxdlc, mxdnw)

! extracted from pseud2. 8 July 2021. JLM
! indd, indu -> indv, indx. 15 September 2021. JLM
! new structure indx. 2 November 2021. JLM

 implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdorb                           !<  dimension of the number of orbitals
  integer, intent(in)               ::  mxdlc                            !<  dimension of scattering channel ang. mom.
  integer, intent(in)               ::  mxdnw                            !<  dimension of number of wave-functions same l

  integer, intent(in)               ::  norb                             !<  number of orbitals
  integer, intent(in)               ::  ncore                            !<  number of orbitals treated as core
  integer, intent(in)               ::  lo(mxdorb)                       !<  angular quantum number l
  integer, intent(in)               ::  iso(mxdorb)                      !<  2*spin or 2*(j-l)

! output

  integer, intent(out)              ::  indv(0:mxdlc,-1:1)               !<  main orbital for scattering channel l (legacy)
  integer, intent(out)              ::  indx(mxdnw,0:mxdlc,-1:1)         !<  orbitals for scattering channel l
  integer, intent(out)              ::  nindx(0:mxdlc,-1:1)              !<  number of orbitals for scattering channel l
  integer, intent(out)              ::  lcmax                            !<  maximum ang mom in scattering channels

! counters

  integer        ::  i, l

  do i = -1,1
  do l = 0,mxdlc
    indv(l,i) = 0
    nindx(l,i) = 0
    indx(1,l,i) = 0
    indx(2,l,i) = 0
  enddo
  enddo

  lcmax = 0
  do i = ncore+1,norb

    if(lcmax < lo(i)) lcmax = lo(i)

    nindx(lo(i),iso(i)) = nindx(lo(i),iso(i)) + 1

    if(nindx(lo(i),iso(i)) == 1) then
      indv(lo(i),iso(i)) = i
    endif
    if(nindx(lo(i),iso(i)) < mxdnw+1) then
      indx(nindx(lo(i),iso(i)),lo(i),iso(i)) = i
    else
      write(6,*)
      write(6,'("  warning in atom_psd_ind - more than ",i3,             &
                    " valence orbitals")') mxdnw
      write(6,'("  with angular momentum ",i2," and 2*spin ",i2," exist")') lo(i), iso(i)
      write(6,*)
    endif

  enddo

  return

end subroutine atom_psd_ind

