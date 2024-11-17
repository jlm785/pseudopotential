!>  Calculates several quantities for the Quantum Espresso UPF format
!>
!>  \author       J.L.Martins
!>  \version      6.0.9
!>  \date         17 November 2024.
!>  \copyright    GNU Public License v2

subroutine atom_kb_psd_out_upf_config(lmax_pot, psdtitle, lso,           &
      llocal, inorm,                                                     &
      n_conf, l_skip, o_conf, r_core,                                    &
      nproj, l_vnl, is_vnl, nchi, l_chi, is_chi,                         &
      mxdl)

! Extracted from the new atom_kb_psd_out_upf. 17 November 2024. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! dimensions:

  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum components

! input:

  integer, intent(in)               ::  lmax_pot                         !<  maximum angular momentum in potential
  character(len=10), intent(in)     ::  psdtitle(20)                     !<  pseudopotential parameters
  logical, intent(in)               ::  lso                              !<  indicates if spin-orbit is used (irel == 'rel')

  integer, intent(in)               ::  llocal                           !<  angular momentum of local pseudopotential (<0 supremum)
  integer, intent(in)               ::  inorm(0:mxdl,-1:1)               !<  sign of denominator of KB operator

! output

  integer, intent(out)              ::  n_conf(0:lmax_pot)               !<  principal quantum number
  logical, intent(out)              ::  l_skip(0:lmax_pot)               !<  skip as information is missing
  real(REAL64), intent(out)         ::  o_conf(0:lmax_pot,2)             !<  occupation of orbital
  real(REAL64), intent(out)         ::  r_core(0:lmax_pot)               !<  core radius

  integer, intent(out)              ::  nproj                            !<  number of projectors
  integer, intent(out)              ::  nchi                             !<  number of pseudo wave-functions

  integer, intent(out)              ::  l_vnl(2*lmax_pot+1)              !<  angular momentum for projector i
  integer, intent(out)              ::  is_vnl(2*lmax_pot+1)             !<  2(j-l) for projector i

  integer, intent(out)              ::  l_chi(2*lmax_pot+1)              !<  angular momentum for wave-function i
  integer, intent(out)              ::  is_chi(2*lmax_pot+1)             !<  2(j-l) for wave-function i

! local variables

  character(len=10)       ::  ctemp
  character(len=2)        ::  nl
  integer                 ::  ierr

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64, UM = 1.0_REAL64
  character(len=1), parameter  ::  IL(0:6) = (/'s','p','d','f','g','h','i'/)

! counters

  integer       ::   l, np


! recover information from psdtitle

  do l = 0,lmax_pot
    n_conf(l) = 0
    o_conf(l,1) = ZERO
    o_conf(l,2) = ZERO
    r_core(l) = ZERO
    ctemp = psdtitle(l*2+1)
    nl = '  '
    read(ctemp(1:2),'(a2)',iostat=ierr) nl
    if(nl == '  ' .or. ierr /= 0) then
      l_skip(l) = .TRUE.
    else
      l_skip(l) = .FALSE.
      read(ctemp(1:1),'(i1)',iostat=ierr) n_conf(l)
      if(ctemp(2:2) /= IL(l)) then
        write(6,*) '   inconsitent psdtitle, l, IL, n_conf',l,IL(l),ctemp(2:2)
      endif
      read(ctemp(4:9),*,iostat=ierr) o_conf(l,1)
      ctemp = psdtitle(l*2+2)
      if(ctemp(1:1) == ' ' .or. ctemp(1:1) == 'r') then
        read(ctemp(6:10),*,iostat=ierr) r_core(l)
      else
        read(ctemp(1:4),*,iostat=ierr) o_conf(l,2)
        read(ctemp(7:10),*,iostat=ierr) r_core(l)
      endif

    endif
  enddo


! attributes of projector  (the first part is just paranoia)

  if(lso) then

    if(llocal == 0) then
      nproj = 2*lmax_pot
    else
      nproj = 2*lmax_pot+1
    endif

    np = 0
    if(inorm(0, 1) /= 0) then
      np = 1
    endif

    if(lmax_pot > 0) then
      do l = 1,lmax_pot
        if(inorm(l,-1) /= 0) then
          np = np + 1
        endif
        if(inorm(l, 1) /= 0) then
          np = np + 1
        endif

      enddo
    endif

  else

    if(llocal < 0) then
      nproj = lmax_pot+1
    else
      nproj = lmax_pot
    endif

    np = 0
    do l = 0,lmax_pot
      if(inorm(l,0) /= 0) then
        np = np + 1
      endif
    enddo

  endif

  if(np /= nproj) then
    write(6,*) '  stopped in kb_psd_out_upf'
    write(6,*) '  inconsistent number of projectors',np,nproj

    STOP

  endif


  if(lso) then

    np = 0
    if(inorm(0, 1) /= 0) then
      np = np + 1
      l_vnl(np) = 0
      is_vnl(np) = 1
    endif

    if(lmax_pot > 0) then
      do l = 1,lmax_pot
        if(inorm(l,-1) /= 0) then
          np = np + 1
          l_vnl(np) = l
          is_vnl(np) =-1
        endif
        if(inorm(l, 1) /= 0) then
          np = np + 1
          l_vnl(np) = l
          is_vnl(np) = 1
        endif

      enddo
    endif

  else

    np = 0
    do l = 0,lmax_pot
      if(inorm(l,0) /= 0) then
        np = np + 1
        l_vnl(np) = l
        is_vnl(np) = 0
      endif
    enddo

  endif

! attributes of the wave-functions

  if(lso) then
    nchi = 2*lmax_pot+1
  else
    nchi = lmax_pot+1
  endif

  if(lso) then

    l_chi(1) = 0
    is_chi(1) = 1
    if(lmax_pot > 0) then
      do l = 1,lmax_pot
        l_chi(2*l  ) = l
        l_chi(2*l+1) = l
        is_chi(2*l  ) = -1
        is_chi(2*l+1) =  1
      enddo
    endif

  else

    do l = 0,lmax_pot
      l_chi(l+1) = l
      is_chi(l+1) = 0
    enddo

  endif

  return

end subroutine atom_kb_psd_out_upf_config
