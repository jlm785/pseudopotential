!>  Writes the fourier transform of the KB pseudopotential in ascii format
!>  to file fname
!>
!>  \author       N. Troullier, J.L.Martins
!>  \version      6.0.8
!>  \date         November 1990, April-June 2012, 19 May 2022.
!>  \copyright    GNU Public License v2

subroutine atom_kb_psd_out_four(iotape, fname,                           &
      nameat, icorr, irel, nicore, irdate, irvers, irayps, psdtitle,     &
      nql, nqnl, delql, nqbas, delqbas, zion, vql0,                      &
      npot, lo, ev, inorm, vkbprft, vlocft, cdcft, cdvft,                &
      n_bsets, norbas, lo_b, basft,                                      &
      mxdl, nqmax, nqmaxbas, mxdbas, mxdset)

! modified by JLMartins 24/4/2012 and 12/6/2012
! wv -> bas.  18 September 2021. JLM
! psdtitle. 19 May 2022. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input:

  integer, intent(in)               ::  mxdl                             !<  dimension of angular momentum components
  integer, intent(in)               ::  nqmax                            !<  dimension of fourier grid points for potential
  integer, intent(in)               ::  nqmaxbas                         !<  dimension of fourier grid points.
  integer, intent(in)               ::  mxdbas                           !<  dimension of basis functions
  integer, intent(in)               ::  mxdset                           !<  dimension for number of atomic basis sets

  integer, intent(in)               ::  iotape                           !<  io tape number
  character(len=*), intent(in)      ::  fname                            !<  file name of the output of the atomic program, parsec style

  character(len=2), intent(in)      ::  nameat                           !<  chemical symbol of the element
  character(len=2), intent(in)      ::  icorr                            !<  correlation used in the calculation
  character(len=3), intent(in)      ::  irel                             !<  flag for relativistic (r) or spin-polarized (s) original calculation
  character(len=4), intent(in)      ::  nicore                           !<  flag for core correction
  character(len=10), intent(in)     ::  irdate, irvers                   !<  date and version of original calculation
  character(len=10), intent(in)     ::  irayps(4)                        !<  type of pseudopotential
  character(len=10), intent(in)     ::  psdtitle(20)                     !<  pseudopotential parameters

  integer, intent(in)               ::  npot(-1:1)                       !<  number of orbitals (s,p,d,...).  -1:  j=l-1/2.  0:  average.  1:  j=l+1/2
  integer, intent(in)               ::  nql                              !<  number of points for the local Fourier grid
  integer, intent(in)               ::  nqnl                             !<  number of points for the non-local Fourier grid
  real(REAL64), intent(in)          ::  delql                            !<  spacing of points in Fourier grid
  integer, intent(in)               ::  nqbas                            !<  number of points in the fourier grid
  real(REAL64) , intent(in)         ::  delqbas                          !<  spacing of the fourier grid

  real(REAL64), intent(in)          ::  zion                             !<  ionic charge

  integer, intent(in)               ::  lo(mxdl+1,-1:1)                  !<  angular momentum of orbital for each spin
  real(REAL64), intent(in)          ::  ev(0:mxdl,-1:1)                  !<  eigenvalues (l,2j-2l).

  real(REAL64), intent(in)          ::  vlocft(0:nqmax)                  !<  Fourier transform of local pseudopotential
  real(REAL64), intent(in)          ::  vql0                             !<  zero frequency of the integral without the Coulomb part

  integer, intent(in)               ::  inorm (0:mxdl,-1:1)              !<  sign of denominator of KB operator
  real(REAL64), intent(in)          ::  vkbprft(0:nqmax,0:mxdl,-1:1)     !<  Fourier transform of kb-projector
  real(REAL64), intent(in)          ::  cdcft(0:nqmax)                   !<  Fourier transform of 4*pi*r**2 charge density of core
  real(REAL64), intent(in)          ::  cdvft (0:nqmax)                  !<  Fourier transform of 4*pi*r**2 charge density of valence.  kb_conv is not spin polarized...

  integer, intent(in)               ::  n_bsets                          !<  number of atomic basis sets
  integer, intent(in)               ::  norbas(mxdset)                   !<  number of basis functions
  integer, intent(in)               ::  lo_b(mxdbas,mxdset)              !<  angular momentum of basis function
  real(REAL64), intent(in)          ::  basft(0:nqmaxbas,mxdbas,mxdset)  !< Fourier (Bessel) transform of the radial wavefunction

! local variables

  integer          ::  l

! constants

  real(REAL64), parameter   ::  HALF = 0.5_REAL64

! counters

  integer                        :: i, k, nb


  if(iotape == 0) RETURN

  open (unit = iotape, file = trim(fname), form='FORMATTED', status='UNKNOWN')

  write(iotape,'(1x,a2,1x,a2,1x,a3,1x,a4,1x,6a10,1x,20a10)')             &
       nameat,icorr,irel,nicore,irdate,irvers,irayps,psdtitle
  write(iotape,*) nint(zion),nql,delql,vql0

  if(irel == 'rel') then
    write(iotape,'(3i5)') npot(0),npot(-1),npot(1)
    write(iotape,'(20i5)') (lo(i,0),i=1,npot(0)),                        &
        (lo(i,-1),i=1,npot(-1)),(lo(i,1),i=1,npot(1))
    write(iotape,'(20i5)') (inorm(lo(i,0),0),i=1,npot(0)),               &
        (inorm(lo(i,-1),-1),i=1,npot(-1)),                               &
        (inorm(lo(i,1),1),i=1,npot(1))
    write(iotape,*) (ev(lo(i,0),0)*HALF,i=1,npot(0)),                    &
        (ev(lo(i,-1),-1)*HALF,i=1,npot(-1)),                             &
        (ev(lo(i,1),1)*HALF,i=1,npot(1))
  else
    write(iotape,'(i5)') npot(0)
    write(iotape,'(20i5)') (lo(i,0),i=1,npot(0))
    write(iotape,'(20i5)') (inorm(lo(i,0),0),i=1,npot(0))
    write(iotape,*) (ev(lo(i,0),0)*0.5d0,i=1,npot(0))
  endif

  do k=1,nql
    write(iotape,*) vlocft(k)
  enddo

  if(irel == 'rel') then
    do i = 1,npot(0)
      l = lo(i,0)
      if(l == 0) then
        do k = 0,nqnl
          write(iotape,*) (l*vkbprft(k,l,-1)+(l+1)*vkbprft(k,l,1))/(2*l+1)
        enddo
      else
        do k  =0,nqnl
          write(iotape,*) (l*vkbprft(k,l,-1)+(l+1)*vkbprft(k,l,1))/(2*l+1),  &
                      2*(vkbprft(k,l,1)- vkbprft(k,l,-1))/(2*l+1)
        enddo
      endif
    enddo
  else
    do i = 1,npot(0)
      l = lo(i,0)
      do k = 0,nqnl
        write(iotape,*) vkbprft(k,l,0)
      enddo
    enddo
  endif

  do k = 1,nql
    write(iotape,*) cdcft(k)
  enddo
  do k = 1,nql
    write(iotape,*) cdvft(k)
  enddo


  write(iotape,*) nqbas+1,delqbas,norbas(1),n_bsets
  do nb = 1,n_bsets
    do i = 1,norbas(nb)
      l = lo_b(i,nb)
      write(iotape,*) l, ev(l,0)/2, nb
      do k = 0,nqbas
        write(iotape,*) basft(k,i,nb)
      enddo
    enddo
  enddo

  close(unit=iotape)

  return

end subroutine atom_kb_psd_out_four
