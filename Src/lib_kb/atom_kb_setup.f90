!>  Modifies the parameters for the non-local KB pseudopotential.
!>
!>  \author       J.L.Martins
!>  \version      6.0.2
!>  \date         30 August 2021, 18 September 2021.
!>  \copyright    GNU Public License v2

subroutine atom_kb_setup(tbasis, llocal, lmax_pot, nql, delql,           &
         lmax_bas, n_bas, r_bas, nz_bas, r_siesta, r_99,                 &
         mxdl)

! mxdl, 18 September 2021. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdl                             !<  array dimension for basis angular momentum

  integer,intent(in)                ::  lmax_pot                         !<  maximum angular momentum for potential

  real(REAL64), intent(in)          ::  r_siesta(0:mxdl)                 !<  cutoff using the SIESTA recipe
  real(REAL64), intent(in)          ::  r_99(0:mxdl)                     !<  radius with 99% of charge

! input and output

  character(len=3), intent(inout)   ::  tbasis                           !<  type of basis

  integer, intent(inout)            ::  llocal                           !<  angular momentum for local potential (negative: maximum of l-dependent)
  integer, intent(inout)            ::  nql                              !<  number of points for the Fourier grid of local potential and densities
  real(REAL64), intent(inout)       ::  delql                            !<  spacing of points in Fourier grid

  integer,intent(inout)             ::  lmax_bas                         !<  maximum angular momentum in basis
  integer,intent(inout)             ::  n_bas(0:mxdl)                    !<  basis functions for angular momentum l
  real(REAL64), intent(inout)       ::  r_bas(3,0:mxdl)                  !<  cutoff for the basis (up to triple zeta and l = 4)
  integer, intent(inout)            ::  nz_bas(3,0:mxdl)                 !<  number of non-trivial zeroes in basis function

! local variables

  character(len=1)        ::  yesno
  integer                 ::  n_in
  real(REAL64)            ::  r_in

! constants

  real(REAL64), parameter           ::  ZERO = 0.0_REAL64

! counters

  integer                           ::  i, l


  write(6,*)
  write(6,*) '  The current parameters for the Kleiman-Bylander'
  write(6,*) '  non-local pseudopotential transformation are:'
  write(6,*)
  if(llocal < 0) then
    write(6,*) '  The maximum of all potentials will be used as local'
  else
    write(6,*) '  The local potential will be  l = ',llocal
  endif
  write(6,*)
  write(6,'("   The Fourier grid will range from q = 0.0000 to",         &
       &  f10.4, "  with steps of ",f10.4)') nql*delql, delql
  write(6,'("   Enough for a cutoff of ",f10.3," Ha.")')                 &
             (nql*delql)*(nql*delql)/8
  write(6,*)

  write(6,*)
  write(6,*) '  The basis type is: ',tbasis
  write(6,*)
  write(6,*) '  l       r_siesta       r_99            r_bas   n_nodes'
  write(6,*)
  do l = 0,lmax_bas
    write(6,'(i5,2(3x,f10.3),3x,10(3x,f10.3,2x,i2))') l, r_siesta(l),     &
              r_99(l), (r_bas(i,l),nz_bas(i,l),i=1,n_bas(l))
  enddo
  write(6,*)

  write(6,*)
  write(6,*) '  Do you want to change any of those parameters? (y/n)'
  write(6,*)

  read(5,*) yesno

  if(yesno == 'Y' .or. yesno == 'y') then

    write(6,*)
    write(6,*) '  Do you want to change the local potential? (y/n)'
    write(6,*)

    read(5,*) yesno

    if(yesno == 'Y' .or. yesno == 'y') then

      write(6,*)
      write(6,*) '  Enter the angular momentum for local potential'
      write(6,*) '  between 0 and ',lmax_pot,'.  Enter negative for maximum'
      write(6,*)

      read(5,*) n_in

      if(n_in > lmax_pot) then
        llocal = lmax_pot
        write(6,*)
        write(6,*) '  Wrong value. Will use l = ',llocal
        write(6,*)
      else
        llocal = n_in
      endif

    endif

    write(6,*)
    write(6,*) '  Do you want to change the Fourier grid? (y/n)'
    write(6,*)

    read(5,*) yesno

    if(yesno == 'Y' .or. yesno == 'y') then

      write(6,*)
      write(6,*) '  Enter the number of steps and step value (nql, delql)'
      write(6,*)

      read(5,*) n_in, r_in

      if(n_in < 0) then
        write(6,*)
        write(6,*) '  Wrong value. Will use nql = ',nql
        write(6,*)
      else
        if(n_in < 100 .or. n_in > 50000) write(6,*) '  Strange value of nql...'
        nql = n_in
      endif

      if(r_in < ZERO) then
        write(6,*)
        write(6,*) '  Wrong value. Will use delql = ',delql
        write(6,*)
      else
        if(r_in < 0.0001 .or. r_in > 0.5) write(6,*) '  Strange value of delql...'
        delql = r_in
      endif
    endif

    write(6,*)
    write(6,*) '  Do you want to change the basis parameters? (y/n)'
    write(6,*)

    read(5,*) yesno

    if(yesno == 'Y' .or. yesno == 'y') then

      write(6,*)
      write(6,'("  Enter the maximum angular momentum (0,...,",i2,")")') mxdl
      write(6,*)

      read(5,*) n_in

      if(n_in < 0 .or. n_in > mxdl) then

        write(6,*)
        write(6,'("  Enter the maximum angular momentum (0,...,",i2,     &
            &    ") again")') mxdl
        write(6,*)
        read(5,*) n_in
        if(n_in < 0 .or. n_in > mxdl) then
          write(6,*)
          write(6,*) '  Wrong value using max l = ',lmax_bas
          write(6,*)
        else
          lmax_bas = n_in
        endif

      else

        lmax_bas = n_in

      endif

      do l = 0,lmax_bas
        write(6,*)
        write(6,'("  Enter the number of basis functions for l = ",      &
            &   i2,"(1,2,3)")') l
        write(6,*)

        read(5,*) n_in

        if(n_in < 1 .or. n_in > 3) then
          write(6,*)
          write(6,*) '  Enter the number of basis functions (1,2,3) again'
          write(6,*)
          read(5,*) n_in
          if(n_in < 1 .or. n_in > 3) then
            write(6,*)
            write(6,*) '  Wrong value using ',n_bas(l),' functions'
            write(6,*)
          else
            n_bas(l) = n_in
          endif
        else
          n_bas(l) = n_in
        endif

        do i = 1,n_bas(l)
          write(6,*)
          write(6,'("  Enter the cutoff radius for basis functions l = ",    &
            &    i2," n = ",i3)') l, i
          write(6,*)
          read(5,*) r_in
          if(r_in < 0) then
            write(6,*)
            write(6,*) '  Enter the cutoff radius again (> 0) '
            write(6,*)
            read(5,*) r_in
            if(r_in < 0) then
              write(6,*)
              write(6,*) '  using previous value ', r_bas(i,l)
              write(6,*)
            else
              if(r_in < 1.0 .or. r_in > 15.0) write(6,*) '  Strange value...'
              r_bas(i,l) = r_in
            endif
          else
            if(r_in < 1.0 .or. r_in > 15.0) write(6,*) '  Strange value...'
            r_bas(i,l) = r_in
          endif

        enddo
      enddo

      write(6,*)
      write(6,*) '  Do you want to change the number of non-trivial nodes? (y/n)'
      write(6,*)

      read(5,*) yesno

      if(yesno == 'Y' .or. yesno == 'y') then

        do l = 0,lmax_bas
          do i = 1,n_bas(l)
            write(6,*)
            write(6,*) '  Enter the number of non-trivial nodes for ',   &
                 'basis function l = ',l,' n = ',i
            write(6,*)
            read(5,*) n_in

            if(n_in < 0) then
              write(6,*)
              write(6,*) '  Wrong value. Will use nz_bas = 0'
              write(6,*)
              nz_bas(i,l) = 0
            elseif(n_in - (i-1) > 0) then
              write(6,*)
              write(6,*) '  Strange value. Enter it again'
              write(6,*)
              read(5,*) n_in
              nz_bas(i,l) = n_in
            else
              nz_bas(i,l) = n_in
            endif
          enddo
        enddo
      endif

    endif


  else

    write(6,*)
    write(6,*) '  Continuing with the above parameters'
    write(6,*)

  endif


  return

end subroutine atom_kb_setup
