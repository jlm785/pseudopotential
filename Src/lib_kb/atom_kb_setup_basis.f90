!>  Modifies the parameters for the atomic basis sets
!>
!>  \author       J.L.Martins+claude
!>  \version      6.1.0
!>  \date         14 January 2026.
!>  \copyright    GNU Public License v2

subroutine atom_kb_setup_basis(n_bsets, tbasis,                          &
         lmax_bas, n_bas, r_bas, nz_bas, r_siesta, r_99,                 &
         mxdl, mxdset)

! split from atom_kb_setup. 14 January 2026. JLM
! up to 4 zeta, r_bas = 0 for relaxed energy derivative of previous function. 24 October 2026. JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdl                             !<  array dimension for basis angular momentum
  integer, intent(in)               ::  mxdset                           !<  dimension for number of atomic basis sets

  real(REAL64), intent(in)          ::  r_siesta(0:mxdl)                 !<  cutoff using the SIESTA recipe
  real(REAL64), intent(in)          ::  r_99(0:mxdl)                     !<  radius with 99% of charge

! input and output

  integer, intent(inout)            ::  n_bsets                          !<  number of atomic basis sets
  character(len=3), intent(inout)   ::  tbasis(mxdset)                   !<  type of basis
  integer,intent(inout)             ::  lmax_bas(mxdset)                 !<  maximum angular momentum in basis

  integer,intent(inout)             ::  n_bas(0:mxdl,mxdset)             !<  basis functions for angular momentum l
  real(REAL64), intent(inout)       ::  r_bas(4,0:mxdl,mxdset)           !<  cutoff for the basis (up to quadruple zeta), 0: energy derivative of previous
  integer, intent(inout)            ::  nz_bas(4,0:mxdl,mxdset)          !<  number of non-trivial zeroes in basis function

! local variables

  character(len=1)        ::  yesno
  integer                 ::  n_in
  real(REAL64)            ::  r_in

  integer                 ::  n_bsets_old

  logical                 ::  lvalid                                     !  valid cutoff radius
  integer                 ::  ntry

! constants

  real(REAL64), parameter           ::  ZERO = 0.0_REAL64
  real(REAL64), parameter           ::  SMALL = 1.0E-6_REAL64
  real(REAL64), parameter           ::  EPS = 1.0E-9_REAL64

  integer, parameter                ::  MXDZETA = 4                      !  maximum number of basis functions for each l

! counters

  integer                           ::  i, j, l

  write(6,*)
  write(6,'("   There are ",i4," default atomic basis sets")') n_bsets

  n_bsets_old = n_bsets

  write(6,*)
  write(6,*) '  Do you want to change the number of basis sets? (y/n)'
  write(6,*)

  read(5,*) yesno

  if(yesno == 'Y' .or. yesno == 'y') then

    write(6,*)
    write(6,'("  Enter the number of atomic basis sets (1,...,",i2,")")') n_bsets_old
    write(6,*)

    read(5,*) n_in

    if(n_in < 1 .or. n_in > mxdset) then

      write(6,*)
      write(6,*) "  wrong value, enter number of basis again"
      write(6,*)

      read(5,*) n_in

      if(n_in < 1 .or. n_in > mxdset) then
        write(6,*)
        write(6,*) "  wrong value again, using ", n_bsets_old
        write(6,*)

        n_in = n_bsets_old

      endif

    endif

    n_bsets = n_in

  endif

! loop over basis sets

  do j = 1,n_bsets

    yesno = 'Y'

    if(j <= n_bsets_old) then

      write(6,*)
      write(6,'("   The type of basis ",i4," is : ",a3)') j, tbasis(j)
      write(6,*)
      write(6,*) '  l       r_siesta       r_99            r_bas   n_nodes'
      write(6,*)
      do l = 0,lmax_bas(j)
        write(6,'(i5,2(3x,f10.3),3x)', advance='no') l, r_siesta(l), r_99(l)
        do i = 1,n_bas(l,j)
          if(r_bas(i,l,j) > EPS) then
            write(6,'(3x,f10.3,2x,i2)', advance='no') r_bas(i,l,j), nz_bas(i,l,j)
          else
            write(6,'(a)', advance='no') '          phi    '
          endif
        enddo
        write(6,*)
      enddo
      write(6,*)
      write(6,*) '  phi: relaxed energy derivative of the previous function'
      write(6,*)

      write(6,*)
      write(6,*) '  Do you want to change this parameters? (y/n)'
      write(6,*)

      read(5,*) yesno

    else

      write(6,*)
      write(6,*) '   New basis set'
      write(6,*)

    endif

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
          write(6,*) '  Wrong value using max l = ',lmax_bas(j)
          write(6,*)
        else
          lmax_bas(j) = n_in
        endif

      else

        lmax_bas(j) = n_in

      endif

      do l = 0,lmax_bas(j)
        write(6,*)
        write(6,'("  Enter the number of basis functions for l = ",      &
            &   i2," (1,...,",i2,")")') l, MXDZETA
        write(6,*)

        read(5,*) n_in

        if(n_in < 1 .or. n_in > MXDZETA) then
          write(6,*)
          write(6,'("  Enter the number of basis functions (1,...,",i2,   &
              &   ") again")') MXDZETA
          write(6,*)
          read(5,*) n_in
          if(n_in < 1 .or. n_in > MXDZETA) then
            write(6,*)
            write(6,*) '  Wrong value using ',n_bas(l,j),' functions'
            write(6,*)
          else
            n_bas(l,j) = n_in
          endif
        else
          n_bas(l,j) = n_in
        endif

!       a cutoff radius of 0 asks for the relaxed energy derivative (phi)
!       of the previous function, which must not be itself a derivative.
!       The recommended double zeta is a cutoff followed by phi.

        do i = 1,n_bas(l,j)

          lvalid = .FALSE.
          do ntry = 1,2
            write(6,*)
            write(6,'("  Enter the cutoff radius for basis functions l = ",  &
              &    i2," n = ",i3)') l, i
            if(i > 1) then
              if(r_bas(i-1,l,j) > EPS) then
                write(6,*) '  (or 0 for the energy derivative of the ',     &
                           'previous function, recommended)'
              endif
            endif
            write(6,*)
            read(5,*) r_in

            if(r_in > SMALL) then
              lvalid = .TRUE.
              if(r_in < 1.0 .or. r_in > 15.0) write(6,*) '  Strange value...'
            elseif(abs(r_in) <= SMALL .and. i > 1) then
              if(r_bas(i-1,l,j) > ZERO) then
                lvalid = .TRUE.
                r_in = ZERO
              endif
            endif

            if(lvalid) exit

            write(6,*)
            if(i == 1) then
              write(6,*) '  The cutoff radius must be > 0'
            else
              write(6,*) '  The cutoff radius must be > 0, or 0 if the ',   &
                         'previous function is not an energy derivative'
            endif
          enddo

          if(lvalid) then
            r_bas(i,l,j) = r_in
          else
            if(i == 1) then
              if(r_bas(i,l,j) <= ZERO) r_bas(i,l,j) = max(r_siesta(l),r_siesta(0))
            elseif(r_bas(i-1,l,j) > EPS) then
              r_bas(i,l,j) = ZERO
            else
              r_bas(i,l,j) = 1.2*r_bas(i-2,l,j)
            endif
            write(6,*)
            if(r_bas(i,l,j) > EPS) then
              write(6,*) '  Using the cutoff radius ', r_bas(i,l,j)
            else
              write(6,*) '  Using the energy derivative of the previous function'
            endif
            write(6,*)
          endif

          if(r_bas(i,l,j) <= EPS) nz_bas(i,l,j) = 0

        enddo
      enddo

      write(6,*)
      write(6,*) '  Do you want to change the number of non-trivial nodes? (y/n)'
      write(6,*)

      read(5,*) yesno

      if(yesno == 'Y' .or. yesno == 'y') then

        do l = 0,lmax_bas(j)
          do i = 1,n_bas(l,j)

!           not used for energy derivatives

            if(r_bas(i,l,j) > EPS) then

              write(6,*)
              write(6,*) '  Enter the number of non-trivial nodes for ',   &
                   'basis function l = ',l,' n = ',i
              write(6,*)
              read(5,*) n_in

              if(n_in < 0) then
                write(6,*)
                write(6,*) '  Wrong value. Will use nz_bas = 0'
                write(6,*)
                nz_bas(i,l,j) = 0
              elseif(n_in - (i-1) > 0) then
                write(6,*)
                write(6,*) '  Strange value. Enter it again'
                write(6,*)
                read(5,*) n_in
                nz_bas(i,l,j) = n_in
              else
                nz_bas(i,l,j) = n_in
              endif

            endif

          enddo
        enddo
      endif


    else

      write(6,*)
      write(6,*) '  Continuing with the above parameters'
      write(6,*)

    endif

  enddo


  return

end subroutine atom_kb_setup_basis
