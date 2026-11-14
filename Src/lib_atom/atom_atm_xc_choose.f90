!>  Lets the user choose (interactive) the exchange and correlation
!>  functional written in atom.dat.
!>  The default (empty line or end of input) is the Ceperley-Alder LDA.
!>
!>  \author       claude
!>  \version      6.2.0
!>  \date         14 November 2026.
!>  \copyright    GNU Public License v2

subroutine atom_atm_xc_choose(icorr)

! Written 14 November 2026. JLM+claude

  implicit none

  integer, parameter          ::  NXC = 3                                !  number of functionals

! output

  character(len=2), intent(out)     ::  icorr                            !<  exchange and correlation functional (ca, pw, pb)

! local variables

  character(len=2)                  ::  xcode(NXC)                       !  code in atom.dat
  character(len=3)                  ::  xtype(NXC)                       !  LDA or GGA
  character(len=48)                 ::  xname(NXC)                       !  description

  character(len=80)                 ::  line
  integer                           ::  ichoice, ios

! counters

  integer     ::  k


  xcode(1) = 'ca'
  xtype(1) = 'LDA'
  xname(1) = 'Ceperley-Alder, Perdew-Zunger parametrization'

  xcode(2) = 'pw'
  xtype(2) = 'LDA'
  xname(2) = 'Perdew-Wang 92'

  xcode(3) = 'pb'
  xtype(3) = 'GGA'
  xname(3) = 'Perdew-Burke-Ernzerhof (PBE)'

  write(6,*)
  write(6,'("  Exchange and correlation functional:")')
  write(6,*)
  do k = 1,NXC
    write(6,'(4x,i2,3x,a3,3x,a2,3x,a)') k, xtype(k), xcode(k),            &
         trim(xname(k))
  enddo
  write(6,*)
  write(6,'("  Enter your choice (default 1)")')
  write(6,*)

  ichoice = 1
  read(5,'(a)',iostat=ios) line
  if(ios == 0) then
    read(line,*,iostat=ios) ichoice
    if(ios /= 0) ichoice = 1
  endif

  if(ichoice < 1 .or. ichoice > NXC) then
    write(6,'("  Choice ",i5," out of range, using 1")') ichoice
    ichoice = 1
  endif

  icorr = xcode(ichoice)

  return

end subroutine atom_atm_xc_choose
