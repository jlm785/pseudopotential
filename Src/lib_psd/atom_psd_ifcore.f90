!>  Finds the core correction recommendations
!>
!>  \author       Jose Luis Martins
!>  \version      6.0.3
!>  \date         23 September 2021.
!>  \copyright    GNU Public License v2

subroutine atom_psd_ifcore(ioread, filein, ioae, fileae,                 &
           ifcore_atm, ifcore_dat, ifcore_tab)


  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)


! input

  integer, intent(in)               ::  ioread                           !  default tape for reading
  character(len=*), intent(in)      ::  filein                           !  name of default tape for reading

  integer, intent(in)               ::  ioae                             !<  default tape for writing
  character(len=*), intent(in)      ::  fileae                           !<  name of default tape for all-electron results

! output

  integer, intent(out)              ::  ifcore_atm                       !<  0: no partial core correction; 1 : xc partial core correction; 2 : not recommended. From filein
  integer, intent(out)              ::  ifcore_dat                       !<  0: no partial core correction; 1 : xc partial core correction; 2 : not recommended. From fileae
  integer, intent(out)              ::  ifcore_tab                       !<  0: no partial core correction; 1 : xc partial core correction; 2 : not recommended. From p_tbl_psd_pcc

! local variables

  integer              ::  itype                                         !  type of calculation (-1 signals end of calculation)
  character(len=2)     ::  ctype                                         !  type of calculation flag
  character(len=2)     ::  pcc                                           !  no pcc: pg; pcc: pe.
  character(len=2)     ::  nameat, name                                  !  chemical symbot of the atom
  logical              ::  lex
  integer              ::  lmaxp1


  open (unit=ioae, file=trim(fileae), status='old', form='unformatted')

  read(ioae) itype
  ifcore_dat = itype - 1
  read(ioae)
  read(ioae)
  read(ioae) lmaxp1, nameat

  close (unit=ioae)


  inquire(file=filein, exist = lex)

  if(lex) then

    open(unit=ioread, file=trim(filein), status='OLD', form='FORMATTED')

    read(ioread,'(3x,a2)') ctype
    if(ctype /= 'ae') then
      read(ioread,*)
    endif
    read(ioread,'(3x,a2)') name

    ifcore_atm = 0
    if (ctype == 'pe') then
      ifcore_atm = 1
    elseif (ctype == 'ph') then
      ifcore_atm = 2
    endif

  else

    ifcore_atm =-1

  endif

  close(unit=ioread)


  call atom_p_tbl_psd_pcc(nameat, pcc, 0)

  ifcore_tab = 0
  if(pcc == 'pe') ifcore_tab = 1

  return

end subroutine atom_psd_ifcore
