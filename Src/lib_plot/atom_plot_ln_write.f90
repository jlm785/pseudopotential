!>  Writes a file with the logarithmic derivatives for gnuplot and matplotlib.
!>  Points where the calculation failed (stored as zero) are written
!>  as NaN, so the plotted line has a gap instead of a dip to zero.
!>
!>  \author       claude+J.L.Martins
!>  \version      6.1.0
!>  \date         25 October 2026.
!>  \copyright    GNU Public License v2

subroutine atom_plot_ln_write(ioplot, fname, c_val, nsc, numtotal,      &
       ehist, dlog, mxdsc, mxdpts)

! extracted from atom_plot_ln_outpdf and atom_plot_ln_outdis. 25 October 2026. claude

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdsc                            !<  dimension of maximum scattering channels
  integer, intent(in)               ::  mxdpts                           !<  dimension of number of histogram points

  integer, intent(in)               ::  ioplot                           !<  tape number
  character(len=*), intent(in)      ::  fname                            !<  file name

  integer, intent(in)               ::  nsc                              !<  number of scattering channels
  character(len=4), intent(in)      ::  c_val(mxdsc)                     !<  column labels

  integer, intent(in)               ::  numtotal                         !<  number of energies
  real(REAL64), intent(in)          ::  ehist(mxdpts)                    !<  energy
  real(REAL64), intent(in)          ::  dlog(mxdpts,mxdsc)               !<  log derivatives (zero: failed)

! constants

  real(REAL64), parameter    ::  ZERO = 0.0_REAL64

! counters

  integer                  :: i, k


  open (unit=ioplot,file=fname,status='unknown',form='formatted')

  write(ioplot,'("#    E(Ry)",10x,20(a4,11x))') (c_val(k),k=1,nsc)
  do i = 1,numtotal
    write(ioplot,'(2x,e13.6)', advance='no') ehist(i)
    do k = 1,nsc
!     exactly zero is the marker of a failed point
      if(dlog(i,k) == ZERO) then
        write(ioplot,'(2x,a13)', advance='no') '          NaN'
      else
        write(ioplot,'(2x,e13.6)', advance='no') dlog(i,k)
      endif
    enddo
    write(ioplot,'(a)') ''
  enddo
  close(ioplot)

  return

end subroutine atom_plot_ln_write
