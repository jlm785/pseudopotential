!>  Prints one set of xy data to itplot and identifies it with a marker,
!>  as atom_plot_one but in exponential format for data spanning many
!>  orders of magnitude (densities close to the nucleus, Fourier tails),
!>  and without angular momentum in the marker.
!>
!>  \author       Jose Luis Martins, claude
!>  \version      6.2.0
!>  \date         16 November 2026.
!>  \copyright    GNU Public License v2

subroutine atom_plot_one_es(nfirst, nlast, r, v, rstep, marker,          &
            itplot, mxdnr)

! Written from atom_plot_one. 16 November 2026. claude+JLM

  implicit none

  integer, parameter          :: REAL64 = selected_real_kind(12)

! input

  integer, intent(in)               ::  mxdnr                            !<  dimension  of the number of points
  integer, intent(in)               ::  itplot                           !<  default plot file

  integer, intent(in)               ::  nfirst, nlast                    !<  first and last points to plot
  real(REAL64), intent(in)          ::  r(mxdnr)                         !<  radial or reciprocal grid
  real(REAL64), intent(in)          ::  v(mxdnr)                         !<  function on the grid

  real(REAL64), intent(in)          ::  rstep                            !<  approximate step for plot
  character(len=*), intent(in)      ::  marker                           !<  marker (at most 10 characters)

! local variables

  real(REAL64)           ::  step
  character(len=17)      ::  markline                                    !  marker line, left justified

! constants

  real(REAL64), parameter       ::  ZERO = 0.0_REAL64

! counter

  integer             ::  j


! the first point is always written (it may be very close to r = 0)

  step = ZERO
  do j = nfirst,nlast
    if (r(j) >= step) then
      write(itplot,'(1x,es18.10,3x,es20.12)') r(j), v(j)
      step = r(j) + rstep
    endif
  enddo

  markline = 'marker '//trim(adjustl(marker))
  write(itplot,'(1x,a17)') markline

  return

end subroutine atom_plot_one_es
