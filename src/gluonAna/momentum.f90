!! momentum calculations module
module FLUE_mom
   use FLUE_constants, only: WP, pi
   implicit none

contains

   pure subroutine get_qhat(coord, shape, qhat, a)
    !!Calculates q^_\mu = (2*pi/a) * coord[i] / coordMax[i]
      ! Arguments
      !integer, dimension(4), intent(in) :: coord
      real(kind=WP), dimension(4), intent(in) :: coord
    !! The coordinate we are at
      integer, dimension(4), intent(in) :: shape
    !! The lattice dims
      real(kind=WP), optional, intent(in) :: a
    !! lattice spacing
      real(kind=WP) :: latSpace
      real(kind=WP), dimension(4), intent(out) :: qhat
    !! output q^_\mu
      latSpace = 1.0_wp
      if (present(a)) latSpace = a
      qhat = (2.0_wp * pi / a) * real(coord, kind=WP) / real(shape, kind=WP)
   end subroutine get_qhat

end module FLUE_mom
