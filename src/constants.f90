!! Constants module
module FLUE_constants
   use, intrinsic :: iso_c_binding, only: c_double, c_double_complex
   implicit none

   integer, parameter :: WP = c_double
   integer, parameter :: WC = c_double
   !real(WP) :: pi = 3.1415926535897932384626433_WP
   real(kind=WP), parameter :: PI = acos(-1.0_wp)
end module FLUE_constants
