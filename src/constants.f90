module FLUE_constants
  !! Constants module
  use, intrinsic :: ISO_C_BINDING, only: C_DOUBLE, C_DOUBLE_COMPLEX
  implicit none

  integer, parameter :: WP = C_DOUBLE
  integer, parameter :: WC = C_DOUBLE
  !real(WP) :: pi = 3.1415926535897932384626433_WP
  real(kind=WP), parameter :: PI = acos(-1.0_WP)
end module FLUE_constants
