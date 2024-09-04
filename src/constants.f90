module FLUE_constants
  !! Constants module
  use, intrinsic :: ISO_C_BINDING, only: C_DOUBLE
  implicit none

  integer, parameter :: WP = C_DOUBLE
  real(WP) :: pi = 3.1415926535897932384626433_WP
end module FLUE_constants
