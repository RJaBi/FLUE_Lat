!! Constants module
module FLUE_constants
   use, intrinsic :: ISO_C_BINDING, only: C_DOUBLE, C_DOUBLE_COMPLEX, C_FLOAT
   implicit none(external)
   public

   ! These are used
   integer, parameter :: WP = C_DOUBLE
   integer, parameter :: WC = C_DOUBLE
   !real(WP) :: pi = 3.1415926535897932384626433_WP
   real(kind=WP), parameter :: PI = ACOS(-1.0_WP)

   integer, parameter :: SP = C_FLOAT
   
end module FLUE_constants
