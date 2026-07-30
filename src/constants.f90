!! Constants module
module FLUE_constants
!< Module FLUE_constants
!< Numerical constants and precision definitions for the FLUE codebase.
!<
!< Purpose:
!<   This module provides a centralised and portable definition of working
!<   precisions and commonly used mathematical constants using
!<   `ISO_C_BINDING`.
!<
!< Design:
!<   - `WP` defines the working real precision (double precision).
!<   - `WC` defines the corresponding complex precision.
!<   - `SP` provides single precision support where needed.
!<
!< Mathematical constants:
!<   - `PI` is defined via `acos(-1)` for portability.
!<   - `TWO_PI` is provided to avoid repeated computation and magic numbers.
!<
!< Example:
!<```
!< use FLUE_constants, only: WP, PI
!< real(kind=WP) :: theta
!< theta = 2.0_WP * PI
!<```
   use, intrinsic :: ISO_C_BINDING, only: C_DOUBLE, C_DOUBLE_COMPLEX, C_FLOAT, C_INT
   implicit none(type, external)
   public
   integer, parameter :: WP = C_DOUBLE
   !< Working precision for real-valued computations [double precision]
   integer, parameter :: WC = C_DOUBLE_COMPLEX
   !< Working precision for complex-valued computations [double precision complex].
   real(kind=WP), parameter :: PI = ACOS(-1.0_WP)
   !< Circular constant π ~ 3.14.....
   real(kind=WP), parameter :: TWO_PI = 2.0_WP * PI
   !< 2π.
   integer, parameter :: SP = C_FLOAT
   !< Single precision kind
end module FLUE_constants

