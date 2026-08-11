module FLUE_SU3_random
  !< Contains a (unused) function to construct an SU3 matrix
  !< From 3 SU2 matrices
  !< i.e. Eqn 4.31 of Gattringer Lang textbook
   use FLUE_constants, only: WP, WC
   use FLUE_SU2_random, only: constructSU2Matrix
   use stdlib_intrinsics, only: stdlib_matmul
   implicit none(type, external)
   private

   public :: constructSU3Matrix

contains

  pure function constructSU3Matrix(R, S, T) result(U)
    !< Constructs an SU3 matrix from 3 SU2 matrices
    !< eqn 4.31 of Gattringer Lang textbook
     complex(kind=WC), dimension(2, 2), intent(IN) :: R, S, T
     !< Input SU2 matrices
      complex(kind=WC), dimension(3, 3) :: U
      !< output SU3 matrix
      complex(kind=WC), dimension(3, 3) :: RE, SE, TE
      !< embedded SU2 matrices in 3x3 complex matrices
    !! in eqn 4.31 of Gattringer Lang
    !! embed r
      RE = (0.0_WP, 0.0_WP)
      RE(1, 1) = r(1, 1)
      RE(1, 2) = r(1, 2)
      RE(2, 1) = r(2, 1)
      RE(2, 2) = r(2, 2)
      RE(3, 3) = (1.0_WP, 0.0_WP)
    !! embed s
      SE = (0.0_WP, 0.0_WP)
      SE(1, 1) = s(1, 1)
      SE(1, 3) = s(1, 2)
      SE(3, 1) = s(2, 1)
      SE(3, 3) = s(2, 2)
      SE(2, 2) = (1.0_WP, 0.0_WP)
    !! embed t
      TE = (0.0_WP, 0.0_WP)
      TE(1, 1) = (1.0_WP, 0.0_WP)
      TE(2, 2) = t(1, 1)
      TE(2, 3) = t(1, 2)
      TE(3, 2) = t(2, 1)
      TE(3, 3) = t(2, 2)
      U = stdlib_matmul(RE, SE, TE)
   end function constructSU3Matrix

end module FLUE_SU3_random
