MODULE FLUE_SU3_random
   USE FLUE_constants, ONLY: WP, WC
   USE FLUE_SU2_random, ONLY: constructSU2Matrix
   USE stdlib_intrinsics, ONLY: stdlib_matmul
   IMPLICIT NONE(TYPE, EXTERNAL)
   PRIVATE

   PUBLIC :: constructSU3Matrix

CONTAINS

   PURE FUNCTION constructSU3Matrix(R, S, T) RESULT(U)
      COMPLEX(kind=WC), DIMENSION(2, 2), INTENT(IN) :: R, S, T
      COMPLEX(kind=WC), DIMENSION(3, 3) :: U
    !! embed r,s,t into 3x3 matrices
      COMPLEX(kind=WC), DIMENSION(3, 3) :: RE, SE, TE
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
   END FUNCTION constructSU3Matrix

END MODULE FLUE_SU3_random
