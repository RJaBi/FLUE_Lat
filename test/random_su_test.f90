MODULE test_random_su
   USE FLUE_constants, ONLY : WP, WC
   USE FLUE_matrixConstants, ONLY : Ident2x2, Ident3x3
   USE FLUE_SU3_random, ONLY : constructSU3Matrix
   USE test_helpers, ONLY : assert_close_mat3, assert_is_unitary3, assert_det_one3
   USE testdrive, ONLY : new_unittest, unittest_type, error_type
   IMPLICIT NONE(TYPE, EXTERNAL)
   PRIVATE
   PUBLIC :: collect_random_su

CONTAINS

   SUBROUTINE collect_random_su(testsuite)
      TYPE(unittest_type), ALLOCATABLE, INTENT(OUT) :: testsuite(:)
      testsuite = [ &
         new_unittest("construct_su3_identity", test_construct_su3_identity) &
      ]
   END SUBROUTINE collect_random_su

   SUBROUTINE test_construct_su3_identity(error)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC) :: R(2,2), S(2,2), T(2,2), U(3,3)
      REAL(WP), PARAMETER :: tol = 1.0e-14_WP

      R = Ident2x2
      S = Ident2x2
      T = Ident2x2
      U = constructSU3Matrix(R, S, T)

      CALL assert_close_mat3(error, U, Ident3x3, tol, "constructSU3Matrix(I,I,I) should equal identity")
      IF (allocated(error)) RETURN
      CALL assert_is_unitary3(error, U, tol, "constructSU3Matrix(I,I,I) should be unitary")
      IF (allocated(error)) RETURN
      CALL assert_det_one3(error, U, tol, "constructSU3Matrix(I,I,I) should have determinant one")
   END SUBROUTINE test_construct_su3_identity

END MODULE test_random_su
