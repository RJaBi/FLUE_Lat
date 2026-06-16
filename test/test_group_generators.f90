MODULE test_group_generators
   USE FLUE_constants, ONLY : WP, WC
   USE FLUE_SU2_random, ONLY : constructSU2Matrix, randomNumbers
   USE FLUE_SU3_random, ONLY : constructSU3Matrix
   USE test_helpers, ONLY : assert_is_unitary2, assert_det_one2, &
                            assert_is_unitary3, assert_det_one3, seed_rng_fixed
   USE testdrive, ONLY : new_unittest, unittest_type, error_type
   IMPLICIT NONE(TYPE, EXTERNAL)
   PRIVATE
   PUBLIC :: collect_group_generators

CONTAINS

   SUBROUTINE collect_group_generators(testsuite)
      TYPE(unittest_type), ALLOCATABLE, INTENT(OUT) :: testsuite(:)
      testsuite = [ &
         new_unittest("su2_generator_properties", test_su2_generator_properties), &
         new_unittest("su3_generator_properties", test_su3_generator_properties)  &
      ]
   END SUBROUTINE collect_group_generators

   SUBROUTINE test_su2_generator_properties(error)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      REAL(WP) :: r(4)
      COMPLEX(WC) :: U(2,2)
      INTEGER :: i
      REAL(WP), PARAMETER :: tol = 1.0e-12_WP

      CALL seed_rng_fixed(12345)

      DO i = 1, 100
         r = randomNumbers()
         U = constructSU2Matrix(r)
         CALL assert_is_unitary2(error, U, tol, "constructSU2Matrix should return a unitary matrix")
         IF (allocated(error)) RETURN
         CALL assert_det_one2(error, U, tol, "constructSU2Matrix should return determinant one")
         IF (allocated(error)) RETURN
      END DO
   END SUBROUTINE test_su2_generator_properties

   SUBROUTINE test_su3_generator_properties(error)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      REAL(WP) :: r1(4), r2(4), r3(4)
      COMPLEX(WC) :: R(2,2), S(2,2), T(2,2), U(3,3)
      INTEGER :: i
      REAL(WP), PARAMETER :: tol = 1.0e-12_WP

      CALL seed_rng_fixed(12345)

      DO i = 1, 100
         r1 = randomNumbers()
         r2 = randomNumbers()
         r3 = randomNumbers()
         R = constructSU2Matrix(r1)
         S = constructSU2Matrix(r2)
         T = constructSU2Matrix(r3)
         U = constructSU3Matrix(R, S, T)
         CALL assert_is_unitary3(error, U, tol, "constructSU3Matrix should return a unitary matrix")
         IF (allocated(error)) RETURN
         CALL assert_det_one3(error, U, tol, "constructSU3Matrix should return determinant one")
         IF (allocated(error)) RETURN
      END DO
   END SUBROUTINE test_su3_generator_properties

END MODULE test_group_generators
