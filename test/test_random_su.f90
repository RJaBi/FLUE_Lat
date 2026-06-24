module test_random_su
   use FLUE_constants, only: WP, WC
   use FLUE_matrixConstants, only: Ident2x2, Ident3x3
   use FLUE_SU3_random, only: constructSU3Matrix
   use test_helpers, only: assert_close_mat3, assert_is_unitary3, assert_det_one3
   use testdrive, only: new_unittest, unittest_type, error_type
   implicit none(type, external)
   private
   public :: collect_random_su

contains

   subroutine collect_random_su(testsuite)
      type(unittest_type), allocatable, intent(OUT) :: testsuite(:)
      testsuite = [ &
                  new_unittest("construct_su3_identity", test_construct_su3_identity) &
                  ]
   end subroutine collect_random_su

   subroutine test_construct_su3_identity(error)
     !! This just tests that an SU3 matrix of 3 SU2 should be identity
     !! Which by definition is unitary and has determinant 1
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC) :: R(2, 2), S(2, 2), T(2, 2), U(3, 3)
      real(WP), parameter :: tol = 1.0E-14_WP
      R = Ident2x2
      S = Ident2x2
      T = Ident2x2
      U = constructSU3Matrix(R, S, T)
      call assert_close_mat3(error, U, Ident3x3, tol, "constructSU3Matrix(I,I,I) should equal identity")
      if (ALLOCATED(error)) return
      call assert_is_unitary3(error, U, tol, "constructSU3Matrix(I,I,I) should be unitary")
      if (ALLOCATED(error)) return
      call assert_det_one3(error, U, tol, "constructSU3Matrix(I,I,I) should have determinant one")
   end subroutine test_construct_su3_identity

end module test_random_su
