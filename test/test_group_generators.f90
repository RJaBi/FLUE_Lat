module test_group_generators
   use FLUE_constants, only: WP, WC
   use FLUE_SU2_random, only: constructSU2Matrix
   use FLUE_SU3_random, only: constructSU3Matrix
   use philox, only: philox_fill_uniform_cc, C64
   use stdlib_random, only: dist_rand
   use test_helpers, only: assert_is_unitary2, assert_det_one2, &
                           assert_is_unitary3, assert_det_one3, seed_rng_fixed
   use testdrive, only: new_unittest, unittest_type, error_type
   implicit none(type, external)
   private
   public :: collect_group_generators

contains

   subroutine collect_group_generators(testsuite)
      type(unittest_type), allocatable, intent(OUT) :: testsuite(:)
      testsuite = [ &
                  new_unittest("su2_generator_properties", test_su2_generator_properties), &
                  new_unittest("su3_generator_properties", test_su3_generator_properties) &
                  ]
   end subroutine collect_group_generators

   subroutine test_su2_generator_properties(error)
     !! Generate a bunch of SU2 matrices
     !! Check that they are indeed SU2 matrices
      type(error_type), allocatable, intent(OUT) :: error
      real(kind=WP) :: r(4)
      complex(kind=WC) :: U(2, 2)
      integer(kind=C64), dimension(4) :: counter
      integer(kind=C64), dimension(2) :: key
      integer(kind=C64) :: randInt
      integer :: i
      real(kind=WP), parameter :: tol = 1.0E-12_WP
      call seed_rng_fixed(12345)
      do i = 1, 100
         randInt = dist_rand(C64)
         counter(1) = randInt
         randInt = dist_rand(C64)
         counter(2) = randInt
         randInt = dist_rand(C64)
         counter(3) = randInt
         randInt = dist_rand(C64)
         counter(4) = randInt
         randInt = dist_rand(C64)
         key(1) = randInt
         randInt = dist_rand(C64)
         key(2) = randInt
         call philox_fill_uniform_cc(counter, key, 1_C64, 1.0_WP, 1.0_WP, r)
         !r = randomNumbers()
         U = constructSU2Matrix(r)
         call assert_is_unitary2(error, U, tol, "constructSU2Matrix should return a unitary matrix")
         if (ALLOCATED(error)) return
         call assert_det_one2(error, U, tol, "constructSU2Matrix should return determinant one")
         if (ALLOCATED(error)) return
      end do
   end subroutine test_su2_generator_properties

   subroutine test_su3_generator_properties(error)
     !! Construct SU3 matrices
     !! Check that they are indeed unitary with det 1
      type(error_type), allocatable, intent(OUT) :: error
      real(kind=WP) :: r1(4), r2(4), r3(4)
      complex(kind=WC) :: R(2, 2), S(2, 2), T(2, 2), U(3, 3)
      integer :: i
      real(kind=WP), parameter :: tol = 1.0E-12_WP
      integer(kind=C64), dimension(4) :: counter
      integer(kind=C64), dimension(2) :: key
      integer(kind=C64) :: randInt
      call seed_rng_fixed(12345)

      do i = 1, 100
         randInt = dist_rand(C64)
         counter(1) = randInt
         randInt = dist_rand(C64)
         counter(2) = randInt
         randInt = dist_rand(C64)
         counter(3) = randInt
         randInt = dist_rand(C64)
         counter(4) = randInt
         randInt = dist_rand(C64)
         key(1) = randInt
         randInt = dist_rand(C64)
         key(2) = randInt
         call philox_fill_uniform_cc(counter, key, 1_C64, 1.0_WP, 1.0_WP, r1)
         randInt = dist_rand(C64)
         counter(1) = randInt
         randInt = dist_rand(C64)
         counter(2) = randInt
         randInt = dist_rand(C64)
         counter(3) = randInt
         randInt = dist_rand(C64)
         counter(4) = randInt
         randInt = dist_rand(C64)
         key(1) = randInt
         randInt = dist_rand(C64)
         key(2) = randInt
         call philox_fill_uniform_cc(counter, key, 1_C64, 1.0_WP, 1.0_WP, r2)
         randInt = dist_rand(C64)
         counter(1) = randInt
         randInt = dist_rand(C64)
         counter(2) = randInt
         randInt = dist_rand(C64)
         counter(3) = randInt
         randInt = dist_rand(C64)
         counter(4) = randInt
         randInt = dist_rand(C64)
         key(1) = randInt
         randInt = dist_rand(C64)
         key(2) = randInt
         call philox_fill_uniform_cc(counter, key, 1_C64, 1.0_WP, 1.0_WP, r3)
         !r1 = randomNumbers()
         !r2 = randomNumbers()
         !r3 = randomNumbers()
         R = constructSU2Matrix(r1)
         S = constructSU2Matrix(r2)
         T = constructSU2Matrix(r3)
         U = constructSU3Matrix(R, S, T)
         call assert_is_unitary3(error, U, tol, "constructSU3Matrix should return a unitary matrix")
         if (ALLOCATED(error)) return
         call assert_det_one3(error, U, tol, "constructSU3Matrix should return determinant one")
         if (ALLOCATED(error)) return
      end do
   end subroutine test_su3_generator_properties

end module test_group_generators
