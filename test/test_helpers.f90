
module test_helpers
   use FLUE_constants, only: WP, WC
   use FLUE_matrixConstants, only: Ident2x2, Ident3x3
   use stdlib_random, only: random_seed
   use testdrive, only: error_type, check
   implicit none(type, external)
   private

   ! Get's the max(abs(Matrix1-Matrix2))
   public :: maxabs_mat2, maxabs_mat3
   public :: det2x2, det3x3
   ! Various testdrive check's that two values
   ! are within some tolerance
   public :: assert_close_real
   public :: assert_close_complex
   ! These work using maxabs_mat
   public :: assert_close_mat2
   public :: assert_close_mat3
   ! These compare to identity matrices
   public :: assert_is_identity2
   public :: assert_is_identity3
   ! These construct U^dag * U
   ! and check it is identity
   public :: assert_is_unitary2
   public :: assert_is_unitary3
   ! These check the determinant is close to 1
   public :: assert_det_one2
   public :: assert_det_one3
   ! These check that the trace is 0 within tolerance
   public :: assert_traceless3
   ! Sets the rng seed for stdlib random_seed
   public :: seed_rng_fixed
   ! fills rank 7 complex arrays with the identity
   ! in the two leftmost indices
   ! zero elsewhere
   public :: fill_identity_su2, fill_identity_su3

contains

   pure real(WP) function maxabs_mat2(a, b) result(v)
      complex(WC), intent(IN) :: a(2, 2), b(2, 2)
      v = MAXVAL(ABS(a - b))
   end function maxabs_mat2

   pure real(WP) function maxabs_mat3(a, b) result(v)
      complex(WC), intent(IN) :: a(3, 3), b(3, 3)
      v = MAXVAL(ABS(a - b))
   end function maxabs_mat3

   pure complex(WC) function det2x2(a) result(d)
      complex(WC), intent(IN) :: a(2, 2)
      d = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)
   end function det2x2

   pure complex(WC) function det3x3(a) result(d)
      complex(WC), intent(IN) :: a(3, 3)
      d = a(1, 1) * (a(2, 2) * a(3, 3) - a(2, 3) * a(3, 2)) &
          - a(1, 2) * (a(2, 1) * a(3, 3) - a(2, 3) * a(3, 1)) &
          + a(1, 3) * (a(2, 1) * a(3, 2) - a(2, 2) * a(3, 1))
   end function det3x3

   subroutine assert_close_real(error, actual, expected, tol, message)
      type(error_type), allocatable, intent(OUT) :: error
      real(WP), intent(IN) :: actual, expected, tol
      character(len=*), intent(IN) :: message
      call check(error, ABS(actual - expected) < tol, message)
   end subroutine assert_close_real

   subroutine assert_close_complex(error, actual, expected, tol, message)
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC), intent(IN) :: actual, expected
      real(WP), intent(IN) :: tol
      character(len=*), intent(IN) :: message
      call check(error, ABS(actual - expected) < tol, message)
   end subroutine assert_close_complex

   subroutine assert_close_mat2(error, actual, expected, tol, message)
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC), intent(IN) :: actual(2, 2), expected(2, 2)
      real(WP), intent(IN) :: tol
      character(len=*), intent(IN) :: message
      call check(error, maxabs_mat2(actual, expected) < tol, message)
   end subroutine assert_close_mat2

   subroutine assert_close_mat3(error, actual, expected, tol, message)
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC), intent(IN) :: actual(3, 3), expected(3, 3)
      real(WP), intent(IN) :: tol
      character(len=*), intent(IN) :: message
      call check(error, maxabs_mat3(actual, expected) < tol, message)
   end subroutine assert_close_mat3

   subroutine assert_is_identity2(error, a, tol, message)
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC), intent(IN) :: a(2, 2)
      real(WP), intent(IN) :: tol
      character(len=*), intent(IN) :: message
      call assert_close_mat2(error, a, Ident2x2, tol, message)
   end subroutine assert_is_identity2

   subroutine assert_is_identity3(error, a, tol, message)
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC), intent(IN) :: a(3, 3)
      real(WP), intent(IN) :: tol
      character(len=*), intent(IN) :: message
      call assert_close_mat3(error, a, Ident3x3, tol, message)
   end subroutine assert_is_identity3

   subroutine assert_is_unitary2(error, a, tol, message)
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC), intent(IN) :: a(2, 2)
      real(WP), intent(IN) :: tol
      character(len=*), intent(IN) :: message
      complex(WC) :: udagu(2, 2)
      udagu = MATMUL(CONJG(TRANSPOSE(a)), a)
      call assert_is_identity2(error, udagu, tol, message)
   end subroutine assert_is_unitary2

   subroutine assert_is_unitary3(error, a, tol, message)
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC), intent(IN) :: a(3, 3)
      real(WP), intent(IN) :: tol
      character(len=*), intent(IN) :: message
      complex(WC) :: udagu(3, 3)
      udagu = MATMUL(CONJG(TRANSPOSE(a)), a)
      call assert_is_identity3(error, udagu, tol, message)
   end subroutine assert_is_unitary3

   subroutine assert_det_one2(error, a, tol, message)
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC), intent(IN) :: a(2, 2)
      real(WP), intent(IN) :: tol
      character(len=*), intent(IN) :: message
      call assert_close_complex(error, det2x2(a), CMPLX(1.0_WP, 0.0_WP, kind=WC), tol, message)
   end subroutine assert_det_one2

   subroutine assert_det_one3(error, a, tol, message)
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC), intent(IN) :: a(3, 3)
      real(WP), intent(IN) :: tol
      character(len=*), intent(IN) :: message
      call assert_close_complex(error, det3x3(a), CMPLX(1.0_WP, 0.0_WP, kind=WC), tol, message)
   end subroutine assert_det_one3

   subroutine assert_traceless3(error, a, tol, message)
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC), intent(IN) :: a(3, 3)
      real(WP), intent(IN) :: tol
      character(len=*), intent(IN) :: message
      complex(WC) :: tr
      tr = a(1, 1) + a(2, 2) + a(3, 3)
      call check(error, ABS(tr) < tol, message)
   end subroutine assert_traceless3

   subroutine seed_rng_fixed(seed_value)
      integer, intent(IN) :: seed_value
      integer :: seed_get
      call RANDOM_SEED(put=seed_value, get=seed_get)
   end subroutine seed_rng_fixed

   !=========================================================
   ! Fillers
   !=========================================================
   subroutine fill_identity_su3(U)
      complex(WC), intent(OUT) :: U(:, :, :, :, :, :, :)
      integer :: mu, nt, nx, ny, nz
      U = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      do concurrent(mu=1:SIZE(U, 3), nt=1:SIZE(U, 4), nx=1:SIZE(U, 5), &
                    ny=1:SIZE(U, 6), nz=1:SIZE(U, 7))
         U(:, :, mu, nt, nx, ny, nz) = Ident3x3
      end do
   end subroutine fill_identity_su3

   subroutine fill_identity_su2(U)
      complex(WC), intent(OUT) :: U(:, :, :, :, :, :, :)
      integer :: mu, nt, nx, ny, nz
      U = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      do concurrent(mu=1:SIZE(U, 3), nt=1:SIZE(U, 4), nx=1:SIZE(U, 5), &
                    ny=1:SIZE(U, 6), nz=1:SIZE(U, 7))
         U(:, :, mu, nt, nx, ny, nz) = Ident2x2
      end do
   end subroutine fill_identity_su2

end module test_helpers
