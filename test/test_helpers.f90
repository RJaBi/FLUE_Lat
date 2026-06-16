
MODULE test_helpers
   USE FLUE_constants, ONLY : WP, WC
   USE FLUE_matrixConstants, ONLY : Ident2x2, Ident3x3
   USE stdlib_random, ONLY: random_seed
   USE testdrive, ONLY : error_type, check
   IMPLICIT NONE(TYPE, EXTERNAL)
   PRIVATE

   ! Get's the max(abs(Matrix1-Matrix2))
   PUBLIC :: maxabs_mat2, maxabs_mat3
   PUBLIC :: det2x2, det3x3
   ! Various testdrive check's that two values
   ! are within some tolerance
   PUBLIC :: assert_close_real
   PUBLIC :: assert_close_complex
   ! These work using maxabs_mat
   PUBLIC :: assert_close_mat2
   PUBLIC :: assert_close_mat3
   ! These compare to identity matrices
   PUBLIC :: assert_is_identity2
   PUBLIC :: assert_is_identity3
   ! These construct U^dag * U
   ! and check it is identity
   PUBLIC :: assert_is_unitary2
   PUBLIC :: assert_is_unitary3
   ! These check the determinant is close to 1
   PUBLIC :: assert_det_one2
   PUBLIC :: assert_det_one3
   ! These check that the trace is 0 within tolerance
   PUBLIC :: assert_traceless3
   ! Sets the rng seed for stdlib random_seed
   PUBLIC :: seed_rng_fixed
   ! fills rank 7 complex arrays with the identity
   ! in the two leftmost indices
   ! zero elsewhere
   PUBLIC :: fill_identity_su2, fill_identity_su3

CONTAINS

   PURE REAL(WP) FUNCTION maxabs_mat2(a, b) RESULT(v)
      COMPLEX(WC), INTENT(IN) :: a(2,2), b(2,2)
      v = maxval(abs(a - b))
   END FUNCTION maxabs_mat2

   PURE REAL(WP) FUNCTION maxabs_mat3(a, b) RESULT(v)
      COMPLEX(WC), INTENT(IN) :: a(3,3), b(3,3)
      v = maxval(abs(a - b))
   END FUNCTION maxabs_mat3

   PURE COMPLEX(WC) FUNCTION det2x2(a) RESULT(d)
      COMPLEX(WC), INTENT(IN) :: a(2,2)
      d = a(1,1)*a(2,2) - a(1,2)*a(2,1)
   END FUNCTION det2x2

   PURE COMPLEX(WC) FUNCTION det3x3(a) RESULT(d)
      COMPLEX(WC), INTENT(IN) :: a(3,3)
      d = a(1,1)*(a(2,2)*a(3,3) - a(2,3)*a(3,2)) &
        - a(1,2)*(a(2,1)*a(3,3) - a(2,3)*a(3,1)) &
        + a(1,3)*(a(2,1)*a(3,2) - a(2,2)*a(3,1))
   END FUNCTION det3x3

   SUBROUTINE assert_close_real(error, actual, expected, tol, message)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      REAL(WP), INTENT(IN) :: actual, expected, tol
      CHARACTER(len=*), INTENT(IN) :: message
      CALL check(error, abs(actual - expected) < tol, message)
   END SUBROUTINE assert_close_real

   SUBROUTINE assert_close_complex(error, actual, expected, tol, message)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC), INTENT(IN) :: actual, expected
      REAL(WP), INTENT(IN) :: tol
      CHARACTER(len=*), INTENT(IN) :: message
      CALL check(error, abs(actual - expected) < tol, message)
   END SUBROUTINE assert_close_complex

   SUBROUTINE assert_close_mat2(error, actual, expected, tol, message)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC), INTENT(IN) :: actual(2,2), expected(2,2)
      REAL(WP), INTENT(IN) :: tol
      CHARACTER(len=*), INTENT(IN) :: message
      CALL check(error, maxabs_mat2(actual, expected) < tol, message)
   END SUBROUTINE assert_close_mat2

   SUBROUTINE assert_close_mat3(error, actual, expected, tol, message)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC), INTENT(IN) :: actual(3,3), expected(3,3)
      REAL(WP), INTENT(IN) :: tol
      CHARACTER(len=*), INTENT(IN) :: message
      CALL check(error, maxabs_mat3(actual, expected) < tol, message)
   END SUBROUTINE assert_close_mat3

   SUBROUTINE assert_is_identity2(error, a, tol, message)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC), INTENT(IN) :: a(2,2)
      REAL(WP), INTENT(IN) :: tol
      CHARACTER(len=*), INTENT(IN) :: message
      CALL assert_close_mat2(error, a, Ident2x2, tol, message)
   END SUBROUTINE assert_is_identity2

   SUBROUTINE assert_is_identity3(error, a, tol, message)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC), INTENT(IN) :: a(3,3)
      REAL(WP), INTENT(IN) :: tol
      CHARACTER(len=*), INTENT(IN) :: message
      CALL assert_close_mat3(error, a, Ident3x3, tol, message)
   END SUBROUTINE assert_is_identity3

   SUBROUTINE assert_is_unitary2(error, a, tol, message)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC), INTENT(IN) :: a(2,2)
      REAL(WP), INTENT(IN) :: tol
      CHARACTER(len=*), INTENT(IN) :: message
      COMPLEX(WC) :: udagu(2,2)
      udagu = matmul(conjg(transpose(a)), a)
      CALL assert_is_identity2(error, udagu, tol, message)
   END SUBROUTINE assert_is_unitary2

   SUBROUTINE assert_is_unitary3(error, a, tol, message)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC), INTENT(IN) :: a(3,3)
      REAL(WP), INTENT(IN) :: tol
      CHARACTER(len=*), INTENT(IN) :: message
      COMPLEX(WC) :: udagu(3,3)
      udagu = matmul(conjg(transpose(a)), a)
      CALL assert_is_identity3(error, udagu, tol, message)
   END SUBROUTINE assert_is_unitary3

   SUBROUTINE assert_det_one2(error, a, tol, message)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC), INTENT(IN) :: a(2,2)
      REAL(WP), INTENT(IN) :: tol
      CHARACTER(len=*), INTENT(IN) :: message
      CALL assert_close_complex(error, det2x2(a), cmplx(1.0_WP, 0.0_WP, kind=WC), tol, message)
   END SUBROUTINE assert_det_one2

   SUBROUTINE assert_det_one3(error, a, tol, message)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC), INTENT(IN) :: a(3,3)
      REAL(WP), INTENT(IN) :: tol
      CHARACTER(len=*), INTENT(IN) :: message
      CALL assert_close_complex(error, det3x3(a), cmplx(1.0_WP, 0.0_WP, kind=WC), tol, message)
   END SUBROUTINE assert_det_one3

   SUBROUTINE assert_traceless3(error, a, tol, message)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC), INTENT(IN) :: a(3,3)
      REAL(WP), INTENT(IN) :: tol
      CHARACTER(len=*), INTENT(IN) :: message
      COMPLEX(WC) :: tr
      tr = a(1,1) + a(2,2) + a(3,3)
      CALL check(error, abs(tr) < tol, message)
    END SUBROUTINE assert_traceless3


    SUBROUTINE seed_rng_fixed(seed_value)
      INTEGER, INTENT(IN) :: seed_value
      INTEGER :: seed_get
      CALL random_seed(put=seed_value, get=seed_get)
    END SUBROUTINE seed_rng_fixed

    !=========================================================
    ! Fillers
    !=========================================================
    SUBROUTINE fill_identity_su3(U)
      COMPLEX(WC), INTENT(OUT) :: U(:,:,:,:,:,:,:)
      INTEGER :: mu, nt, nx, ny, nz
      U = cmplx(0.0_WP, 0.0_WP, kind=WC)
      DO CONCURRENT (mu = 1:size(U,3), nt = 1:size(U,4), nx = 1:size(U,5), &
           ny = 1:size(U,6), nz = 1:size(U,7))
         U(:,:,mu,nt,nx,ny,nz) = Ident3x3
      END DO
    END SUBROUTINE fill_identity_su3

    SUBROUTINE fill_identity_su2(U)
      COMPLEX(WC), INTENT(OUT) :: U(:,:,:,:,:,:,:)
      INTEGER :: mu, nt, nx, ny, nz
      U = cmplx(0.0_WP, 0.0_WP, kind=WC)
      DO CONCURRENT (mu = 1:size(U,3), nt = 1:size(U,4), nx = 1:size(U,5), &
           ny = 1:size(U,6), nz = 1:size(U,7))
         U(:,:,mu,nt,nx,ny,nz) = Ident2x2
      END DO
    END SUBROUTINE fill_identity_su2

END MODULE test_helpers
