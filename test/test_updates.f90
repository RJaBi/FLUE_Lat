MODULE test_updates
  USE FLUE_constants,       ONLY : WP, WC
  USE FLUE_heatbath,        ONLY : updateLinks
  USE FLUE_SU2_heatbath,    ONLY : SU2_updateLinks, constructXMatrix
  USE test_helpers,         ONLY : assert_close_mat2, &
       assert_is_unitary2, assert_is_unitary3, &
       assert_det_one2, assert_det_one3, seed_rng_fixed, &
       fill_identity_su2, fill_identity_su3
  USE testdrive,            ONLY : new_unittest, unittest_type, error_type, check
  IMPLICIT NONE(TYPE, EXTERNAL)
  PRIVATE

   PUBLIC :: collect_updates

CONTAINS

   !=========================================================
   ! Collect tests
   !=========================================================
   SUBROUTINE collect_updates(testsuite)
      TYPE(unittest_type), ALLOCATABLE, INTENT(OUT) :: testsuite(:)

      testsuite = (/ &
         new_unittest("updateLinks_preserves_su3",      test_updateLinks_preserves_su3), &
         new_unittest("su2_updateLinks_preserves_su2",  test_su2_updateLinks_preserves_su2), &
         new_unittest("constructXMatrix_returns_su2",   test_constructXMatrix_returns_su2) /)
   END SUBROUTINE collect_updates


   !=========================================================
   ! Batch invariant checks
   !=========================================================
   SUBROUTINE assert_all_links_su3(error, U, tol, message)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC), INTENT(IN) :: U(:,:,:,:,:,:,:)
      REAL(WP),    INTENT(IN) :: tol
      CHARACTER(len=*), INTENT(IN) :: message
      INTEGER :: mu, nt, nx, ny, nz

      DO mu = 1, size(U,3)
         DO nt = 1, size(U,4)
            DO nx = 1, size(U,5)
               DO ny = 1, size(U,6)
                  DO nz = 1, size(U,7)
                     CALL assert_is_unitary3(error, U(:,:,mu,nt,nx,ny,nz), tol, trim(message)//": unitarity")
                     IF (allocated(error)) RETURN
                     CALL assert_det_one3(error, U(:,:,mu,nt,nx,ny,nz), tol, trim(message)//": determinant")
                     IF (allocated(error)) RETURN
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE assert_all_links_su3

   SUBROUTINE assert_all_links_su2(error, U, tol, message)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC), INTENT(IN) :: U(:,:,:,:,:,:,:)
      REAL(WP),    INTENT(IN) :: tol
      CHARACTER(len=*), INTENT(IN) :: message
      INTEGER :: mu, nt, nx, ny, nz

      DO mu = 1, size(U,3)
         DO nt = 1, size(U,4)
            DO nx = 1, size(U,5)
               DO ny = 1, size(U,6)
                  DO nz = 1, size(U,7)
                     CALL assert_is_unitary2(error, U(:,:,mu,nt,nx,ny,nz), tol, trim(message)//": unitarity")
                     IF (allocated(error)) RETURN
                     CALL assert_det_one2(error, U(:,:,mu,nt,nx,ny,nz), tol, trim(message)//": determinant")
                     IF (allocated(error)) RETURN
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE assert_all_links_su2


   !=========================================================
   ! 1) SU(3) update keeps all links in SU(3)
   !=========================================================
   SUBROUTINE test_updateLinks_preserves_su3(error)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC) :: U(3,3,4,2,2,2,2)
      COMPLEX(WC) :: UUpdated(3,3,4,2,2,2,2)
      REAL(WP), PARAMETER :: beta = 6.0_WP
      REAL(WP), PARAMETER :: xi   = 1.0_WP
      REAL(WP), PARAMETER :: tol  = 1.0e-10_WP

      CALL fill_identity_su3(U)
      UUpdated = U

      CALL seed_rng_fixed(777)
      CALL updateLinks(U, beta, xi, UUpdated)

      CALL assert_all_links_su3(error, UUpdated, tol, "updateLinks should keep all links in SU(3)")
   END SUBROUTINE test_updateLinks_preserves_su3


   !=========================================================
   ! 3) SU(2) update keeps all links in SU(2)
   !=========================================================
   SUBROUTINE test_su2_updateLinks_preserves_su2(error)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC) :: U(2,2,4,2,2,2,2)
      COMPLEX(WC) :: UUpdated(2,2,4,2,2,2,2)
      REAL(WP), PARAMETER :: beta = 2.3_WP
      REAL(WP), PARAMETER :: tol  = 1.0e-10_WP

      CALL fill_identity_su2(U)
      UUpdated = U

      CALL seed_rng_fixed(888)
      CALL SU2_updateLinks(U, beta, UUpdated)

      CALL assert_all_links_su2(error, UUpdated, tol, "SU2_updateLinks should keep all links in SU(2)")
   END SUBROUTINE test_su2_updateLinks_preserves_su2


   !=========================================================
   ! 4) constructXMatrix(alpha,beta) should return SU(2)
   !=========================================================
   SUBROUTINE test_constructXMatrix_returns_su2(error)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC) :: X(2,2)
      REAL(WP), PARAMETER :: tol = 1.0e-12_WP
      REAL(WP), DIMENSION(4) :: alpha_vals, beta_vals
      INTEGER :: i

      alpha_vals = [0.10_WP, 0.50_WP, 1.00_WP, 2.00_WP]
      beta_vals  = [1.50_WP, 2.30_WP, 4.00_WP, 6.00_WP]

      DO i = 1, size(alpha_vals)
         X = constructXMatrix(alpha_vals(i), beta_vals(i))
         CALL assert_is_unitary2(error, X, tol, "constructXMatrix should return a unitary SU(2) matrix")
         IF (allocated(error)) RETURN
         CALL assert_det_one2(error, X, tol, "constructXMatrix should return determinant one")
         IF (allocated(error)) RETURN
      END DO
   END SUBROUTINE test_constructXMatrix_returns_su2

END MODULE test_updates
