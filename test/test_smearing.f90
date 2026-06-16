MODULE test_smearing
  USE FLUE_constants,       ONLY : WP, WC
  USE FLUE_matrixConstants, ONLY : Ident3x3
  USE FLUE_stoutSmearing,   ONLY : StoutSmearLinks
  USE test_helpers,         ONLY : assert_close_mat3, &
                                    assert_is_unitary3, assert_det_one3, fill_identity_su3
  USE testdrive,            ONLY : new_unittest, unittest_type, error_type, check

  IMPLICIT NONE(TYPE, EXTERNAL)
  PRIVATE
  PUBLIC :: collect_smearing

CONTAINS

  !=========================================================
  ! Collect tests
  !=========================================================
  SUBROUTINE collect_smearing(testsuite)
    TYPE(unittest_type), ALLOCATABLE, INTENT(OUT) :: testsuite(:)
    testsuite = [ &
         new_unittest("stout_identity_preserved",      test_stout_identity_preserved), &
         new_unittest("stout_preserves_su3",           test_stout_preserves_su3), &
         new_unittest("stout_zero_sweeps_is_identity", test_stout_zero_sweeps_is_identity) &
         ]
   END SUBROUTINE collect_smearing
   !=========================================================
   ! Fillers
   !=========================================================
   SUBROUTINE fill_commuting_su3(U, theta)
     ! A simple exact SU(3) background:
     ! A = diag(exp(i theta), exp(-i theta), 1)
     COMPLEX(WC), INTENT(OUT) :: U(:,:,:,:,:,:,:)
     REAL(WP),    INTENT(IN)  :: theta
     COMPLEX(WC) :: A(3,3)
     INTEGER :: mu, nt, nx, ny, nz

     A = cmplx(0.0_WP, 0.0_WP, kind=WC)
     A(1,1) = exp(cmplx(0.0_WP,  theta, kind=WC))
     A(2,2) = exp(cmplx(0.0_WP, -theta, kind=WC))
     A(3,3) = cmplx(1.0_WP, 0.0_WP, kind=WC)

     DO CONCURRENT (mu = 1:size(U,3), nt = 1:size(U,4), nx = 1:size(U,5), &
          ny = 1:size(U,6), nz = 1:size(U,7))
        U(:,:,mu,nt,nx,ny,nz) = A
     END DO
   END SUBROUTINE fill_commuting_su3


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


    !=========================================================
    ! 1) Identity field should stay identity under stout smearing
    !=========================================================
    SUBROUTINE test_stout_identity_preserved(error)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC) :: U(3,3,4,2,2,2,2)
      COMPLEX(WC) :: USmear(3,3,4,2,2,2,2)
      INTEGER :: mu, nt, nx, ny, nz
      REAL(WP), PARAMETER :: rho = 0.10_WP
      INTEGER,  PARAMETER :: nSweeps = 1
      REAL(WP), PARAMETER :: tol = 1.0e-13_WP

      CALL fill_identity_su3(U)
      CALL StoutSmearLinks(U, rho, nSweeps, USmear)

      DO mu = 1, size(USmear,3)
         DO nt = 1, size(USmear,4)
            DO nx = 1, size(USmear,5)
               DO ny = 1, size(USmear,6)
                  DO nz = 1, size(USmear,7)
                     CALL assert_close_mat3(error, USmear(:,:,mu,nt,nx,ny,nz), Ident3x3, tol, &
                          "StoutSmearLinks should preserve the identity field")
                     IF (allocated(error)) RETURN
                  END DO
               END DO
            END DO
         END DO
      END DO
    END SUBROUTINE test_stout_identity_preserved


    !=========================================================
    ! 2) Smearing output should remain in SU(3)
    !=========================================================
    SUBROUTINE test_stout_preserves_su3(error)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC) :: U(3,3,4,2,2,2,2)
      COMPLEX(WC) :: USmear(3,3,4,2,2,2,2)
      REAL(WP), PARAMETER :: rho = 0.10_WP
      INTEGER,  PARAMETER :: nSweeps = 2
      REAL(WP), PARAMETER :: tol = 1.0e-10_WP

      CALL fill_commuting_su3(U, 0.23_WP)
      CALL StoutSmearLinks(U, rho, nSweeps, USmear)

      CALL assert_all_links_su3(error, USmear, tol, "StoutSmearLinks should preserve SU(3)")
    END SUBROUTINE test_stout_preserves_su3


    !=========================================================
   ! 3) Zero sweeps should return the input unchanged
    !=========================================================
    SUBROUTINE test_stout_zero_sweeps_is_identity(error)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC) :: U(3,3,4,2,2,2,2)
      COMPLEX(WC) :: USmear(3,3,4,2,2,2,2)
      REAL(WP), PARAMETER :: rho = 0.10_WP
      INTEGER,  PARAMETER :: nSweeps = 0
      REAL(WP), PARAMETER :: tol = 1.0e-14_WP

      CALL fill_commuting_su3(U, 0.31_WP)
      CALL StoutSmearLinks(U, rho, nSweeps, USmear)

      CALL assert_close_field_su3(error, USmear, U, tol, &
           "StoutSmearLinks with nSweeps=0 should return the input unchanged")
    END SUBROUTINE test_stout_zero_sweeps_is_identity


    SUBROUTINE assert_close_field_su3(error, A, B, tol, message)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC), INTENT(IN) :: A(:,:,:,:,:,:,:), B(:,:,:,:,:,:,:)
      REAL(WP),    INTENT(IN) :: tol
      CHARACTER(len=*), INTENT(IN) :: message

      CALL testdrive_check_close_field(error, maxval(abs(A - B)) < tol, message)
    END SUBROUTINE assert_close_field_su3


   SUBROUTINE testdrive_check_close_field(error, ok, message)
     TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
     LOGICAL, INTENT(IN) :: ok
     CHARACTER(len=*), INTENT(IN) :: message
     CALL check(error, ok, message)
   END SUBROUTINE testdrive_check_close_field

 END MODULE test_smearing
