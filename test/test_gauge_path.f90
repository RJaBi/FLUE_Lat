MODULE test_gauge_path
  USE FLUE_constants, ONLY : WP, WC
  USE FLUE_matrixConstants, ONLY : Ident3x3
  USE FLUE_wloops, ONLY : genericPath, periodCoord
  USE test_helpers, ONLY : assert_close_mat3, fill_identity_su3
  USE testdrive, ONLY : new_unittest, unittest_type, error_type, check
  IMPLICIT NONE(TYPE, EXTERNAL)
  PRIVATE
  PUBLIC :: collect_gauge_path

CONTAINS

  SUBROUTINE collect_gauge_path(testsuite)
    TYPE(unittest_type), ALLOCATABLE, INTENT(OUT) :: testsuite(:)
    testsuite = [ &
         new_unittest("periodcoord_no_wrap_interior",        test_periodcoord_no_wrap_interior), &
         new_unittest("periodcoord_wrap_forward_t",          test_periodcoord_wrap_forward_t), &
         new_unittest("periodcoord_wrap_backward_t",         test_periodcoord_wrap_backward_t), &
         new_unittest("periodcoord_wrap_all_axes",           test_periodcoord_wrap_all_axes), &
         new_unittest("genericpath_identity",                test_genericpath_identity), &
         new_unittest("genericpath_backtracking_su3",        test_genericpath_backtracking_su3), &
         new_unittest("genericpath_two_step_product",        test_genericpath_two_step_product), &
         new_unittest("genericpath_periodic_single_wrap",    test_genericpath_periodic_single_wrap), &
         new_unittest("genericpath_periodic_wrap_backtrack", test_genericpath_periodic_wrap_backtrack) &
         ]
  END SUBROUTINE collect_gauge_path


  PURE FUNCTION datashape_from_field(U) RESULT(ds)
    COMPLEX(WC), INTENT(IN) :: U(:,:,:,:,:,:,:)
    INTEGER :: ds(7)
    ds = shape(U)
  END FUNCTION datashape_from_field
  !=========================================================
  ! Direct periodCoord tests
  !=========================================================
  SUBROUTINE test_periodcoord_no_wrap_interior(error)
    !! Check that periodcoord does not effect the interior points
    TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
    COMPLEX(WC) :: data(3,3,4,2,2,2,2)
    INTEGER :: coord(4), got(4), expected(4), ds(7)
    CALL fill_identity_su3(data)
    ds = datashape_from_field(data)
    coord    = [1, 2, 1, 2]
    expected = [1, 2, 1, 2]
    got = periodCoord(coord, ds)
    CALL check(error, all(got == expected), &
         "periodCoord should leave interior coordinates unchanged")
  END SUBROUTINE test_periodcoord_no_wrap_interior

  SUBROUTINE test_periodcoord_wrap_forward_t(error)
    !! Check forward wrapping
    TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
    COMPLEX(WC) :: data(3,3,4,2,2,2,2)
    INTEGER :: coord(4), got(4), expected(4), ds(7)
    CALL fill_identity_su3(data)
    ds = datashape_from_field(data)
    ! NT = ds(4) = 2, so stepping to 3 should wrap to 1
    coord    = [3, 1, 1, 1]
    expected = [1, 1, 1, 1]
    got = periodCoord(coord, ds)
    CALL check(error, all(got == expected), &
         "periodCoord should wrap forward across the temporal boundary")
  END SUBROUTINE test_periodcoord_wrap_forward_t

  SUBROUTINE test_periodcoord_wrap_backward_t(error)
    !! Check backwards wrapping
    TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
    COMPLEX(WC) :: data(3,3,4,2,2,2,2)
    INTEGER :: coord(4), got(4), expected(4), ds(7)
    CALL fill_identity_su3(data)
    ds = datashape_from_field(data)
    ! 0 should wrap to NT = ds(4) = 2
    coord    = [0, 1, 1, 1]
    expected = [2, 1, 1, 1]
    got = periodCoord(coord, ds)
    CALL check(error, all(got == expected), &
         "periodCoord should wrap backward across the temporal boundary")
  END SUBROUTINE test_periodcoord_wrap_backward_t

  SUBROUTINE test_periodcoord_wrap_all_axes(error)
    ! Do a wrap on all axes
    TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
    COMPLEX(WC) :: data(3,3,4,2,2,2,2)
    INTEGER :: coord(4), got(4), expected(4), ds(7)
    CALL fill_identity_su3(data)
    ds = datashape_from_field(data)
    ! On a 2x2x2x2 lattice:
    ! nt=3 -> 1, nx=0 -> 2, ny=3 -> 1, nz=0 -> 2
    coord    = [3, 0, 3, 0]
    expected = [1, 2, 1, 2]
    got = periodCoord(coord, ds)
    CALL check(error, all(got == expected), &
         "periodCoord should wrap correctly on all four lattice axes")
  END SUBROUTINE test_periodcoord_wrap_all_axes
  !=========================================================
  ! Existing / strengthened genericPath tests
  !=========================================================
  SUBROUTINE test_genericpath_identity(error)
    !! Generic path on an identity field multiplies identity together
    !! I.e. should return identity
    TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
    COMPLEX(WC) :: data(3,3,4,2,2,2,2)
    INTEGER :: coord(4), path(4)
    COMPLEX(WC) :: P(3,3)
    REAL(WP), PARAMETER :: tol = 1.0e-14_WP
    CALL fill_identity_su3(data)
    coord = [1,1,1,1]
    path  = [1, -1, 1, -1]
    P = genericPath(data, coord, path)
    CALL assert_close_mat3(error, P, Ident3x3, tol, &
         "genericPath on identity field should return identity")
  END SUBROUTINE test_genericpath_identity

  SUBROUTINE test_genericpath_backtracking_su3(error)
    !! Make non-identity field
    !! Make sure that going forward then back cancels to give identity
    !! As gauge links are unitary
    TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
    COMPLEX(WC) :: data(3,3,4,2,2,2,2)
    COMPLEX(WC) :: A(3,3), P(3,3)
    INTEGER :: coord(4), path(2)
    INTEGER :: i,j,k,l
    REAL(WP), PARAMETER :: tol = 1.0e-12_WP
    CALL fill_identity_su3(data)
    A = cmplx(0.0_WP, 0.0_WP, kind=WC)
    A(1,1) = exp(cmplx(0.0_WP,  0.3_WP, kind=WC))
    A(2,2) = exp(cmplx(0.0_WP, -0.3_WP, kind=WC))
    A(3,3) = cmplx(1.0_WP,  0.0_WP, kind=WC)
    DO i = 1, size(data,4)
       DO j = 1, size(data,5)
          DO k = 1, size(data,6)
             DO l = 1, size(data,7)
                data(:,:,1,i,j,k,l) = A
             END DO
          END DO
       END DO
    END DO
    coord = [1,1,1,1]
    path  = [1, -1]
    P = genericPath(data, coord, path)
    CALL assert_close_mat3(error, P, Ident3x3, tol, &
         "genericPath(U_mu then -U_mu) should backtrack to identity")
  END SUBROUTINE test_genericpath_backtracking_su3

  SUBROUTINE test_genericpath_two_step_product(error)
    !! Check that going forward is working as intended
    TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
    COMPLEX(WC) :: data(3,3,4,2,2,2,2)
    COMPLEX(WC) :: A(3,3), P(3,3), ref(3,3)
    INTEGER :: coord(4), path(2)
    INTEGER :: i,j,k,l
    REAL(WP), PARAMETER :: tol = 1.0e-12_WP
    CALL fill_identity_su3(data)
    A = cmplx(0.0_WP, 0.0_WP, kind=WC)
    A(1,1) = exp(cmplx(0.0_WP,  0.2_WP, kind=WC))
    A(2,2) = exp(cmplx(0.0_WP, -0.2_WP, kind=WC))
    A(3,3) = cmplx(1.0_WP,  0.0_WP, kind=WC)
    DO i = 1, size(data,4)
       DO j = 1, size(data,5)
          DO k = 1, size(data,6)
             DO l = 1, size(data,7)
                data(:,:,1,i,j,k,l) = A
             END DO
          END DO
       END DO
    END DO
    coord = [1,1,1,1]
    path  = [1, 1]
    P = genericPath(data, coord, path)
    ref = matmul(A, A)
    CALL assert_close_mat3(error, P, ref, tol, &
         "genericPath on two identical forward links should equal A*A")
  END SUBROUTINE test_genericpath_two_step_product

  SUBROUTINE test_genericpath_periodic_single_wrap(error)
    !! Test generic path across boundary forward
    TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
    COMPLEX(WC) :: data(3,3,4,2,2,2,2)
    COMPLEX(WC) :: A(3,3), P(3,3)
    INTEGER :: coord(4), path(1)
    INTEGER :: i,j,k,l
    REAL(WP), PARAMETER :: tol = 1.0e-12_WP
    CALL fill_identity_su3(data)
    A = cmplx(0.0_WP, 0.0_WP, kind=WC)
    A(1,1) = exp(cmplx(0.0_WP,  0.17_WP, kind=WC))
    A(2,2) = exp(cmplx(0.0_WP, -0.17_WP, kind=WC))
    A(3,3) = cmplx(1.0_WP, 0.0_WP, kind=WC)
    DO i = 1, size(data,4)
       DO j = 1, size(data,5)
          DO k = 1, size(data,6)
             DO l = 1, size(data,7)
                data(:,:,1,i,j,k,l) = A
             END DO
          END DO
       END DO
    END DO
    ! Start at boundary in t-direction on a 2-site lattice and step forward once.
    coord = [2,1,1,1]
    path  = [1]
    P = genericPath(data, coord, path)
    CALL assert_close_mat3(error, P, A, tol, &
         "genericPath should wrap correctly across the periodic boundary for one forward step")
  END SUBROUTINE test_genericpath_periodic_single_wrap

  SUBROUTINE test_genericpath_periodic_wrap_backtrack(error)
    !! Test generic path across boundary backwards
    TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
    COMPLEX(WC) :: data(3,3,4,2,2,2,2)
    COMPLEX(WC) :: A(3,3), P(3,3)
    INTEGER :: coord(4), path(2)
    INTEGER :: i,j,k,l
    REAL(WP), PARAMETER :: tol = 1.0e-12_WP
    CALL fill_identity_su3(data)
    A = cmplx(0.0_WP, 0.0_WP, kind=WC)
    A(1,1) = exp(cmplx(0.0_WP,  0.29_WP, kind=WC))
    A(2,2) = exp(cmplx(0.0_WP, -0.29_WP, kind=WC))
    A(3,3) = cmplx(1.0_WP, 0.0_WP, kind=WC)
    DO i = 1, size(data,4)
       DO j = 1, size(data,5)
          DO k = 1, size(data,6)
             DO l = 1, size(data,7)
                data(:,:,1,i,j,k,l) = A
             END DO
          END DO
       END DO
    END DO
    coord = [2,1,1,1]
    path  = [1, -1]
    P = genericPath(data, coord, path)
    CALL assert_close_mat3(error, P, Ident3x3, tol, &
         "genericPath should backtrack correctly across a periodic boundary")
  END SUBROUTINE test_genericpath_periodic_wrap_backtrack
END MODULE test_gauge_path
