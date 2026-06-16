MODULE test_matrix_ops
   USE FLUE_constants, ONLY : WP, WC
   USE FLUE_matrixConstants, ONLY : Ident3x3
   USE FLUE_SU3MatrixOps, ONLY : MultiplyMatMat, MultiplyMatMatDag, &
                                 TraceMultMatMat, RealTraceMultMatMat, &
                                 TracelessConjgSubtract, RealTraceMat, &
                                 FixSU3Matrix, ExpIQ
   USE test_helpers, ONLY : assert_close_real, assert_close_complex, &
                            assert_close_mat3, assert_is_unitary3, &
                            assert_det_one3, assert_traceless3
   USE testdrive, ONLY : new_unittest, unittest_type, error_type
   IMPLICIT NONE(TYPE, EXTERNAL)
   PRIVATE
   PUBLIC :: collect_matrix_ops

CONTAINS

   SUBROUTINE collect_matrix_ops(testsuite)
      TYPE(unittest_type), ALLOCATABLE, INTENT(OUT) :: testsuite(:)
      testsuite = [ &
         new_unittest("multiply_matmat_full",         test_multiply_matmat_full), &
         new_unittest("trace_mult_matmat",            test_trace_mult_matmat), &
         new_unittest("real_trace_mat",               test_real_trace_mat), &
         new_unittest("traceless_conjg_subtract",     test_traceless_conjg_subtract), &
         new_unittest("fix_su3_matrix_full",          test_fix_su3_matrix_full), &
         new_unittest("fix_su3_matrix_idempotent",    test_fix_su3_matrix_idempotent), &
         new_unittest("exp_iq_zero",                  test_exp_iq_zero), &
         new_unittest("exp_iq_unitary_and_inverse",   test_exp_iq_unitary_and_inverse) &
      ]
   END SUBROUTINE collect_matrix_ops

   SUBROUTINE test_multiply_matmat_full(error)
     !! Test multiplyMatMat by comparing to intrinsic matmul
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC) :: a(3,3), b(3,3), got(3,3), ref(3,3)
      REAL(WP), PARAMETER :: tol = 1.0e-14_WP

      a = cmplx(0.0_WP, 0.0_WP, kind=WC)
      b = cmplx(0.0_WP, 0.0_WP, kind=WC)

      a(1,1) = cmplx(1.0_WP, 0.0_WP, kind=WC)
      a(1,2) = cmplx(2.0_WP, 1.0_WP, kind=WC)
      a(2,2) = cmplx(-1.0_WP, 0.5_WP, kind=WC)
      a(3,1) = cmplx(0.0_WP, -2.0_WP, kind=WC)
      a(3,3) = cmplx(4.0_WP, 0.0_WP, kind=WC)

      b(1,1) = cmplx(2.0_WP, 0.0_WP, kind=WC)
      b(1,3) = cmplx(-1.0_WP, 1.0_WP, kind=WC)
      b(2,1) = cmplx(3.0_WP, -1.0_WP, kind=WC)
      b(2,2) = cmplx(0.5_WP, 0.0_WP, kind=WC)
      b(3,3) = cmplx(1.0_WP, -2.0_WP, kind=WC)

      CALL MultiplyMatMat(got, a, b)
      ref = matmul(a, b)

      CALL assert_close_mat3(error, got, ref, tol, "MultiplyMatMat should match intrinsic MATMUL")
   END SUBROUTINE test_multiply_matmat_full

   SUBROUTINE test_trace_mult_matmat(error)
     !! Test the trace multiply subroutines
     !! By comparing to doing it explicitly
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC) :: a(3,3), b(3,3), prod(3,3), tr_ref, tr_got
      REAL(WP) :: rtr_ref, rtr_got
      REAL(WP), PARAMETER :: tol = 1.0e-14_WP

      a = cmplx(0.0_WP, 0.0_WP, kind=WC)
      b = cmplx(0.0_WP, 0.0_WP, kind=WC)
      a(1,1) = cmplx(1.0_WP,  0.5_WP, kind=WC)
      a(2,2) = cmplx(2.0_WP, -0.5_WP, kind=WC)
      a(3,3) = cmplx(-1.0_WP, 0.0_WP, kind=WC)
      b(1,1) = cmplx(0.5_WP, 1.0_WP, kind=WC)
      b(2,2) = cmplx(3.0_WP, 0.0_WP, kind=WC)
      b(3,3) = cmplx(2.0_WP, -1.0_WP, kind=WC)

      prod = matmul(a, b)
      tr_ref = prod(1,1) + prod(2,2) + prod(3,3)
      rtr_ref = real(tr_ref, kind=WP)

      CALL TraceMultMatMat(tr_got, a, b)
      CALL RealTraceMultMatMat(rtr_got, a, b)

      CALL assert_close_complex(error, tr_got, tr_ref, tol, "TraceMultMatMat should match explicit trace")
      IF (allocated(error)) RETURN
      CALL assert_close_real(error, rtr_got, rtr_ref, tol, "RealTraceMultMatMat should match explicit real(trace)")
   END SUBROUTINE test_trace_mult_matmat

   SUBROUTINE test_real_trace_mat(error)
     !! Test taking the trace by knowing expected
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC) :: a(3,3)
      REAL(WP), PARAMETER :: tol = 1.0e-14_WP
      REAL(WP) :: expected

      a = cmplx(0.0_WP, 0.0_WP, kind=WC)
      a(1,1) = cmplx(1.0_WP,  0.5_WP, kind=WC)
      a(2,2) = cmplx(2.0_WP, -0.25_WP, kind=WC)
      a(3,3) = cmplx(-0.5_WP, 0.0_WP, kind=WC)
      expected = 2.5_WP

      CALL assert_close_real(error, RealTraceMat(a), expected, tol, "RealTraceMat should return real(trace)")
   END SUBROUTINE test_real_trace_mat

   SUBROUTINE test_traceless_conjg_subtract(error)
     !! Tests taking the trace and making traceless matrix
     !! By asserting that it returns a traceless matrix
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC) :: a(3,3), b(3,3), c(3,3)
      REAL(WP), PARAMETER :: tol = 1.0e-12_WP

      a = cmplx(0.0_WP, 0.0_WP, kind=WC)
      b = cmplx(0.0_WP, 0.0_WP, kind=WC)
      a(1,2) = cmplx(1.0_WP,  2.0_WP, kind=WC)
      a(2,3) = cmplx(-0.5_WP, 0.5_WP, kind=WC)
      b(2,1) = cmplx(3.0_WP, -1.0_WP, kind=WC)
      b(3,2) = cmplx(0.25_WP, 2.0_WP, kind=WC)

      CALL TracelessConjgSubtract(c, a, b)
      CALL assert_traceless3(error, c, tol, "TracelessConjgSubtract should return a traceless matrix")
   END SUBROUTINE test_traceless_conjg_subtract

   SUBROUTINE test_fix_su3_matrix_full(error)
     !! Tests that fixSU3 does project to SU3 properly
     !! By checking unitarity and determinant
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC) :: u(3,3)
      REAL(WP), PARAMETER :: tol = 1.0e-8_WP

      u = Ident3x3
      u(1,1) = cmplx(1.1_WP, 0.1_WP, kind=WC)
      u(2,1) = cmplx(0.02_WP, -0.01_WP, kind=WC)
      u(3,2) = cmplx(-0.03_WP, 0.04_WP, kind=WC)

      CALL FixSU3Matrix(u)

      CALL assert_is_unitary3(error, u, tol, "FixSU3Matrix should restore unitarity")
      IF (allocated(error)) RETURN
      CALL assert_det_one3(error, u, tol, "FixSU3Matrix should restore determinant one")
   END SUBROUTINE test_fix_su3_matrix_full

   SUBROUTINE test_fix_su3_matrix_idempotent(error)
     !! Check that fixSU3 doesn't change result if already in SU3
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC) :: u(3,3), u1(3,3), u2(3,3)
      REAL(WP), PARAMETER :: tol = 1.0e-10_WP

      u = Ident3x3
      u(1,1) = cmplx(0.95_WP, 0.02_WP, kind=WC)
      u(1,3) = cmplx(0.01_WP, 0.00_WP, kind=WC)
      u(2,1) = cmplx(-0.03_WP, 0.01_WP, kind=WC)

      u1 = u
      CALL FixSU3Matrix(u1)
      u2 = u1
      CALL FixSU3Matrix(u2)

      CALL assert_close_mat3(error, u2, u1, tol, "FixSU3Matrix should be idempotent up to roundoff")
   END SUBROUTINE test_fix_su3_matrix_idempotent

   SUBROUTINE test_exp_iq_zero(error)
     !! Testing that exp(i*Q) where Q = 0 gives identity
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC) :: q(3,3), v(3,3)
      REAL(WP), PARAMETER :: tol = 1.0e-14_WP

      q = cmplx(0.0_WP, 0.0_WP, kind=WC)
      v = ExpIQ(q)

      CALL assert_close_mat3(error, v, Ident3x3, tol, "ExpIQ(0) should be identity")
   END SUBROUTINE test_exp_iq_zero

   SUBROUTINE test_exp_iq_unitary_and_inverse(error)
     !! Testing that expIQ works properly
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(WC) :: q(3,3), v(3,3), vinv(3,3), ref(3,3)
      REAL(WP), PARAMETER :: tol = 1.0e-10_WP

      q = cmplx(0.0_WP, 0.0_WP, kind=WC)
      ! Hermitian and traceless
      q(1,1) = cmplx(-0.2_WP, 0.0_WP, kind=WC)
      q(2,2) = cmplx( 0.1_WP, 0.0_WP, kind=WC)
      q(3,3) = cmplx( 0.1_WP, 0.0_WP, kind=WC)
      q(1,2) = cmplx( 0.03_WP, 0.04_WP, kind=WC)
      q(2,1) = conjg(q(1,2))

      v = ExpIQ(q)
      vinv = ExpIQ(-q)
      ref = conjg(transpose(v))

      CALL assert_is_unitary3(error, v, tol, "ExpIQ(Hermitian traceless Q) should be unitary")
      IF (allocated(error)) RETURN
      CALL assert_det_one3(error, v, tol, "ExpIQ(Hermitian traceless Q) should have determinant one")
      IF (allocated(error)) RETURN
      CALL assert_close_mat3(error, vinv, ref, tol, "ExpIQ(-Q) should equal ExpIQ(Q)^\u2020")
   END SUBROUTINE test_exp_iq_unitary_and_inverse

END MODULE test_matrix_ops
