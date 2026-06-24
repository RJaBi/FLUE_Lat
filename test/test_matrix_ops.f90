module test_matrix_ops
   use FLUE_constants, only: WP, WC
   use FLUE_matrixConstants, only: Ident3x3
   use FLUE_SU3MatrixOps, only: MultiplyMatMat, MultiplyMatMatDag, &
                                TraceMultMatMat, RealTraceMultMatMat, &
                                TracelessConjgSubtract, RealTraceMat, &
                                FixSU3Matrix, ExpIQ
   use test_helpers, only: assert_close_real, assert_close_complex, &
                           assert_close_mat3, assert_is_unitary3, &
                           assert_det_one3, assert_traceless3
   use testdrive, only: new_unittest, unittest_type, error_type
   implicit none(type, external)
   private
   public :: collect_matrix_ops

contains

   subroutine collect_matrix_ops(testsuite)
      type(unittest_type), allocatable, intent(OUT) :: testsuite(:)
      testsuite = [ &
                  new_unittest("multiply_matmat_full", test_multiply_matmat_full), &
                  new_unittest("trace_mult_matmat", test_trace_mult_matmat), &
                  new_unittest("real_trace_mat", test_real_trace_mat), &
                  new_unittest("traceless_conjg_subtract", test_traceless_conjg_subtract), &
                  new_unittest("fix_su3_matrix_full", test_fix_su3_matrix_full), &
                  new_unittest("fix_su3_matrix_idempotent", test_fix_su3_matrix_idempotent), &
                  new_unittest("exp_iq_zero", test_exp_iq_zero), &
                  new_unittest("exp_iq_unitary_and_inverse", test_exp_iq_unitary_and_inverse) &
                  ]
   end subroutine collect_matrix_ops

   subroutine test_multiply_matmat_full(error)
     !! Test multiplyMatMat by comparing to intrinsic matmul
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC) :: a(3, 3), b(3, 3), got(3, 3), ref(3, 3)
      real(WP), parameter :: tol = 1.0E-14_WP

      a = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      b = CMPLX(0.0_WP, 0.0_WP, kind=WC)

      a(1, 1) = CMPLX(1.0_WP, 0.0_WP, kind=WC)
      a(1, 2) = CMPLX(2.0_WP, 1.0_WP, kind=WC)
      a(2, 2) = CMPLX(-1.0_WP, 0.5_WP, kind=WC)
      a(3, 1) = CMPLX(0.0_WP, -2.0_WP, kind=WC)
      a(3, 3) = CMPLX(4.0_WP, 0.0_WP, kind=WC)

      b(1, 1) = CMPLX(2.0_WP, 0.0_WP, kind=WC)
      b(1, 3) = CMPLX(-1.0_WP, 1.0_WP, kind=WC)
      b(2, 1) = CMPLX(3.0_WP, -1.0_WP, kind=WC)
      b(2, 2) = CMPLX(0.5_WP, 0.0_WP, kind=WC)
      b(3, 3) = CMPLX(1.0_WP, -2.0_WP, kind=WC)

      call MultiplyMatMat(got, a, b)
      ref = MATMUL(a, b)

      call assert_close_mat3(error, got, ref, tol, "MultiplyMatMat should match intrinsic MATMUL")
   end subroutine test_multiply_matmat_full

   subroutine test_trace_mult_matmat(error)
     !! Test the trace multiply subroutines
     !! By comparing to doing it explicitly
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC) :: a(3, 3), b(3, 3), prod(3, 3), tr_ref, tr_got
      real(WP) :: rtr_ref, rtr_got
      real(WP), parameter :: tol = 1.0E-14_WP

      a = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      b = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      a(1, 1) = CMPLX(1.0_WP, 0.5_WP, kind=WC)
      a(2, 2) = CMPLX(2.0_WP, -0.5_WP, kind=WC)
      a(3, 3) = CMPLX(-1.0_WP, 0.0_WP, kind=WC)
      b(1, 1) = CMPLX(0.5_WP, 1.0_WP, kind=WC)
      b(2, 2) = CMPLX(3.0_WP, 0.0_WP, kind=WC)
      b(3, 3) = CMPLX(2.0_WP, -1.0_WP, kind=WC)

      prod = MATMUL(a, b)
      tr_ref = prod(1, 1) + prod(2, 2) + prod(3, 3)
      rtr_ref = real(tr_ref, kind=WP)

      call TraceMultMatMat(tr_got, a, b)
      call RealTraceMultMatMat(rtr_got, a, b)

      call assert_close_complex(error, tr_got, tr_ref, tol, "TraceMultMatMat should match explicit trace")
      if (ALLOCATED(error)) return
      call assert_close_real(error, rtr_got, rtr_ref, tol, "RealTraceMultMatMat should match explicit real(trace)")
   end subroutine test_trace_mult_matmat

   subroutine test_real_trace_mat(error)
     !! Test taking the trace by knowing expected
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC) :: a(3, 3)
      real(WP), parameter :: tol = 1.0E-14_WP
      real(WP) :: expected

      a = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      a(1, 1) = CMPLX(1.0_WP, 0.5_WP, kind=WC)
      a(2, 2) = CMPLX(2.0_WP, -0.25_WP, kind=WC)
      a(3, 3) = CMPLX(-0.5_WP, 0.0_WP, kind=WC)
      expected = 2.5_WP

      call assert_close_real(error, RealTraceMat(a), expected, tol, "RealTraceMat should return real(trace)")
   end subroutine test_real_trace_mat

   subroutine test_traceless_conjg_subtract(error)
     !! Tests taking the trace and making traceless matrix
     !! By asserting that it returns a traceless matrix
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC) :: a(3, 3), b(3, 3), c(3, 3)
      real(WP), parameter :: tol = 1.0E-12_WP

      a = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      b = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      a(1, 2) = CMPLX(1.0_WP, 2.0_WP, kind=WC)
      a(2, 3) = CMPLX(-0.5_WP, 0.5_WP, kind=WC)
      b(2, 1) = CMPLX(3.0_WP, -1.0_WP, kind=WC)
      b(3, 2) = CMPLX(0.25_WP, 2.0_WP, kind=WC)

      call TracelessConjgSubtract(c, a, b)
      call assert_traceless3(error, c, tol, "TracelessConjgSubtract should return a traceless matrix")
   end subroutine test_traceless_conjg_subtract

   subroutine test_fix_su3_matrix_full(error)
     !! Tests that fixSU3 does project to SU3 properly
     !! By checking unitarity and determinant
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC) :: u(3, 3)
      real(WP), parameter :: tol = 1.0E-8_WP

      u = Ident3x3
      u(1, 1) = CMPLX(1.1_WP, 0.1_WP, kind=WC)
      u(2, 1) = CMPLX(0.02_WP, -0.01_WP, kind=WC)
      u(3, 2) = CMPLX(-0.03_WP, 0.04_WP, kind=WC)

      call FixSU3Matrix(u)

      call assert_is_unitary3(error, u, tol, "FixSU3Matrix should restore unitarity")
      if (ALLOCATED(error)) return
      call assert_det_one3(error, u, tol, "FixSU3Matrix should restore determinant one")
   end subroutine test_fix_su3_matrix_full

   subroutine test_fix_su3_matrix_idempotent(error)
     !! Check that fixSU3 doesn't change result if already in SU3
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC) :: u(3, 3), u1(3, 3), u2(3, 3)
      real(WP), parameter :: tol = 1.0E-10_WP

      u = Ident3x3
      u(1, 1) = CMPLX(0.95_WP, 0.02_WP, kind=WC)
      u(1, 3) = CMPLX(0.01_WP, 0.00_WP, kind=WC)
      u(2, 1) = CMPLX(-0.03_WP, 0.01_WP, kind=WC)

      u1 = u
      call FixSU3Matrix(u1)
      u2 = u1
      call FixSU3Matrix(u2)

      call assert_close_mat3(error, u2, u1, tol, "FixSU3Matrix should be idempotent up to roundoff")
   end subroutine test_fix_su3_matrix_idempotent

   subroutine test_exp_iq_zero(error)
     !! Testing that exp(i*Q) where Q = 0 gives identity
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC) :: q(3, 3), v(3, 3)
      real(WP), parameter :: tol = 1.0E-14_WP

      q = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      v = ExpIQ(q)

      call assert_close_mat3(error, v, Ident3x3, tol, "ExpIQ(0) should be identity")
   end subroutine test_exp_iq_zero

   subroutine test_exp_iq_unitary_and_inverse(error)
     !! Testing that expIQ works properly
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC), dimension(3, 3) :: q, v, vinv, ref, negQ
      real(WP), parameter :: tol = 1.0E-10_WP

      q = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      ! Hermitian and traceless
      q(1, 1) = CMPLX(-0.2_WP, 0.0_WP, kind=WC)
      q(2, 2) = CMPLX(0.1_WP, 0.0_WP, kind=WC)
      q(3, 3) = CMPLX(0.1_WP, 0.0_WP, kind=WC)
      q(1, 2) = CMPLX(0.03_WP, 0.04_WP, kind=WC)
      q(2, 1) = CONJG(q(1, 2))

      v = ExpIQ(q)
      negQ = -q
      vinv = ExpIQ(negQ)
      ref = CONJG(TRANSPOSE(v))

      call assert_is_unitary3(error, v, tol, "ExpIQ(Hermitian traceless Q) should be unitary")
      if (ALLOCATED(error)) return
      call assert_det_one3(error, v, tol, "ExpIQ(Hermitian traceless Q) should have determinant one")
      if (ALLOCATED(error)) return
      call assert_close_mat3(error, vinv, ref, tol, "ExpIQ(-Q) should equal ExpIQ(Q)^\u2020")
   end subroutine test_exp_iq_unitary_and_inverse

end module test_matrix_ops
