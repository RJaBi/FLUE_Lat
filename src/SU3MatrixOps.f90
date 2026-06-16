! Some subroutines that act upon and create SU(3) (3x3 complex) matrices
MODULE FLUE_SU3MatrixOps
   USE FLUE_constants, ONLY: WC, WP
   USE FLUE_matrixConstants, ONLY: Ident3x3
   IMPLICIT NONE(TYPE, EXTERNAL)
   PRIVATE
   PUBLIC :: MultiplyMatMat, MultiplyMatDagMatDag, &
             TraceMultMatMat, RealTraceMultMatMat, TracelessConjgSubtract, &
             colourDecomp, RealTraceMat, MultiplyMatMatDag, TraceMat
   PUBLIC :: FixSU3Matrix
   PUBLIC :: orthogonalise_vectors, vector_product
   PUBLIC :: ExpIQ

CONTAINS

   ! stripped from cola and de-colour vectored
   PURE SUBROUTINE orthogonalise_vectors(w, v)
      COMPLEX(kind=WC), DIMENSION(3), INTENT(INOUT) :: w
      COMPLEX(kind=WC), DIMENSION(3), INTENT(IN) :: v
      COMPLEX(kind=WC) :: vdotw
      vdotw = SUM(CONJG(v) * w)
      w = w - v * vdotw
   END SUBROUTINE orthogonalise_vectors
   ! stripped from cola and de-colour vectored
   PURE SUBROUTINE vector_product(x, v, w)
      COMPLEX(kind=WC), DIMENSION(3), INTENT(OUT) :: x
      COMPLEX(kind=WC), DIMENSION(3), INTENT(IN) :: v, w
      INTEGER :: ic, jc, kc
      INTEGER, PARAMETER :: nc = 3
      DO ic = 1, nc
         jc = MODULO(ic, 3) + 1
         kc = MODULO(jc, 3) + 1
         x(ic) = CONJG(v(jc) * w(kc) - v(kc) * w(jc))
      END DO
   END SUBROUTINE vector_product
   ! stripped from cola and de-colour vectored
   PURE SUBROUTINE FixSU3Matrix(U_x)
      COMPLEX(kind=WC), DIMENSION(3, 3), INTENT(INOUT) :: U_x
      COMPLEX(kind=WC), DIMENSION(3) :: v1, v2, v3
      v1 = U_x(1, :)
      CALL normalise_vector(v1)
      v2(:) = U_x(2, :)
      CALL orthogonalise_vectors(v2, v1)
      CALL normalise_vector(v2)
      CALL vector_product(v3, v1, v2)
      CALL normalise_vector(v3)
      U_x(1, :) = v1(:)
      U_x(2, :) = v2(:)
      U_x(3, :) = v3(:)
   END SUBROUTINE FixSU3Matrix
   ! stripped from cola and de-colour vectored
   PURE SUBROUTINE normalise_vector(v)
      COMPLEX(kind=WC), DIMENSION(3), INTENT(INOUT) :: v
      REAL(WP) :: norm
      norm = SQRT(SUM(real(v)**2 + AIMAG(v)**2))
      v = v / norm
   END SUBROUTINE normalise_vector
   ! explicitly unroll 3x3 matrix mult
   PURE SUBROUTINE MultiplyMatMat(MM, left, right)
      COMPLEX(kind=WC), DIMENSION(3, 3), INTENT(IN) :: left, right
      COMPLEX(kind=WC), DIMENSION(3, 3), INTENT(OUT) :: MM
      !"""
      !Multiple left by right. Assumes 3x3 (colour) (complex) matrices
      !"""
      !# do the maths for the colour matrices
      MM(1, 1) = left(1, 1) * right(1, 1) + left(1, 2) * right(2, 1) + left(1, 3) * right(3, 1)
      MM(2, 1) = left(2, 1) * right(1, 1) + left(2, 2) * right(2, 1) + left(2, 3) * right(3, 1)
      MM(3, 1) = left(3, 1) * right(1, 1) + left(3, 2) * right(2, 1) + left(3, 3) * right(3, 1)
      !# second index
      MM(1, 2) = left(1, 1) * right(1, 2) + left(1, 2) * right(2, 2) + left(1, 3) * right(3, 2)
      MM(2, 2) = left(2, 1) * right(1, 2) + left(2, 2) * right(2, 2) + left(2, 3) * right(3, 2)
      MM(3, 2) = left(3, 1) * right(1, 2) + left(3, 2) * right(2, 2) + left(3, 3) * right(3, 2)
      !# third index
      MM(1, 3) = left(1, 1) * right(1, 3) + left(1, 2) * right(2, 3) + left(1, 3) * right(3, 3)
      MM(2, 3) = left(2, 1) * right(1, 3) + left(2, 2) * right(2, 3) + left(2, 3) * right(3, 3)
      MM(3, 3) = left(3, 1) * right(1, 3) + left(3, 2) * right(2, 3) + left(3, 3) * right(3, 3)
   END SUBROUTINE MultiplyMatMat

   PURE SUBROUTINE MultiplyMatMatDag(MM, left, right)
      ! A B^\dag
      COMPLEX(kind=WC), DIMENSION(3, 3), INTENT(IN) :: left, right
      COMPLEX(kind=WC), DIMENSION(3, 3), INTENT(OUT) :: MM
      !"""
      !Multiple left by right^\dagger. Assumes 3x3 (colour) (complex) matrices
      !"""
      !# do the maths for the colour matrices
      MM(1, 1) = left(1, 1) * CONJG(right(1, 1)) + left(1, 2) * CONJG(right(1, 2)) + left(1, 3) * CONJG(right(1, 3))
      MM(2, 1) = left(2, 1) * CONJG(right(1, 1)) + left(2, 2) * CONJG(right(1, 2)) + left(2, 3) * CONJG(right(1, 3))
      MM(3, 1) = left(3, 1) * CONJG(right(1, 1)) + left(3, 2) * CONJG(right(1, 2)) + left(3, 3) * CONJG(right(1, 3))
      !# second index
      MM(1, 2) = left(1, 1) * CONJG(right(2, 1)) + left(1, 2) * CONJG(right(2, 2)) + left(1, 3) * CONJG(right(2, 3))
      MM(2, 2) = left(2, 1) * CONJG(right(2, 1)) + left(2, 2) * CONJG(right(2, 2)) + left(2, 3) * CONJG(right(2, 3))
      MM(3, 2) = left(3, 1) * CONJG(right(2, 1)) + left(3, 2) * CONJG(right(2, 2)) + left(3, 3) * CONJG(right(2, 3))
      !# third index
      MM(1, 3) = left(1, 1) * CONJG(right(3, 1)) + left(1, 2) * CONJG(right(3, 2)) + left(1, 3) * CONJG(right(3, 3))
      MM(2, 3) = left(2, 1) * CONJG(right(3, 1)) + left(2, 2) * CONJG(right(3, 2)) + left(2, 3) * CONJG(right(3, 3))
      MM(3, 3) = left(3, 1) * CONJG(right(3, 1)) + left(3, 2) * CONJG(right(3, 2)) + left(3, 3) * CONJG(right(3, 3))
   END SUBROUTINE MultiplyMatMatDag

   PURE SUBROUTINE MultiplyMatdagMatdag(MM, left, right)
      COMPLEX(kind=WC), DIMENSION(3, 3), INTENT(IN) :: left, right
      COMPLEX(kind=WC), DIMENSION(3, 3), INTENT(OUT) :: MM
      !"""
      !#Multiplies two (3,3) complex matrices together. Takes conjugate
      !Does (left*right)^dagger
      !"""
      !# take transpose manually
      MM(1, 1) = CONJG(left(1, 1) * right(1, 1) + left(2, 1) * right(1, 2) + left(3, 1) * right(1, 3))
      MM(2, 1) = CONJG(left(1, 2) * right(1, 1) + left(2, 2) * right(1, 2) + left(3, 2) * right(1, 3))
      MM(3, 1) = CONJG(left(1, 3) * right(1, 1) + left(2, 3) * right(1, 2) + left(3, 3) * right(1, 3))
      !# but take conjugate using np
      MM(1, 2) = CONJG(left(1, 1) * right(2, 1) + left(2, 1) * right(2, 2) + left(3, 1) * right(2, 3))
      MM(2, 2) = CONJG(left(1, 2) * right(2, 1) + left(2, 2) * right(2, 2) + left(3, 2) * right(2, 3))
      MM(3, 2) = CONJG(left(1, 3) * right(2, 1) + left(2, 3) * right(2, 2) + left(3, 3) * right(2, 3))
      !# last index
      MM(1, 3) = CONJG(left(1, 1) * right(3, 1) + left(2, 1) * right(3, 2) + left(3, 1) * right(3, 3))
      MM(2, 3) = CONJG(left(1, 2) * right(3, 1) + left(2, 2) * right(3, 2) + left(3, 2) * right(3, 3))
      MM(3, 3) = CONJG(left(1, 3) * right(3, 1) + left(2, 3) * right(3, 2) + left(3, 3) * right(3, 3))
   END SUBROUTINE MultiplyMatdagMatdag

   PURE FUNCTION RealTraceMat(left) RESULT(trMM)
      COMPLEX(kind=WC), DIMENSION(3, 3), INTENT(IN) :: left
      REAL(kind=WP) :: TrMM
      !"""
      !# !Takes the real part of the trace of (3,3) complex numbers left
      ! Tr(left)
      !"""
      TrMM = real(left(1, 1) + left(2, 2) + left(3, 3), kind=WP)
   END FUNCTION RealTraceMat

   PURE FUNCTION TraceMat(left) RESULT(trMM)
      COMPLEX(kind=WC), DIMENSION(3, 3), INTENT(IN) :: left
      COMPLEX(kind=WC) :: TrMM
      !"""
      !# !Takes the trace of (3,3) complex numbers left
      ! Tr(left)
      !"""
      TrMM = left(1, 1) + left(2, 2) + left(3, 3)
   END FUNCTION TraceMat

   PURE SUBROUTINE TraceMultMatMat(TrMM, left, right)
      COMPLEX(kind=WC), DIMENSION(3, 3), INTENT(IN) :: left, right
      COMPLEX(kind=WC), INTENT(OUT) :: TrMM
      !"""
      !# !Takes the trace of (3,3) complex numbers left, right multiplied together
      !Tr(left*right)
      !"""
      TrMM = left(1, 1) * right(1, 1) + left(1, 2) * right(2, 1) + left(1, 3) * right(3, 1) + &
             left(2, 1) * right(1, 2) + left(2, 2) * right(2, 2) + left(2, 3) * right(3, 2) + &
             left(3, 1) * right(1, 3) + left(3, 2) * right(2, 3) + left(3, 3) * right(3, 3)
   END SUBROUTINE TraceMultMatMat

   PURE SUBROUTINE RealTraceMultMatMat(RTrMM, left, right)
      COMPLEX(kind=WC), DIMENSION(3, 3), INTENT(IN) :: left, right
      REAL(kind=WP), INTENT(OUT) :: RTrMM
      COMPLEX(kind=WC) :: TrMM
      !"""
      !# !Takes the real trace of (3,3) complex numbers left, right multiplied together
      !Real(Tr(left*right))
      !"""
      CALL TraceMultMatMat(TrMM, left, right)
      RTrMM = real(TrMM, kind=WP)
   END SUBROUTINE RealTraceMultMatMat

   PURE SUBROUTINE TraceLessConjgSubtract(TrSub, left, right)
      COMPLEX(kind=WC), DIMENSION(3, 3), INTENT(IN) :: left, right
      COMPLEX(kind=WC), DIMENSION(3, 3), INTENT(OUT) :: TrSub
      COMPLEX(kind=WC) :: trMM
      !"""
      !# Takes the traceless conjugate subtraction of A and B
      ! TrSub = A - B^dagger - Tr(A-B^dagger) / 3.0
      !
      TrSub = left - CONJG(TRANSPOSE(right))
      CALL TraceMultMatMat(TrMM, TrSub, Ident3x3)
      TrSub = TrSub - trMM / 3.0_WP
   END SUBROUTINE TraceLessConjgSubtract

   PURE SUBROUTINE colourDecomp(com, A)
      COMPLEX(kind=WC), DIMENSION(3, 3), INTENT(IN) :: A
      COMPLEX(kind=WC), DIMENSION(8), INTENT(OUT) :: com
      !"""
      !# Decompose the matrix M[][] into the SU(3) Gell-Mann components
      !"""
      com = (0.0_WP, 0.0_WP)
      com(1) = real(A(1, 2), kind=WP) + real(A(2, 1), kind=WP)
      com(2) = -AIMAG(A(1, 2)) + AIMAG(A(2, 1))
      com(3) = A(1, 1) - A(2, 2)
      com(4) = real(A(1, 3), kind=WP) + real(A(3, 1), kind=WP)
      com(5) = -AIMAG(A(1, 3)) + AIMAG(A(3, 1))
      com(6) = real(A(2, 3), kind=WP) + real(A(3, 2), kind=WP)
      com(7) = -AIMAG(A(2, 3)) + AIMAG(A(3, 2))
      com(8) = real(A(1, 1), kind=WP) + real(A(2, 2), kind=WP) - (2.0_WP / (3.0_WP**0.5_WP)) * real(A(3, 3), kind=WP)
   END SUBROUTINE colourDecomp

   PURE FUNCTION ExpIQ(Q) RESULT(V)
      COMPLEX(kind=WC), DIMENSION(3, 3), INTENT(IN) :: Q
      COMPLEX(kind=WC), DIMENSION(3, 3) :: V
      !"""
      !Compute the matrix exponential $V=exp(iQ)$ where Q is hermitian and traceless
      !"""
      REAL(kind=WP), PARAMETER :: eps = EPSILON(1.0_WP)
      COMPLEX(kind=WC), DIMENSION(3, 3) :: Q2, Q3
      REAL(kind=WP) :: c0, c1, c0Max
      REAL(kind=WP) :: w, u, w2, u2
      REAL(kind=WP) :: xi0, pm, theta
      COMPLEX(kind=WC), DIMENSION(3) :: h, f
      INTEGER :: jj  ! a counter
      ! First compute Q^2 and Q^3
      CALL MultiplyMatMat(Q2, Q, Q)
      CALL MultiplyMatMat(Q3, Q2, Q)
      ! Get the terms from the characteristic polynomial
      c0 = (1.0_WP / 3.0_WP) * RealTraceMat(Q3)
      c1 = 0.5_WP * RealTraceMat(Q2)
      ! handle sign of c0
      pm = SIGN(1.0_WP, real(c0, kind=WP))
      c0 = ABS(c0)
      ! Compute auxillary angle theta
      c0Max = 2.0_WP * (c1 / 3.0_WP)**(1.5_WP)
      IF (c0Max < eps) THEN
         theta = 0.0_WP
      ELSE
         ! Force the argument of acos to between [-1, 1]
         theta = ACOS(MIN(1.0_WP, MAX(-1.0_WP, c0 / c0Max)))
         !theta = acos(c0/c0Max)
      END IF
      ! Compute u and w
      u = SQRT(c1 / 3.0_WP) * COS(theta / 3.0_WP)
      w = SQRT(c1) * SIN(theta / 3.0_WP)
      u2 = u**2.0_WP
      w2 = w**2.0_WP
      ! Compute auxillary function xi0
      ! note this is not anisotropy
      ! USe Taylor expansion for small w to avoid cancellation errors
      IF (ABS(w) .LE. 0.05_WP) THEN
         xi0 = 1.0_WP - (1.0_WP / 6.0_WP) * w2 * (1.0_WP - (1.0_WP / 20.0_WP) * w2 * (1.0_WP - (1.0_WP / 42.0_WP) * w2))
      ELSE
         xi0 = SIN(w) / w
      END IF
      ! Compute the h-coefficients which are intermediate coefficients
      h(1) = (u2 - w2) * EXP(CMPLX(0.0_WP, 2.0_WP, kind=WC) * u) + &
             (8.0_WP * u2 * COS(w) + CMPLX(0.0_WP, 2.0_WP, kind=WC) * &
              u * (3.0_WP * u2 + w2) * xi0) * EXP(-CMPLX(0.0_WP, 1.0_WP, kind=WC) * u)
      h(2) = 2.0_WP * u * EXP(CMPLX(0.0_WP, 2.0_WP, kind=WC) * u) - &
             (2.0_WP * u * COS(w) - CMPLX(0.0_WP, 1.0_WP, kind=WC) * (3.0_WP * u2 - w2) * xi0) &
             * EXP(-CMPLX(0.0, 1.0_WP, kind=WC) * u)
      h(3) = EXP(CMPLX(0.0_WP, 2.0_WP, kind=WC) * u) - &
             (COS(w) + CMPLX(0.0_WP, 3.0_WP, kind=WC) * u * xi0) * EXP(-CMPLX(0.0_WP, 1.0_WP, kind=WC) * u)
      ! Normallise
      IF (real(c0max, kind=WP) < eps) THEN
         f(1) = 1.0_WP
         f(2:3) = 0.0_WP
      ELSE
         f(:) = h(:) / (9.0_WP * u2 - w2)
      END IF
      IF (pm < 0.0_WP) THEN
         DO jj = 1, 3
            f(jj) = ((-1.0_WP)**(jj - 1)) * CONJG(f(jj))
         END DO
      END IF
      ! Form the result
      V = f(1) * Ident3x3 + f(2) * Q + f(3) * Q2

   END FUNCTION ExpIQ

END MODULE FLUE_SU3MatrixOps
