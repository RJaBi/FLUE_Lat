
MODULE FLUE_SU3MatrixOps
  !< Module FLUE_SU3MatrixOps
  !< Purpose
  !<   Provides high-performance, pure procedures for manipulating SU(3)
  !<   matrices and vectors in the context of lattice gauge theory.
  !<
  !< Physical Context:
  !<   In lattice QCD, the fundamental dynamical variables are SU(3)-valued
  !<   link matrices:
  !<
  !<       U_mu(x) ∈ SU(3)
  !<
  !<   representing parallel transporters of the colour gauge field between
  !<   neighbouring lattice sites. These matrices are unitary and have
  !<   determinant equal to 1.
  !<
  !<   This module supports:
  !<     - Algebra on SU(3) group elements (matrix multiplication, traces)
  !<     - Operations in the Lie algebra su(3) (traceless Hermitian matrices)
  !<     - Projection/reunitarisation of matrices back onto SU(3)
  !<     - Conversion between matrix and generator (Gell-Mann) representations
  !<     - Efficient computation of exp(iQ), used in gauge updates
  !<
  !< Typical Use in Simulations:
  !<   - Gauge field updates (heatbath, overrelaxation)
  !<   - Projection of numerically drifted matrices back to SU(3)
  !<   - Observables such as Wilson loops involving traces of products
  !<
  !< Mathematical Conventions:
  !<   - Matrices are 3x3 complex arrays.
  !<   - Hermitian conjugate: A^\dagger = CONJG(TRANSPOSE(A))
  !<   - Trace: Tr(A)
  !<   - Generators: λ_i (Gell-Mann matrices)
  !<
  !< Dependencies:
  !<   - FLUE_constants: working precisions WC (complex), WP (real)
  !<   - FLUE_matrixConstants: identity matrix Ident3x3
  !<
  !< Notes:
  !<   - All procedures are pure and suitable for parallel execution.
  !<   - Explicitly unrolled operations improve performance on CPUs/GPUs.
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


  pure subroutine orthogonalise_vectors(w, v)
    !< Orthogonalise a colour vector against another.
    !<
    !< Physical Interpretation:
    !<   Used during reunitarisation of SU(3) link matrices, where rows (or columns)
    !<   must form an orthonormal basis in complex colour space.
    !<
    !< Mathematics:
    !<   $w \leftarrow w - v (v^\dagger w)$
    !<
    complex(kind=WC), dimension(3), intent(INOUT) :: w
    !<       Colour vector to be orthogonalised.
    complex(kind=WC), dimension(3), intent(IN) :: v
    !<       Reference vector (typically already normalised).
    complex(kind=WC) :: vdotw
    !< sum(conjg(v) * w)
      vdotw = SUM(CONJG(v) * w)
      w = w - v * vdotw
   END SUBROUTINE orthogonalise_vectors

   pure subroutine vector_product(x, v, w)
     !< Compute SU(3) vector product.
     !<
     !< Physical Interpretation:
     !<   Constructs a third orthogonal colour vector from two existing ones,
     !<   analogous to a cross product, ensuring a complete orthonormal basis.
     !<
     !< Mathematics:
     !<   $x_i = conjugate(v_j w_k - v_k w_j)$, cyclic indices.
     !<
     complex(kind=WC), dimension(3), intent(OUT) :: x
     !<       Resulting orthogonal vector.
     complex(kind=WC), dimension(3), intent(IN) :: v, w
     !<   v,w : input colour vectors.
     integer :: ic, jc, kc
     !< counters
     integer, parameter :: nc = 3
     !< number of colours for SU(3) (3)
      do ic = 1, nc
         jc = MODULO(ic, 3) + 1
         kc = MODULO(jc, 3) + 1
         x(ic) = CONJG(v(jc) * w(kc) - v(kc) * w(jc))
      END DO
   END SUBROUTINE vector_product

   pure subroutine FixSU3Matrix(U_x)
     !< Project a 3x3 matrix back onto SU(3).
     !<
     !< Physical Interpretation:
     !<   Numerical updates can cause
     !<   link variables U_mu(x) to drift away from SU(3). This routine restores:
     !<
     !<     - Unitarity: U U^\dagger = I
     !<     - Determinant: det(U) ≈ 1
     !<
     !< Method:
     !<   - Gram-Schmidt orthonormalisation of rows
     !<   - Third row constructed from SU(3) vector product
     !<
     !< Usage:
     !<   Typically called after each gauge update or smearing step.
     complex(kind=WC), dimension(3, 3), intent(INOUT) :: U_x
     !<         Gauge link matrix, overwritten with SU(3)-projected version.
     complex(kind=WC), dimension(3) :: v1, v2, v3
     !< rows of SU(3) matrix
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

   pure subroutine normalise_vector(v)
     !< Normalise a colour vector.
     !<
     !< Physical Interpretation:
     !<   Ensures unit norm of vectors forming SU(3) matrix rows/columns.
     !<
     !< Mathematics:
     !<   v <- v / sqrt(sum(|v_i|^2))
     complex(kind=WC), dimension(3), intent(INOUT) :: v
     !< vector which is normallised
     real(kind=WP) :: norm
     !< The value to normallise by
      norm = SQRT(SUM(real(v)**2 + AIMAG(v)**2))
      v = v / norm
   END SUBROUTINE normalise_vector

   pure subroutine MultiplyMatMat(MM, left, right)
     !< Multiply two SU(3) matrices.
     !<
     !< Physical Interpretation:
     !<   Represents sequential parallel transport of colour fields:
     !<
     !<     U_total = U_1 * U_2
     !<
     !<   Common in staples, Wilson loops, and path-ordered products.
     !<
     !< Arguments:
     !<   left, right : input SU(3) matrices
     !<
     !<   MM          : result
     complex(kind=WC), dimension(3, 3), intent(IN) :: left, right
     !< input SU(3) matrices
     complex(kind=WC), dimension(3, 3), intent(OUT) :: MM
     !< resulting matrix
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
     !< Multiply matrix by Hermitian conjugate.
     !<
     !< Physical Interpretation:
     !<   Appears in gauge-invariant quantities such as:
     !<
     !<     U(x, μ) U^\dagger(x, μ)
     !<
     !< Computes:
     !<   MM = left * right^\dagger
     complex(kind=WC), dimension(3, 3), intent(IN) :: left, right
     !< Input SU(3) matrices
     complex(kind=WC), dimension(3, 3), intent(OUT) :: MM
     !< Resulting matrix
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

   pure subroutine MultiplyMatdagMatdag(MM, left, right)
     !< Compute dagger of matrix product.
     !<
     !< Physical Interpretation:
     !<   Used when reversing direction of transport paths in gauge loops.
     !<
     !< Computes:
     !<   MM = (left * right)^\dagger
     complex(kind=WC), dimension(3, 3), intent(IN) :: left, right
     !< Input SU(3) matrices
     complex(kind=WC), dimension(3, 3), intent(OUT) :: MM
     !< Resulting Matrix
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

   pure function RealTraceMat(left) result(trMM)
     !< Real(Trace of a matrix.)
     !<
     !< Physical Interpretation:
     !<   Central observable in lattice QCD:
     !<
     !<     - Action terms: Re Tr(U_plaquette)
     !<
     !< Returns:
     !<   Real trace.
     complex(kind=WC), dimension(3, 3), intent(IN) :: left
     !< Input SU(3) matrix
     real(kind=WP) :: TrMM
     !< the output real trace of the input matrix
     !< Output ReTrace
      !"""
      !# !Takes the real part of the trace of (3,3) complex numbers left
      ! Tr(left)
      !"""
      TrMM = real(left(1, 1) + left(2, 2) + left(3, 3), kind=WP)
   END FUNCTION RealTraceMat

   pure function TraceMat(left) result(trMM)
     !< Trace of a matrix.
     !<
     !< Physical Interpretation:
     !<   Central observable in lattice QCD:
     !<
     !<     - Wilson loops: Tr(path-ordered product)
     !<     - Action terms: Re Tr(U_plaquette)
     complex(kind=WC), dimension(3, 3), intent(IN) :: left
     !< Input SU(3) matrix
     complex(kind=WC) :: TrMM
     !< complex trace output
      !"""
      !# !Takes the trace of (3,3) complex numbers left
      ! Tr(left)
      !"""
      TrMM = left(1, 1) + left(2, 2) + left(3, 3)
   END FUNCTION TraceMat

   pure subroutine TraceMultMatMat(TrMM, left, right)
     !< Trace of a product of matrices.
     !<
     !< Physical Interpretation:
     !<   Used in evaluating gauge-invariant observables such as Wilson loops.
     !<
     !< Computes:
     !<   Tr(left * right)
     complex(kind=WC), dimension(3, 3), intent(IN) :: left, right
     !< Input SU(3) matrices
     complex(kind=WC), intent(OUT) :: TrMM
     !< Resulting (complex) trace
      !"""
      !# !Takes the trace of (3,3) complex numbers left, right multiplied together
      !Tr(left*right)
      !"""
      TrMM = left(1, 1) * right(1, 1) + left(1, 2) * right(2, 1) + left(1, 3) * right(3, 1) + &
             left(2, 1) * right(1, 2) + left(2, 2) * right(2, 2) + left(2, 3) * right(3, 2) + &
             left(3, 1) * right(1, 3) + left(3, 2) * right(2, 3) + left(3, 3) * right(3, 3)
   END SUBROUTINE TraceMultMatMat

   pure subroutine RealTraceMultMatMat(RTrMM, left, right)
     !< Real trace of matrix product.
     !<
     !< Physical Interpretation:
     !<   Direct contribution to gauge action and observables.
     !<
     !< Computes:
     !<   Real(Tr(left * right))
     complex(kind=WC), dimension(3, 3), intent(IN) :: left, right
     !< Input SU(3) Matrices
     real(kind=WP), intent(OUT) :: RTrMM
     !< Resulting ReTrace
     complex(kind=WC) :: TrMM
     !< full complex trace of left * right
      !"""
      !# !Takes the real trace of (3,3) complex numbers left, right multiplied together
      !Real(Tr(left*right))
      !"""
      CALL TraceMultMatMat(TrMM, left, right)
      RTrMM = real(TrMM, kind=WP)
   END SUBROUTINE RealTraceMultMatMat

   pure subroutine TraceLessConjgSubtract(TrSub, left, right)
     !< Traceless anti-Hermitian projection.
     !<
     !< Physical Interpretation:
     !<   Maps matrices into the Lie algebra su(3), required for:
     !<     - Gauge updates
     !<
     !< Computes:
     !<   TrSub = left - right^\dagger
     !<   TrSub <- TrSub - (1/3) Tr(TrSub) I
     !<
     !< Result:
     !<   Traceless matrix in su(3).
     complex(kind=WC), dimension(3, 3), intent(IN) :: left, right
     !< Input SU(3) matrices
     complex(kind=WC), dimension(3, 3), intent(OUT) :: TrSub
     !< Resulting traceless conjgugate subtraction matrix
     complex(kind=WC) :: trMM
     !< trace of left - dagger(right)
      !"""
      !# Takes the traceless conjugate subtraction of A and B
      ! TrSub = A - B^dagger - Tr(A-B^dagger) / 3.0
      !
      TrSub = left - CONJG(TRANSPOSE(right))
      CALL TraceMultMatMat(TrMM, TrSub, Ident3x3)
      TrSub = TrSub - trMM / 3.0_WP
   END SUBROUTINE TraceLessConjgSubtract


   pure subroutine colourDecomp(com, A)
     !< Decompose matrix into Gell-Mann basis.
     !<
     !< Physical Interpretation:
     !<   Converts matrix representation into Lie algebra coefficients:
     !<
     !<     A = sum_i com(i) λ_i
     !<
     !<   Useful for:
     !<     - Analysing gauge fields
     !<     - Constructing Lie algebra updates
     !<     - Diagnostics and measurements
     !<
     !< Arguments:
     !<   A   : SU(3) or general complex matrix
     !<
     !<   com : coefficients in Gell-Mann basis (8 generators)
     complex(kind=WC), dimension(3, 3), intent(IN) :: A
     !< Input SU(3) matrice
     complex(kind=WC), dimension(8), intent(OUT) :: com
     !< Output coefficients in Gell-Mann basis
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

   pure function ExpIQ(Q) result(V)
     !< Compute exponential of Lie algebra element.
     !<
     !< Physical Interpretation:
     !<   Core operation in gauge field evolution:
     !<
     !<     U_new = exp(i Q) * U_old
     !<
     !<   where:
     !<     Q ∈ su(3) is Hermitian and traceless (momentum field).
     !<
     !<   Used in:
     !<     - Heatbath
     !<
     !< Algorithm:
     !<   Uses Cayley-Hamilton theorem to avoid explicit diagonalisation:
     !<
     !<     exp(iQ) = f1 * I + f2 * Q + f3 * Q^2
     !<
     !<   with coefficients determined from invariants:
     !<     c0 = (1/3) Tr(Q^3)
     !<     c1 = (1/2) Tr(Q^2)
     !<
     !< Numerical Stability:
     !<   - Uses Taylor expansion for small eigenvalues
     !<   - Clamps acos argument for safety
     !<   - Handles sign structure explicitly
     !<
     !< Arguments:
     !<   Q : Hermitian, traceless matrix (su(3) element)
     !<
     !< Returns:
     !<   V : SU(3) matrix
     !<
     !< Notes:
     !<   - Produces unitary matrix up to numerical precision.
     !<   - Avoids costly eigen-decomposition.
     complex(kind=WC), dimension(3, 3), intent(IN) :: Q
     !< Input Hermitian, traceless (su(3) element)
     complex(kind=WC), dimension(3, 3) :: V
     !< Output matrix exponential
      !"""
      !Compute the matrix exponential $V=exp(iQ)$ where Q is hermitian and traceless
      !"""
     REAL(kind=WP), PARAMETER :: eps = EPSILON(1.0_WP)
     !< smallest number in WP
      COMPLEX(kind=WC), DIMENSION(3, 3) :: Q2, Q3
      !< Matrices for Cayley Hamilton theorem
      REAL(kind=WP) :: c0, c1, c0Max
      !< for Cayley-Hamilton theorem
      REAL(kind=WP) :: w, u, w2, u2
      !< used in coefficents
      REAL(kind=WP) :: xi0, pm, theta
      !< auxillary function xi0, sign of c0, auxillary angle theta
      COMPLEX(kind=WC), DIMENSION(3) :: h, f
      !< intermediate and final coefficients
      INTEGER :: jj
      !< a counter
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
