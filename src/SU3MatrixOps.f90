! Some subroutines that act upon and create SU(3) (3x3 complex) matrices
module FLUE_SU3MatrixOps
  use FLUE_constants, only: WC, WP
  use FLUE_matrixConstants, only: Ident3x3
  implicit none(type, external)
  private
  public :: MultiplyMatMat, MultiplyMatDagMatDag, &
       TraceMultMatMat, RealTraceMultMatMat, TracelessConjgSubtract, &
       colourDecomp, RealTraceMat, MultiplyMatMatDag, TraceMat
  public :: FixSU3Matrix
  public :: orthogonalise_vectors, vector_product
  public :: ExpIQ

contains

     ! stripped from cola and de-colour vectored
   pure subroutine orthogonalise_vectors(w, v)
      complex(kind=WC), dimension(3), intent(inout) :: w
      complex(kind=WC), dimension(3), intent(in) :: v
      complex(kind=WC) :: vdotw
      vdotw = SUM(CONJG(v) * w)
      w = w - v * vdotw
   end subroutine orthogonalise_vectors
   ! stripped from cola and de-colour vectored
   pure subroutine vector_product(x, v, w)
      complex(kind=WC), dimension(3), intent(out) :: x
      complex(kind=WC), dimension(3), intent(in) :: v, w
      integer :: ic, jc, kc
      integer, parameter :: nc = 3
      do ic = 1, nc
         jc = MODULO(ic, 3) + 1
         kc = MODULO(jc, 3) + 1
         x(ic) = CONJG(v(jc) * w(kc) - v(kc) * w(jc))
      end do
   end subroutine vector_product
   ! stripped from cola and de-colour vectored
   pure subroutine FixSU3Matrix(U_x)
      complex(kind=WC), dimension(3, 3), intent(inout) :: U_x
      complex(kind=WC), dimension(3) :: v1, v2, v3
      v1 = U_x(1, :)
      call normalise_vector(v1)
      v2(:) = U_x(2, :)
      call orthogonalise_vectors(v2, v1)
      call normalise_vector(v2)
      call vector_product(v3, v1, v2)
      call normalise_vector(v3)
      U_x(1, :) = v1(:)
      U_x(2, :) = v2(:)
      U_x(3, :) = v3(:)
   end subroutine FixSU3Matrix
   ! stripped from cola and de-colour vectored
   pure subroutine normalise_vector(v)
      complex(kind=WC), dimension(3), intent(inout) :: v
      real(WP) :: norm
      norm = SQRT(SUM(real(v)**2 + AIMAG(v)**2))
      v = v / norm
    end subroutine normalise_vector
    ! explicitly unroll 3x3 matrix mult
   pure subroutine MultiplyMatMat(MM, left, right)
      complex(kind=WC), dimension(3, 3), intent(in) :: left, right
      complex(kind=WC), dimension(3, 3), intent(out) :: MM
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
   end subroutine MultiplyMatMat


   pure subroutine MultiplyMatMatDag(MM, left, right)
   ! A B^\dag
      complex(kind=WC), dimension(3, 3), intent(in) :: left, right
      complex(kind=WC), dimension(3, 3), intent(out) :: MM
      !"""
      !Multiple left by right^\dagger. Assumes 3x3 (colour) (complex) matrices
      !"""
      !# do the maths for the colour matrices
      MM(1, 1) = left(1, 1) * conjg(right(1, 1)) + left(1, 2) * conjg(right(1, 2)) + left(1, 3) * conjg(right(1, 3))
      MM(2, 1) = left(2, 1) * conjg(right(1, 1)) + left(2, 2) * conjg(right(1,2)) + left(2, 3) * conjg(right(1, 3))
      MM(3, 1) = left(3, 1) * conjg(right(1, 1)) + left(3, 2) * conjg(right(1,2)) + left(3, 3) * conjg(right(1, 3))
      !# second index
      MM(1, 2) = left(1, 1) * conjg(right(2, 1)) + left(1, 2) * conjg(right(2, 2)) + left(1, 3) * conjg(right(2, 3))
      MM(2, 2) = left(2, 1) * conjg(right(2, 1)) + left(2, 2) * conjg(right(2, 2)) + left(2, 3) * conjg(right(2, 3))
      MM(3, 2) = left(3, 1) * conjg(right(2, 1)) + left(3, 2) * conjg(right(2, 2)) + left(3, 3) * conjg(right(2, 3))
      !# third index
      MM(1, 3) = left(1, 1) * conjg(right(3, 1)) + left(1, 2) * conjg(right(3,2 )) + left(1, 3) * conjg(right(3, 3))
      MM(2, 3) = left(2, 1) * conjg(right(3, 1)) + left(2, 2) * conjg(right(3,2)) + left(2, 3) * conjg(right(3, 3))
      MM(3, 3) = left(3, 1) * conjg(right(3, 1)) + left(3, 2) * conjg(right(3,2)) + left(3, 3) * conjg(right(3, 3))
   end subroutine MultiplyMatMatDag

   pure subroutine MultiplyMatdagMatdag(MM, left, right)
      complex(kind=WC), dimension(3, 3), intent(in) :: left, right
      complex(kind=WC), dimension(3, 3), intent(out) :: MM
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
   end subroutine MultiplyMatdagMatdag

   pure function RealTraceMat(left) result(trMM)
      complex(kind=WC), dimension(3, 3), intent(in) :: left
      real(kind=WP):: TrMM
      !"""
      !# !Takes the real part of the trace of (3,3) complex numbers left
      ! Tr(left)
      !"""
      TrMM = real(left(1, 1) + left(2, 2) + left(3, 3), kind=WP)
end function RealTraceMat

   pure function TraceMat(left) result(trMM)
      complex(kind=WC), dimension(3, 3), intent(in) :: left
      complex(kind=WC):: TrMM
      !"""
      !# !Takes the trace of (3,3) complex numbers left
      ! Tr(left)
      !"""
      TrMM = left(1, 1) + left(2, 2) + left(3, 3)
end function TraceMat


   pure subroutine TraceMultMatMat(TrMM, left, right)
      complex(kind=WC), dimension(3, 3), intent(in) :: left, right
      complex(kind=WC), intent(out) :: TrMM
      !"""
      !# !Takes the trace of (3,3) complex numbers left, right multiplied together
      !Tr(left*right)
      !"""
      TrMM = left(1, 1) * right(1, 1) + left(1, 2) * right(2, 1) + left(1, 3) * right(3, 1) + &
             left(2, 1) * right(1, 2) + left(2, 2) * right(2, 2) + left(2, 3) * right(3, 2) + &
             left(3, 1) * right(1, 3) + left(3, 2) * right(2, 3) + left(3, 3) * right(3, 3)
   end subroutine TraceMultMatMat

   pure subroutine RealTraceMultMatMat(RTrMM, left, right)
      complex(kind=WC), dimension(3, 3), intent(in) :: left, right
      real(kind=WP), intent(out) :: RTrMM
      complex(kind=WC) :: TrMM
      !"""
      !# !Takes the real trace of (3,3) complex numbers left, right multiplied together
      !Real(Tr(left*right))
      !"""
      call TraceMultMatMat(TrMM, left, right)
      RTrMM = real(TrMM, kind=WP)
   end subroutine RealTraceMultMatMat

   pure subroutine TraceLessConjgSubtract(TrSub, left, right)
      complex(kind=WC), dimension(3, 3), intent(in) :: left, right
      complex(kind=WC), dimension(3, 3), intent(out) :: TrSub
      complex(kind=WC) :: trMM
      !"""
      !# Takes the traceless conjugate subtraction of A and B
      ! TrSub = A - B^dagger - Tr(A-B^dagger) / 3.0
      !
      TrSub = left - CONJG(TRANSPOSE(right))
      call TraceMultMatMat(TrMM, TrSub, Ident3x3)
      TrSub = TrSub - trMM / 3.0_WP
   end subroutine TraceLessConjgSubtract

   pure subroutine colourDecomp(com, A)
      complex(kind=WC), dimension(3, 3), intent(in) :: A
      complex(kind=WC), dimension(8), intent(out) :: com
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
   end subroutine colourDecomp

   pure function ExpIQ(Q) result(V)
      complex(kind=WC), dimension(3,3), intent(in) :: Q
      complex(kind=WC), dimension(3,3) :: V
      !"""
      !Compute the matrix exponential $V=exp(iQ)$ where Q is hermitian
      !"""
      real(kind=WP), parameter :: eps = epsilon(1.0_WP)
      complex(kind=WC), dimension(3,3) :: Q2, Q3
      real(kind=WP) :: c0, c1, c0Max
      real(kind=WP) :: w, u, w2, u2
      real(kind=WP) :: xi0, pm, theta
      complex(kind=WC), dimension(3) :: h, f
      integer :: jj  ! a counter
      ! First compute Q^2 and Q^3
      call MultiplyMatMat(Q2, Q, Q)
      call MultiplyMatMat(Q3, Q2, Q)
      ! Get the terms from the characteristic polynomial
      c0 = (1.0_WP/3.0_WP) * RealTraceMat(Q3)
      c1 = 0.5_WP * RealTraceMat(Q2)
      ! handle sign of c0
      pm = sign(1.0_WP, real(c0, kind=WP))
      c0 = abs(c0)
      ! Compute auxillary angle theta
      c0Max = 2.0_WP * (c1 / 3.0_WP)**(1.5_WP)
      if (c0Max < eps) then
         theta = 0.0_WP
      else
         ! Force the argument of acos to between [-1, 1]
         theta = acos(min(1.0_WP, max(-1.0_WP, c0/c0Max)))
         !theta = acos(c0/c0Max)
      end if
      ! Compute u and w
      u = sqrt(c1/3.0_WP) * cos(theta / 3.0_WP)
      w = sqrt(c1) * sin(theta / 3.0_WP)
      u2 = u**2.0_WP
      w2 = w**2.0_WP
      ! Compute auxillary function xi0
      ! note this is not anisotropy
      ! USe Taylor expansion for small w to avoid cancellation errors
      if (abs(w) .le. 0.05_WP) then
         xi0 = 1.0_WP - (1.0_WP/6.0_WP) * w2 * (1.0_WP - (1.0_WP/20.0_WP) * w2 * (1.0_WP - (1.0_WP / 42.0_WP) * w2))
      else
         xi0 = sin(w) / w
      end if
      ! Compute the h-coefficients which are intermediate coefficients
      h(1) = (u2 - w2)*exp(cmplx(0.0_WP,2.0_WP,kind=WC)*u) + &
          (8.0_WP*u2*cos(w) + cmplx(0.0_WP,2.0_WP,kind=WC)*u*(3.0_WP*u2 + w2)*xi0)*exp(-cmplx(0.0_WP,1.0_WP,kind=WC)*u)
      h(2) = 2.0_WP*u*exp(cmplx(0.0_WP,2.0_WP,kind=WC)*u) - &
          (2.0_WP*u*cos(w) - cmplx(0.0_WP,1.0_WP,kind=WC)*(3.0_WP*u2 - w2)*xi0)*exp(-cmplx(0.0,1.0_WP,kind=WC)*u)
      h(3) = exp(cmplx(0.0_WP,2.0_WP,kind=WC)*u) - &
          (cos(w) + cmplx(0.0_WP,3.0_WP,kind=WC)*u*xi0)*exp(-cmplx(0.0_WP,1.0_WP,kind=WC)*u)
      ! Normallise
      if ( real(c0max, kind=WP) < eps ) then
       f(1) = 1.0_WP
       f(2:3) = 0.0_WP
    else
       f(:) = h(:)/(9.0_WP*u2-w2)
    end if
    if (pm < 0.0_WP) then
      do jj=1, 3
         f(jj) = ((-1.0_WP)**(jj-1)) * conjg(f(jj))
      end do
   end if
   ! Form the result
   V = f(1) * Ident3x3 + f(2) * Q + f(3) * Q2

end function ExpIQ


end module FLUE_SU3MatrixOps
