! Some subroutines that act on SU(3) (3x3 complex) matrices
module FLUE_SU3MatrixOps
   use FLUE_constants, only: WC, WP
   implicit none
   ! Strictly this matrix visually is the tranpose in memory
   ! but as this is the identity it doesnt matter
   complex(kind=WC), dimension(3, 3), parameter :: Ident = reshape(source=[ &
                                                                   (1.0_wp, 0.0_wp), (0.0_wp, 0.0_wp), (0.0_wp, 0.0_wp), &
                                                                   (0.0_wp, 0.0_wp), (1.0_wp, 0.0_wp), (0.0_wp, 0.0_wp), &
                                                                (0.0_wp, 0.0_wp), (0.0_wp, 0.0_wp), (1.0_wp, 0.0_wp)], shape=[3, 3])
   !private

   !public :: TracelessConjgSubtract, colourDecomp

contains

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

   pure subroutine MultiplyMatdagMatdag(MM, left, right)
      complex(kind=WC), dimension(3, 3), intent(in) :: left, right
      complex(kind=WC), dimension(3, 3), intent(out) :: MM
      !"""
      !#Multiplies two (3,3) complex matrices together. Takes conjugate
      !Does (left*right)^dagger
      !"""
      !# take transpose manually
      MM(1, 1) = conjg(left(1, 1) * right(1, 1) + left(2, 1) * right(1, 2) + left(3, 1) * right(1, 3))
      MM(2, 1) = conjg(left(1, 2) * right(1, 1) + left(2, 2) * right(1, 2) + left(3, 2) * right(1, 3))
      MM(3, 1) = conjg(left(1, 3) * right(1, 1) + left(2, 3) * right(1, 2) + left(3, 3) * right(1, 3))
      !# but take conjugate using np
      MM(1, 2) = conjg(left(1, 1) * right(2, 1) + left(2, 1) * right(2, 2) + left(3, 1) * right(2, 3))
      MM(2, 2) = conjg(left(1, 2) * right(2, 1) + left(2, 2) * right(2, 2) + left(3, 2) * right(2, 3))
      MM(3, 2) = conjg(left(1, 3) * right(2, 1) + left(2, 3) * right(2, 2) + left(3, 3) * right(2, 3))
      !# last index
      MM(1, 3) = conjg(left(1, 1) * right(3, 1) + left(2, 1) * right(3, 2) + left(3, 1) * right(3, 3))
      MM(2, 3) = conjg(left(1, 2) * right(3, 1) + left(2, 2) * right(3, 2) + left(3, 2) * right(3, 3))
      MM(3, 3) = conjg(left(1, 3) * right(3, 1) + left(2, 3) * right(3, 2) + left(3, 3) * right(3, 3))
   end subroutine MultiplyMatdagMatdag

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
      TrSub = left - conjg(transpose(right))
      call TraceMultMatMat(TrMM, TrSub, Ident)
      TrSub = TrSub - trMM / 3.0_wp
   end subroutine TraceLessConjgSubtract

   pure subroutine colourDecomp(com, A)
      complex(kind=WC), dimension(3, 3), intent(in) :: A
      complex(kind=WC), dimension(8), intent(out) :: com
      !"""
      !# Decompose the matrix M[][] into the SU(3) Gell-Mann components
      !"""
      com = (0.0_wp, 0.0_wp)
      com(1) = real(A(1, 2), kind=WP) + real(A(2, 1), kind=WP)
      com(2) = -aimag(A(1, 2)) + aimag(A(2, 1))
      com(3) = A(1, 1) - A(2, 2)
      com(4) = real(A(1, 3), kind=WP) + real(A(3, 1), kind=WP)
      com(5) = -aimag(A(1, 3)) + aimag(A(3, 1))
      com(6) = real(A(2, 3), kind=WP) + real(A(3, 2), kind=WP)
      com(7) = -aimag(A(2, 3)) + aimag(A(3, 2))
      com(8) = real(A(0, 0), kind=WP) + real(A(2, 2), kind=WP) - (2.0_wp / (3.0_wp**0.5_wp)) * real(A(3, 3), kind=WP)
   end subroutine colourDecomp
end module FLUE_SU3MatrixOps
