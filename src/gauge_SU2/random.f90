MODULE FLUE_SU2_random
   USE FLUE_constants, ONLY: WP, WC
   USE FLUE_matrixConstants, ONLY: Ident2x2, sigma1, sigma2, sigma3
   IMPLICIT NONE(TYPE, EXTERNAL)
   PRIVATE
   PUBLIC :: constructSU2Matrix
CONTAINS

  PURE FUNCTION constructSU2Matrix(r) RESULT(U)
    REAL(kind=WP), DIMENSION(0:3), INTENT(IN) :: r
    COMPLEX(kind=WC), DIMENSION(2, 2) :: U
    REAL(kind=WP), DIMENSION(0:3) :: x
    REAL(kind=WP) :: normr
    ! U = x0 I + i xvec * sigmavec
    ! eqn 4.24
    !
    normr = SQRT(SUM(r**2))
    IF (normr <= TINY(1.0_WP)) THEN
       x = [1.0_WP, 0.0_WP, 0.0_WP, 0.0_WP]
    ELSE
       x = r / normr
    END IF
    U = x(0) * Ident2x2
    U = U + CMPLX(0.0_WP, 1.0_WP) * x(1) * sigma1
    U = U + CMPLX(0.0_WP, 1.0_WP) * x(2) * sigma2
    U = U + CMPLX(0.0_WP, 1.0_WP) * x(3) * sigma3
  END FUNCTION constructSU2Matrix

END MODULE FLUE_SU2_random
