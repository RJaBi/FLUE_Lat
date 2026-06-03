MODULE FLUE_heatbath
   USE FLUE_constants, ONLY: WP, WC
   USE FLUE_matrixConstants, ONLY: Ident3x3
   USE FLUE_SU2_heatbath, ONLY: constructXMatrix
   USE FLUE_SU2_random, ONLY: constructSU2Matrix
   USE FLUE_SU3MatrixOps, ONLY: MultiplyMatMat, FixSU3Matrix
   USE FLUE_wloops, ONLY: genericPath, periodCoord
   IMPLICIT NONE(TYPE, EXTERNAL)
   PRIVATE
   PUBLIC :: updateLinks
CONTAINS

   SUBROUTINE updateLinks(U, beta, UUpdated)
      COMPLEX(kind=WC), DIMENSION(:, :, :, :, :, :, :), INTENT(IN) :: U
      REAL(kind=WP), INTENT(IN) :: beta
      COMPLEX(kind=WC), DIMENSION(:, :, :, :, :, :, :), INTENT(INOUT) :: UUpdated
    !! lattice geometry
      INTEGER, DIMENSION(7) :: dataShape
      INTEGER :: nt, nx, ny, nz
    !! counters
      INTEGER :: ix, iy, iz, it, mu
      INTEGER, DIMENSION(4) :: coord
    !! matrices
      COMPLEX(kind=WC), DIMENSION(3, 3) :: staple, W, embed, ULink, UTemp
      COMPLEX(kind=WC), DIMENSION(2, 2) :: SU2_M, MfromW, SU2_X, SU2_U
    !! For SU2 quarternion
      REAL(kind=WP), DIMENSION(0:3) :: aQuart
      REAL(kind=WP) :: alpha

      dataShape = SHAPE(U)
      nt = dataShape(4)
      nx = dataShape(5)
      ny = dataShape(6)
      nz = dataShape(7)

      UUpdated = U
      DO it = 1, nt
         DO ix = 1, nx
            DO iy = 1, ny
               DO iz = 1, nz
                  coord = (/it, ix, iy, iz/)
                  DO mu = 1, 4
                     ULink = UUpdated(:, :, mu, it, ix, iy, iz)
                     CALL stapleAt(UUpdated, staple, coord, mu)
                     ! W = U*V
                     CALL MultiplyMatMat(W, ULink, staple)
                     ! Do the first SU2 update in 1,2
                     MfromW(1, 1) = W(1, 1)
                     MfromW(2, 2) = W(2, 2)
                     MfromW(1, 2) = W(1, 2)
                     MfromW(2, 1) = W(2, 1)
                     ! Project this to SU2 proper
                     aQuart(0) = real(MfromW(1, 1) + MfromW(2, 2), kind=WP)
                     aQuart(1) = AIMAG(MfromW(1, 2) + MfromW(2, 1))
                     aQuart(2) = real(MfromW(1, 2) - MfromW(2, 1), kind=WP)
                     aQuart(3) = AIMAG(MfromW(1, 1) - MfromW(2, 2))
                     aQuart = aQuart * 0.5_WP
                     alpha = SQRT(SUM(aQuart**2.0))
                     aQuart = aQuart / alpha
                     SU2_M = constructSU2Matrix(aQuart)
                     ! Do an SU2 Update
                     SU2_X = constructXMatrix(alpha, 2.0_WP * beta / 3.0_WP)
                     ! Take dagger
                     SU2_M = CONJG(TRANSPOSE(SU2_M))
                     ! XV^dag
                     SU2_U = MATMUL(SU2_X, SU2_M)
                     ! Embed the SU2 matrix in SU3
                     embed = Ident3x3
                     embed(1, 1) = SU2_U(1, 1)
                     embed(2, 2) = SU2_U(2, 2)
                     embed(1, 2) = SU2_U(1, 2)
                     embed(2, 1) = SU2_U(2, 1)
                     ! embed * U
                     CALL MultiplyMatMat(Utemp, embed, ULink)
                     ULink = Utemp
                     ! Update W
                     CALL MultiplyMatMat(Utemp, embed, W)
                     W = UTemp
                   !!! Next do 1,3
                     MfromW(1, 1) = W(1, 1)
                     MfromW(2, 2) = W(3, 3)
                     MfromW(1, 2) = W(1, 3)
                     MfromW(2, 1) = W(3, 1)
                     ! Project this to SU2 proper
                     aQuart(0) = real(MfromW(1, 1) + MfromW(2, 2), kind=WP)
                     aQuart(1) = AIMAG(MfromW(1, 2) + MfromW(2, 1))
                     aQuart(2) = real(MfromW(1, 2) - MfromW(2, 1), kind=WP)
                     aQuart(3) = AIMAG(MfromW(1, 1) - MfromW(2, 2))
                     aQuart = aQuart * 0.5_WP
                     alpha = SQRT(SUM(aQuart**2.0))
                     aQuart = aQuart / alpha
                     SU2_M = constructSU2Matrix(aQuart)
                     ! Do an SU2 Update
                     SU2_X = constructXMatrix(alpha, 2.0_WP * beta / 3.0_WP)
                     ! Take dagger
                     SU2_M = CONJG(TRANSPOSE(SU2_M))
                     ! XV^dag
                     SU2_U = MATMUL(SU2_X, SU2_M)
                     ! Embed the SU2 matrix in SU3
                     embed = Ident3x3
                     embed(1, 1) = SU2_U(1, 1)
                     embed(3, 3) = SU2_U(2, 2)
                     embed(1, 3) = SU2_U(1, 2)
                     embed(3, 1) = SU2_U(2, 1)
                     ! embed * U
                     CALL MultiplyMatMat(Utemp, embed, ULink)
                     ULink = Utemp
                     ! Update W
                     CALL MultiplyMatMat(Utemp, embed, W)
                     W = UTemp
                   !!! And now do 2,3
                     MfromW(1, 1) = W(2, 2)
                     MfromW(2, 2) = W(3, 3)
                     MfromW(1, 2) = W(2, 3)
                     MfromW(2, 1) = W(3, 2)
                     ! Project this to SU2 proper
                     aQuart(0) = real(MfromW(1, 1) + MfromW(2, 2), kind=WP)
                     aQuart(1) = AIMAG(MfromW(1, 2) + MfromW(2, 1))
                     aQuart(2) = real(MfromW(1, 2) - MfromW(2, 1), kind=WP)
                     aQuart(3) = AIMAG(MfromW(1, 1) - MfromW(2, 2))
                     aQuart = aQuart * 0.5_WP
                     alpha = SQRT(SUM(aQuart**2.0))
                     aQuart = aQuart / alpha
                     SU2_M = constructSU2Matrix(aQuart)
                     ! Do an SU2 Update
                     SU2_X = constructXMatrix(alpha, 2.0_WP * beta / 3.0_WP)
                     ! Take dagger
                     SU2_M = CONJG(TRANSPOSE(SU2_M))
                     ! XV^dag
                     SU2_U = MATMUL(SU2_X, SU2_M)
                     ! Embed the SU2 matrix in SU3
                     embed = Ident3x3
                     embed(2, 2) = SU2_U(1, 1)
                     embed(3, 3) = SU2_U(2, 2)
                     embed(2, 3) = SU2_U(1, 2)
                     embed(3, 2) = SU2_U(2, 1)
                     ! embed * U
                     CALL MultiplyMatMat(Utemp, embed, ULink)
                     ULink = Utemp
                     ! Update W
                     CALL MultiplyMatMat(Utemp, embed, W)
                     W = UTemp
                     CALL FixSU3Matrix(ULink)
                     ! And finally update link
                     UUpdated(:, :, mu, it, ix, iy, iz) = ULink
                  END DO
               END DO
            END DO
         END DO
      END DO
   END SUBROUTINE updateLinks

   PURE SUBROUTINE stapleAt(U, V, coord, mu)
      COMPLEX(kind=WC), DIMENSION(:, :, :, :, :, :, :), INTENT(IN) :: U
      INTEGER, DIMENSION(4), INTENT(IN) :: coord
      INTEGER, INTENT(IN) :: mu
      COMPLEX(kind=WC), DIMENSION(3, 3), INTENT(OUT) :: V
      INTEGER, DIMENSION(4) :: thisCoord, step
      INTEGER :: nu
      step = 0
      step(mu) = 1
      thisCoord = coord + step
      thisCoord = periodCoord(thisCoord, SHAPE(U))
      V = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      DO nu = 1, 4
         IF (nu == mu) CYCLE
         V = V + genericPath(U, thisCoord, (/nu, -mu, -nu/)) &
             + genericPath(U, thisCoord, (/-nu, -mu, nu/))
      END DO
   END SUBROUTINE stapleAt

END MODULE FLUE_heatbath
