MODULE FLUE_SU2_heatbath
   USE FLUE_constants, ONLY: WP, WC, TWO_PI
   USE FLUE_matrixConstants, ONLY: Ident2x2
   USE FLUE_SU2_random, ONLY: constructSU2Matrix
   USE FLUE_SU2_wloops, ONLY: SU2_genericPath
   USE FLUE_wloops, ONLY: periodCoord
   USE stdlib_stats_distribution_uniform, ONLY: rvs_uniform
   IMPLICIT NONE(TYPE, EXTERNAL)
   PRIVATE
   PUBLIC :: SU2_updateLinks
   PUBLIC :: constructXMatrix
CONTAINS
   FUNCTION getLambda2(alpha, beta) RESULT(lambda2)
      REAL(kind=WP), INTENT(IN) :: alpha, beta
      REAL(kind=WP) :: lambda2, rr
      REAL(kind=WP) :: ri(3)
      REAL(kind=WP), PARAMETER :: eps = TINY(1.0_WP)
      ! Defensive: beta or alpha too small -> distribution becomes Haar-like.
      ! Caller should handle beta~0 separately, but we also guard here.
      IF (alpha <= eps .OR. beta <= eps) THEN
         lambda2 = 0.0_WP
         RETURN
      END IF
      DO
         ri = rvs_uniform(loc=0.0_WP, scale=1.0_WP, array_size=3)
         IF (ANY(ri <= eps)) CYCLE  ! avoid log(0)
         ! Eq. (4.45): proposal for lambda^2
         lambda2 = -(LOG(ri(1)) + LOG(ri(3)) * COS(TWO_PI * ri(2))**2) &
                   / (2.0_WP * alpha * beta)
         ! Make sure inside bounds
         IF (lambda2 < 0.0_WP) CYCLE
         IF (lambda2 > 1.0_WP) CYCLE
         ! Eq. (4.46): accept with prob sqrt(1 - lambda^2)
         rr = rvs_uniform(scale=1.0_WP)
         IF (rr <= eps) CYCLE
         ! done?
         IF (rr * rr <= 1.0_WP - lambda2) EXIT
      END DO
   END FUNCTION getLambda2

   FUNCTION getXVec(x0) RESULT(xvec)
      REAL(kind=WP), INTENT(IN) :: x0
      REAL(kind=WP) :: xvec(3)
      REAL(kind=WP) :: xlen, requiredLen
      REAL(kind=WP), PARAMETER :: eps = TINY(1.0_WP)
      requiredLen = MAX(0.0_WP, 1.0_WP - x0 * x0)          ! requiredLen = 1 - x0^2
      IF (requiredLen <= eps) THEN
         xvec = 0.0_WP
         RETURN
      END IF
      DO
         xvec = rvs_uniform(loc=-1.0_WP, scale=2.0_WP, array_size=3)
         ! Reject exactly +1 only if you insist on [-1,1) — but using eps is safer
         IF (ANY(xvec >= 1.0_WP - eps)) CYCLE
         ! check conditions
         xlen = SUM(xvec**2.0_WP)
         IF (xlen <= 1.0_WP .AND. xlen > eps) EXIT
      END DO
      ! Rescale to |xvec| = sqrt(1 - x0^2)
      xvec = xvec * (SQRT(requiredLen) / SQRT(xlen))
   END FUNCTION getXVec

   FUNCTION constructXMatrix(alpha, beta) RESULT(X)
    !! section 4.3.1
      REAL(kind=WP), INTENT(IN) :: alpha, beta
      COMPLEX(kind=WC), DIMENSION(2, 2) :: X
      REAL(kind=WP) :: lambda2
      REAL(kind=WP), DIMENSION(0:3) :: xAll
    !! Get lambda
      lambda2 = getLambda2(alpha, beta)
      xAll(0) = 1.0_WP - 2.0_WP * lambda2
      xAll(1:3) = getXVec(xAll(0))
      ! debug: check unit quaternion
      !if (abs( xAll(0)**2 + sum(xAll(1:3)**2) - 1.0_WP ) > 1e-10_WP) then
      !   write(*,*) 'xvec', xAll
      !   stop
      !end if
      X = constructSU2Matrix(xAll)
   END FUNCTION constructXMatrix

   SUBROUTINE SU2_updateLinks(U, beta, UUpdated)
      COMPLEX(kind=WC), DIMENSION(:, :, :, :, :, :, :), INTENT(IN) :: U
      REAL(kind=WP), INTENT(IN) :: beta
      COMPLEX(kind=WC), DIMENSION(:, :, :, :, :, :, :), INTENT(INOUT) :: UUpdated
      !complex(kind=WC), allocatable, dimension(:,:,:,:,:,:,:) :: staples
      COMPLEX(kind=WC), DIMENSION(2, 2) :: XMatrix, V, Vdag
      COMPLEX(kind=WC) :: detV
      REAL(kind=WP) :: alpha
    !! lattice geometry
      INTEGER, DIMENSION(7) :: dataShape
      INTEGER :: nt, nx, ny, nz
    !! counters
      INTEGER :: ix, iy, iz, it, mu
      INTEGER, DIMENSION(4) :: coord
      REAL(kind=WP), PARAMETER :: eps = TINY(1.0_WP)
      dataShape = SHAPE(U)
      nt = dataShape(4)
      nx = dataShape(5)
      ny = dataShape(6)
      nz = dataShape(7)
      !allocate(staples(2, 2, 4, nt, nx, ny, nz))
      !call computeStaples(U, staples)
      UUpdated = U
      IF (beta .LE. eps) THEN
         WRITE (*, *) 'beta', beta, ' must not be zero. stopping'
         STOP
      END IF
      DO it = 1, nt
         DO ix = 1, nx
            DO iy = 1, ny
               DO iz = 1, nz
                  coord = (/it, ix, iy, iz/)
                  DO mu = 1, 4
                     CALL stapleAt(UUpdated, V, coord, mu)
                     ! ad - bc
                     detV = V(1, 1) * V(2, 2) - V(1, 2) * V(2, 1)
                     ! For SU(2) staples, det should be real and positive (up to roundoff)
                     !if (abs(aimag(detV)) > 1.0e-8_WP) then
                     !   write(*,*) 'detV', detV
                     !   write(*,*) 'V', V
                     !   !error stop "Staple det has significant imaginary part"
                     !end if
                     ! calculate alpha from the det
                     alpha = SQRT(MAX(real(detV, kind=WP), 0.0_WP))
                     IF (alpha < eps) THEN
                        ! Degenerate staple: fall back to Haar-ish update (or identity)
                        XMatrix = constructXMatrix(1.0_WP, beta)    ! mild fallback
                        UUpdated(:, :, mu, it, ix, iy, iz) = XMatrix
                     ELSE
                        ! Vtilde = V / alpha, and we need Vtilde^dagger on the right
                        Vdag = CONJG(TRANSPOSE(V / alpha))
                        !Vdag = V / alpha
                        XMatrix = constructXMatrix(alpha, beta)
                        ! Update: U' = X * Vtilde^dagger  (Montvay–Münster / Gattringer–Lang)
                        UUpdated(:, :, mu, it, ix, iy, iz) = MATMUL(XMatrix, Vdag)
                     END IF
                     !if (it==1 .and. ix==1 .and. iy==1 .and. iz==1 .and. mu==1) then
                     !   write(*,*) "Example alpha=", alpha, "detV=", detV
                     !   ! write(*,*) 'V / alpha', V / alpha
                     !end if
                     !VDag = matmul( conjg(transpose(UUpdated(:,:,mu,it,ix,iy,iz))), UUpdated(:,:,mu,it,ix,iy,iz))
                     !if (maxval(abs(VDag - Ident2x2)) > 1.0e-8_WP) then
                     !   write(*,*) 'Unitarity check'
                     !   write(*,*) 'link', it, ix, iy, iz, mu, UUpdated(:,:,mu,it,ix,iy,iz), 'not unitary', VDag
                     !   stop
                     !end if
                  END DO
               END DO
            END DO
         END DO
      END DO

   END SUBROUTINE SU2_updateLinks

   PURE SUBROUTINE stapleAt(U, V, coord, mu)
      COMPLEX(kind=WC), DIMENSION(:, :, :, :, :, :, :), INTENT(IN) :: U
      INTEGER, DIMENSION(4), INTENT(IN) :: coord
      INTEGER, INTENT(IN) :: mu
      COMPLEX(kind=WC), DIMENSION(2, 2), INTENT(OUT) :: V
      INTEGER, DIMENSION(4) :: thisCoord, step
      INTEGER :: nu
      step = 0
      step(mu) = 1
      thisCoord = coord + step
      thisCoord = periodCoord(thisCoord, SHAPE(U))
      V = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      DO nu = 1, 4
         IF (nu == mu) CYCLE
         V = V + SU2_genericPath(U, thisCoord, (/nu, -mu, -nu/)) &
             + SU2_genericPath(U, thisCoord, (/-nu, -mu, nu/))
      END DO
   END SUBROUTINE stapleAt

END MODULE FLUE_SU2_heatbath
