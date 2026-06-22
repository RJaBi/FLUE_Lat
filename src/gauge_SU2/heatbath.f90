MODULE FLUE_SU2_heatbath
   !!
   !! SU(2) heatbath update using Philox only.
   !!
   !! Key points:
   !!   * No stdlib random routines are used.
   !!   * getLambda2() uses philox_uniform_oo on (0,1),
   !!   * getXVec() uses philox_uniform_co on [-1,1)
   !!   * The update is staged by direction and colour so the inner site loop can be
   !!     expressed as do concurrent.
   !!
   USE FLUE_constants,      ONLY : WP, WC, TWO_PI
   USE FLUE_philox_helpers, ONLY : derive_stage_key, site_linear_index, site_colour
   USE FLUE_SU2_random,     ONLY : constructSU2Matrix
   USE FLUE_SU2_wloops,     ONLY : SU2_genericPath
   USE FLUE_wloops,         ONLY : periodCoord
   USE Philox,     ONLY : philox_uniform_co, philox_uniform_oo, c64
   IMPLICIT NONE(TYPE, EXTERNAL)
   PRIVATE

   PUBLIC :: SU2_updateLinks
   PUBLIC :: constructXMatrix

CONTAINS

  PURE FUNCTION getLambda2(alpha, beta, key, counter0) RESULT(lambda2)
     !!
    !! Sample lambda^2 for the SU(2) heatbath accept/reject step.
    !!
    !! uses Philox draws directly.
    !!
    !! Random-number conventions:
    !!   * ri(1:3) and rr are strictly in (0,1), so we use philox_uniform_oo.
    !!
    !! Counter usage:
    !!   counter(1) = site id        (passed in through counter0)
    !!   counter(2) = subgroup id    (set in constructXMatrix)
    !!   counter(3) = attempt number (incremented here)
    !!   counter(4) = phase id = 1   (lambda2 phase)
    !!
    REAL(kind=WP), INTENT(IN) :: alpha, beta
    INTEGER(C64), INTENT(IN)  :: key(2), counter0(4)
    REAL(kind=WP) :: lambda2, rr
    REAL(kind=WP) :: ri(3)
    INTEGER(C64) :: c(4), attempt
    REAL(kind=WP), PARAMETER :: eps = tiny(1.0_WP)
    IF (alpha <= eps .or. beta <= eps) THEN
       lambda2 = 0.0_WP
       RETURN
    END IF
    attempt = 0_C64
    DO
       c = counter0
       c(3) = attempt
       c(4) = 1_C64
       ! c = [linearised index, subgroup ID, attempt, 1]
       ri(1) = philox_uniform_oo(c, key, 1_C64, 0.0_WP, 1.0_WP)
       ! counter, key, idx, bounds
       ri(2) = philox_uniform_oo(c, key, 2_C64, 0.0_WP, 1.0_WP)
       ri(3) = philox_uniform_oo(c, key, 3_C64, 0.0_WP, 1.0_WP)
       rr    = philox_uniform_oo(c, key, 4_C64, 0.0_WP, 1.0_WP)
       ! Got random numbers, construct lambda2
       lambda2 = -(log(ri(1)) + log(ri(3)) * cos(TWO_PI * ri(2))**2) / &
            (2.0_WP * alpha * beta)
       ! if lambda2 isn't suitable, try again
       IF (lambda2 < 0.0_WP) THEN
          attempt = attempt + 1_C64
          CYCLE
       END IF
       IF (lambda2 > 1.0_WP) THEN
          attempt = attempt + 1_C64
          CYCLE
       END IF
       ! If it is suitable exit
       IF (rr * rr <= 1.0_WP - lambda2) EXIT
       ! If it isn't, go again
       attempt = attempt + 1_C64
    END DO
  END FUNCTION getLambda2


  PURE FUNCTION getXVec(x0, key, counter0) RESULT(xvec)
    !!
    !! Sample a 3-vector uniformly from the unit ball, then rescale it so that
    !! ||xvec|| = sqrt(1 - x0^2).
    !!
    !! Random-number convention:
    !!   * Each component is drawn from [-1,1), so we use philox_uniform_co.
    !!
    !! Counter usage:
    !!   counter(1) = site id
    !!   counter(2) = subgroup id
    !!   counter(3) = attempt number
    !!   counter(4) = phase id = 2   (x-vector phase)
    !!
    REAL(kind=WP), INTENT(IN) :: x0
    INTEGER(C64), INTENT(IN)  :: key(2), counter0(4)
    REAL(kind=WP) :: xvec(3)
    REAL(kind=WP) :: xlen, requiredLen
    INTEGER(C64) :: c(4), attempt
    REAL(kind=WP), PARAMETER :: eps = tiny(1.0_WP)
    requiredLen = max(0.0_WP, 1.0_WP - x0 * x0)
    IF (requiredLen <= eps) THEN
       xvec = 0.0_WP
       RETURN
    END IF
    attempt = 0_C64
    DO
       c = counter0
       c(3) = attempt
       c(4) = 2_C64
       ! c = [linearised site index, subgroup ID, attempt, 2]
       xvec(1) = philox_uniform_co(c, key, 1_C64, -1.0_WP, 1.0_WP)
       xvec(2) = philox_uniform_co(c, key, 2_C64, -1.0_WP, 1.0_WP)
       xvec(3) = philox_uniform_co(c, key, 3_C64, -1.0_WP, 1.0_WP)
       xlen = sum(xvec**2)
       ! If suitable, exit
       IF (xlen <= 1.0_WP .and. xlen > eps) EXIT
       ! else go again
       attempt = attempt + 1_C64
    END DO
    xvec = xvec * (sqrt(requiredLen) / sqrt(xlen))
  END FUNCTION getXVec


   PURE FUNCTION constructXMatrix(alpha, beta, key, counter0, subgroup_id) RESULT(X)
     !!
     !! Construct the SU(2) heatbath matrix X used in the update.
     !!
     !! This wraps together:
     !!   1) sampling lambda^2,
     !!   2) constructing x0 = 1 - 2 lambda^2,
     !!   3) sampling the spatial part x(1:3),
     !!   4) building the corresponding SU(2) matrix.
     !!
     !! subgroup_id is carried in counter(2) so the three SU(2) sub-updates inside
     !! an SU(3) heatbath update use disjoint Philox substreams.
     !!
     REAL(kind=WP), INTENT(IN) :: alpha, beta
     INTEGER(C64), INTENT(IN)  :: key(2), counter0(4)
     INTEGER, INTENT(IN)       :: subgroup_id
     COMPLEX(kind=WC), DIMENSION(2, 2) :: X
     REAL(kind=WP) :: lambda2
     REAL(kind=WP), DIMENSION(0:3) :: xAll
     INTEGER(C64) :: c(4)
     c = counter0
     c(2) = int(subgroup_id, C64)
     ! c = [linearised index, subgroupID, 0, 0]
     lambda2   = getLambda2(alpha, beta, key, c)
     xAll(0)   = 1.0_WP - 2.0_WP * lambda2
     xAll(1:3) = getXVec(xAll(0), key, c)
     X = constructSU2Matrix(xAll)
   END FUNCTION constructXMatrix


   PURE FUNCTION su2_updated_link(U, beta, coord, mu, key, dims) RESULT(Unew)
      !!
      !! Compute the new SU(2) link at one site and one direction.
      !!
      !! Steps:
      !!   1) Compute the staple V.
      !!   2) Compute alpha from det(V).
      !!   3) Build the heatbath matrix X using Philox.
      !!   4) Return the updated link:
      !!        U' = X * (V / alpha)^dagger
      !!
      !! If alpha is too small, the code falls back to a mild update with alpha=1.
      !!
      COMPLEX(kind=WC), DIMENSION(:, :, :, :, :, :, :), INTENT(IN) :: U
      REAL(kind=WP), INTENT(IN) :: beta
      INTEGER, INTENT(IN) :: coord(4), mu, dims(4)
      INTEGER(C64), INTENT(IN) :: key(2)
      COMPLEX(kind=WC), DIMENSION(2, 2) :: Unew

      COMPLEX(kind=WC), DIMENSION(2, 2) :: XMatrix, V, Vdag
      COMPLEX(kind=WC) :: detV
      REAL(kind=WP) :: alpha
      INTEGER(C64) :: counter0(4)
      REAL(kind=WP), PARAMETER :: eps = tiny(1.0_WP)
      ! Calculate the staple
      CALL stapleAt(U, V, coord, mu)
      ! get dterminant
      detV  = V(1,1) * V(2,2) - V(1,2) * V(2,1)
      ! Calculate alpha
      alpha = sqrt(max(real(detV, kind=WP), 0.0_WP))
      ! Convert the site index into a linear index
      ! So that each site has it's own independent random number space
      counter0 = [ site_linear_index(coord, dims), 0_C64, 0_C64, 0_C64 ]
      IF (alpha < eps) THEN
         XMatrix = constructXMatrix(1.0_WP, beta, key, counter0, 1)
         Unew    = XMatrix
      ELSE
         Vdag    = conjg(transpose(V / alpha))
         XMatrix = constructXMatrix(alpha, beta, key, counter0, 1)
         ! i.e. alpha, beta, random stream, linearised index, SU2 subgroup ID (needed for the 3 SU2 updates in SU3)
         Unew    = matmul(XMatrix, Vdag)
      END IF
   END FUNCTION su2_updated_link


   SUBROUTINE SU2_updateLinks(U, beta, UUpdated, master_key, sweep_id)
      !!
      !! Update all SU(2) links on the lattice.
      !!
      !! The update is staged by:
      !!   * direction mu
      !!   * checkerboard parity colour
      !!
      !! so that each do concurrent block updates an independent set of links.
      !!
      !! Inputs:
      !!   U          : current lattice
      !!   beta       : gauge coupling
      !!   master_key : user-provided Philox master key
      !!   sweep_id   : identifies the sweep so keys vary deterministically
      !!
      !! Output:
      !!   UUpdated   : updated lattice
      !!
      COMPLEX(kind=WC), DIMENSION(:, :, :, :, :, :, :), INTENT(IN)    :: U
      REAL(kind=WP),                                      INTENT(IN)    :: beta
      COMPLEX(kind=WC), DIMENSION(:, :, :, :, :, :, :), INTENT(INOUT) :: UUpdated
      INTEGER(kind=C64), DIMENSION(2),                        INTENT(IN)    :: master_key
      INTEGER,                                           INTENT(IN)    :: sweep_id
      INTEGER, DIMENSION(7) :: dataShape
      INTEGER :: nt, nx, ny, nz
      INTEGER :: it, ix, iy, iz, mu, colour
      INTEGER :: dims4(4)
      INTEGER(C64) :: key(2)
      REAL(kind=WP), PARAMETER :: eps = tiny(1.0_WP)
      dataShape = shape(U)
      nt = dataShape(4)
      nx = dataShape(5)
      ny = dataShape(6)
      nz = dataShape(7)
      dims4 = [nt, nx, ny, nz]
      UUpdated = U
      IF (beta <= eps) THEN
         WRITE(*,*) 'beta =', beta, ' must be positive. stopping'
         STOP
      END IF
      ! Loop over directions, checkerboard
      DO mu = 1, 4
         DO colour = 0, 1
            ! Get a key for the random number for this direction, checker, sweep and 'run' (master)
            CALL derive_stage_key(master_key, sweep_id, stage_tag=2, mu=mu, colour=colour, key=key)
            ! Can now parallelise over
            DO CONCURRENT (it = 1:nt, ix = 1:nx, iy = 1:ny, iz = 1:nz) &
                 DEFAULT(NONE) SHARED(UUpdated, beta) LOCAL_INIT(mu, key, dims4, colour)
               IF (site_colour([it, ix, iy, iz], mu, .FALSE.) == colour) THEN
                  ! Skip this execution cause it's on the other 'colour'
                  CYCLE
               END IF
               UUpdated(:, :, mu, it, ix, iy, iz) = su2_updated_link( &
                    UUpdated, beta, [it, ix, iy, iz], mu, key, dims4 )
            END DO
         END DO
      END DO
    END SUBROUTINE SU2_updateLinks


   PURE SUBROUTINE stapleAt(U, V, coord, mu)
      !!
      !! Compute the standard Wilson-type staple for an SU(2) link.
      !!
      !! The staple is evaluated at the link U_mu(coord), and is the sum of the
      !! forward and backward staples in all transverse directions nu != mu.
      !!
      COMPLEX(kind=WC), DIMENSION(:, :, :, :, :, :, :), INTENT(IN) :: U
      INTEGER, DIMENSION(4), INTENT(IN) :: coord
      INTEGER, INTENT(IN) :: mu
      COMPLEX(kind=WC), DIMENSION(2, 2), INTENT(OUT) :: V

      INTEGER, DIMENSION(4) :: thisCoord, step
      INTEGER :: nu
      step = 0
      step(mu) = 1
      thisCoord = coord + step
      thisCoord = periodCoord(thisCoord, shape(U))
      V = cmplx(0.0_WP, 0.0_WP, kind=WC)
      DO nu = 1, 4
         IF (nu == mu) CYCLE
         V = V + SU2_genericPath(U, thisCoord, [ nu, -mu, -nu ]) &
               + SU2_genericPath(U, thisCoord, [ -nu, -mu,  nu ])
      END DO
   END SUBROUTINE stapleAt

END MODULE FLUE_SU2_heatbath
