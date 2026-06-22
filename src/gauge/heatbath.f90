MODULE FLUE_heatbath
  !!
  !! SU(3) heatbath update built from three embedded SU(2) subgroup updates.
  !!
  !! This module follows the structure of the original code:
  !!   * compute staple
  !!   * form W = U * staple
  !!   * perform subgroup updates in (1,2), (1,3), (2,3)
  !!   * reunitarize
  !!
  !! The full lattice update is staged by direction mu and colour class, allowing
  !! the site loop inside each stage to be expressed as do concurrent.
  !!
  !! Wilson   -> 2 colours
  !! Symanzik -> 4 colours
  !!
  USE FLUE_constants,       ONLY : WP, WC
  USE FLUE_matrixConstants, ONLY : Ident3x3
  USE FLUE_philox_helpers,  ONLY : derive_stage_key, site_linear_index, &
       site_colour, ncolours_for_action
  USE FLUE_SU2_heatbath,    ONLY : constructXMatrix
  USE FLUE_SU2_random,      ONLY : constructSU2Matrix
  USE FLUE_SU3MatrixOps,    ONLY : FixSU3Matrix, MultiplyMatMat
  USE FLUE_wloops,          ONLY : genericPath, periodCoord
  USE philox, ONLY: c64
  USE stdlib_ascii,         ONLY : to_lower
  IMPLICIT NONE(TYPE, EXTERNAL)
  PRIVATE
  PUBLIC :: updateLinks
  PUBLIC :: build_colour_sites

CONTAINS

  ! The force inline seems to help speed sometimes?
  !dir$ forceinline
  PURE FUNCTION su3_updated_link(U, beta, it, ix, iy, iz, mu, key, dims4, use_symanzik) RESULT(Uout)
    COMPLEX(kind=WC), INTENT(IN) :: U(:,:,:,:,:,:,:)
    REAL(kind=WP), INTENT(IN) :: beta
    INTEGER, INTENT(IN) :: it, ix, iy, iz, mu
    INTEGER, INTENT(IN) :: dims4(4)
    LOGICAL, INTENT(IN) :: use_symanzik
    INTEGER(C64), INTENT(IN) :: key(2)
    COMPLEX(kind=WC) :: Uout(3,3)
    !!
    !! Compute the updated SU(3) link at one site and one direction.
    !!
    !! Steps:
    !!   1) Extract the current link U_mu(coord).
    !!   2) Compute either the Wilson or Symanzik staple.
    !!   3) Form W = U * staple.
    !!   4) Apply the three Cabibbo-Marinari SU(2) subgroup updates:
    !!        (1,2), (1,3), (2,3)
    !!   5) Reunitarize the final matrix.
    !!
    COMPLEX(kind=WC), DIMENSION(3, 3) :: staple, W
    INTEGER(C64) :: counter0(4)
    Uout = U(:, :, mu, it, ix, iy, iz)
    IF (use_symanzik) THEN
       CALL stapleSymanzik(U, staple, (/it, ix, iy, iz/), mu)
    ELSE
       CALL stapleWilson(U, staple, (/it, ix, iy, iz/), mu)
    END IF
    CALL MultiplyMatMat(W, Uout, staple)
    counter0 = [ site_linear_index((/it, ix, iy, iz/), dims4), 0_C64, 0_C64, 0_C64 ]
    CALL apply_su2_subgroup_update(Uout, W, 1, 2, 1, beta, key, counter0)
    ! current link, product * staple being updated in line with ULink,
    ! i1, i2 to identify which SU3 subgroup,
    ! subgroup ID used for Philox substream
    ! beta, Philox key
    ! counter with linearised site index
    CALL apply_su2_subgroup_update(Uout, W, 1, 3, 2, beta, key, counter0)
    CALL apply_su2_subgroup_update(Uout, W, 2, 3, 3, beta, key, counter0)
    CALL FixSU3Matrix(Uout)
   END FUNCTION su3_updated_link

   PURE SUBROUTINE apply_su2_subgroup_update(ULink, W, i1, i2, subgroup_id, beta, key, counter0)
     !!
     !! Apply one embedded SU(2) heatbath update inside SU(3).
     !!
     !! Inputs:
     !!   ULink       : current SU(3) link being updated
     !!   W           : current product U * staple, updated in tandem with ULink
     !!   i1, i2      : the two SU(3) row/column indices defining the subgroup
     !!   subgroup_id : 1, 2, or 3, used for the Philox substream
     !!   beta, key   : heatbath parameters and Philox key
     !!   counter0    : counter with site id in counter(1)
     !!
     !! Method:
     !!   1) Extract the 2x2 block from W.
     !!   2) Project it to the SU(2) quaternion representation.
     !!   3) Build the SU(2) heatbath matrix.
     !!   4) Embed the result back into SU(3).
     !!   5) Left-multiply both ULink and W by the embedded matrix.
     !!
     COMPLEX(kind=WC), DIMENSION(3, 3), INTENT(INOUT) :: ULink, W
     INTEGER, INTENT(IN) :: i1, i2, subgroup_id
     REAL(kind=WP), INTENT(IN) :: beta
     INTEGER(C64), INTENT(IN) :: key(2), counter0(4)
     COMPLEX(kind=WC), DIMENSION(2, 2) :: MfromW, SU2_M, SU2_X, SU2_U
     COMPLEX(kind=WC), DIMENSION(3, 3) :: embed, Wtemp
     REAL(kind=WP), DIMENSION(0:3) :: aQuart
     REAL(kind=WP) :: alpha
     REAL(kind=WP), PARAMETER :: eps = tiny(1.0_WP)
     MfromW(1,1) = W(i1, i1)
     MfromW(2,2) = W(i2, i2)
     MfromW(1,2) = W(i1, i2)
     MfromW(2,1) = W(i2, i1)
     aQuart(0) = real(MfromW(1,1) + MfromW(2,2), kind=WP)
     aQuart(1) = aimag(MfromW(1,2) + MfromW(2,1))
     aQuart(2) = real(MfromW(1,2) - MfromW(2,1), kind=WP)
     aQuart(3) = aimag(MfromW(1,1) - MfromW(2,2))
     aQuart    = 0.5_WP * aQuart
     alpha = sqrt(sum(aQuart**2))
     IF (alpha <= eps) THEN
        SU2_U = constructXMatrix(1.0_WP, 2.0_WP * beta / 3.0_WP, key, counter0, subgroup_id)
     ELSE
        aQuart = aQuart / alpha
        SU2_M  = constructSU2Matrix(aQuart)
        SU2_X  = constructXMatrix(alpha, 2.0_WP * beta / 3.0_WP, key, counter0, subgroup_id)
        SU2_M  = conjg(transpose(SU2_M))
        SU2_U  = matmul(SU2_X, SU2_M)
     END IF
     embed = Ident3x3
     embed(i1, i1) = SU2_U(1,1)
     embed(i2, i2) = SU2_U(2,2)
     embed(i1, i2) = SU2_U(1,2)
     embed(i2, i1) = SU2_U(2,1)
     CALL MultiplyMatMat(ULink, embed, ULink)
     WTemp = W
     CALL MultiplyMatMat(W, embed, WTemp)
   END SUBROUTINE apply_su2_subgroup_update


   SUBROUTINE build_colour_sites(nt, nx, ny, nz, use_symanzik, &
        sites_it, sites_ix, sites_iy, sites_iz, counts)
     !! Separates sites out into checkerboard
     !! Separate arrays for each of t, x, y, z
     INTEGER, INTENT(IN) :: nt, nx, ny, nz
     LOGICAL, INTENT(IN) :: use_symanzik
     INTEGER, ALLOCATABLE, INTENT(OUT) :: sites_it(:,:,:)
     INTEGER, ALLOCATABLE, INTENT(OUT) :: sites_ix(:,:,:)
     INTEGER, ALLOCATABLE, INTENT(OUT) :: sites_iy(:,:,:)
     INTEGER, ALLOCATABLE, INTENT(OUT) :: sites_iz(:,:,:)
     INTEGER, ALLOCATABLE, INTENT(OUT) :: counts(:,:)
     INTEGER :: mu, colour, it, ix, iy, iz
     INTEGER :: idx, ncolours, max_sites
     ncolours = ncolours_for_action(use_symanzik)
     max_sites = nt*nx*ny*nz
     ALLOCATE(sites_it(max_sites, ncolours, 4))
     ALLOCATE(sites_ix(max_sites, ncolours, 4))
     ALLOCATE(sites_iy(max_sites, ncolours, 4))
     ALLOCATE(sites_iz(max_sites, ncolours, 4))
     ALLOCATE(counts(ncolours,4))
     counts = 0
     DO mu = 1, 4
        DO it = 1, nt
           DO ix = 1, nx
              DO iy = 1, ny
                 DO iz = 1, nz
                    colour = site_colour([it,ix,iy,iz], mu, use_symanzik)
                    counts(colour+1, mu) = counts(colour+1, mu) + 1
                    idx = counts(colour+1, mu)
                    sites_it(idx,colour+1,mu) = it
                    sites_ix(idx,colour+1,mu) = ix
                    sites_iy(idx,colour+1,mu) = iy
                    sites_iz(idx,colour+1,mu) = iz
                 END DO
              END DO
           END DO
        END DO
     END DO
   END SUBROUTINE build_colour_sites

   SUBROUTINE updateLinks(U, beta, UUpdated, master_key, sweep_id,sites_it, sites_ix, sites_iy, sites_iz, counts, actionTag)
     !------------------------------------
     ! Arguments
     !------------------------------------
     COMPLEX(kind=WC), INTENT(IN)    :: U(:,:,:,:,:,:,:)
     REAL(kind=WP),    INTENT(IN)    :: beta
     COMPLEX(kind=WC), INTENT(INOUT) :: UUpdated(:,:,:,:,:,:,:)
     INTEGER(C64), INTENT(IN) :: master_key(2)
     INTEGER, INTENT(IN)      :: sweep_id
     INTEGER, INTENT(IN) :: sites_it(:,:,:)
     INTEGER, INTENT(IN) :: sites_ix(:,:,:)
     INTEGER, INTENT(IN) :: sites_iy(:,:,:)
     INTEGER, INTENT(IN) :: sites_iz(:,:,:)
     INTEGER, INTENT(IN) :: counts(:,:)
     CHARACTER(len=*), INTENT(IN), OPTIONAL :: actionTag
     !------------------------------------
     ! Locals
     !------------------------------------
     INTEGER :: mu, colour, k
     INTEGER :: it, ix, iy, iz
     INTEGER :: ncolours
     INTEGER :: dims4(4)
     LOGICAL :: use_symanzik
     INTEGER(C64) :: key(2)
     COMPLEX(kind=WC) :: Uloc(3,3)
     INTEGER :: shapeU(7)
     !------------------------------------
     ! Setup
     !------------------------------------
     shapeU = shape(U)
     dims4 = shapeU(4:7)
     use_symanzik = .FALSE.
     IF (present(actionTag)) THEN
        SELECT CASE (to_lower(trim(actionTag)))
        CASE ('symanzik')
           use_symanzik = .TRUE.
        CASE ('wilson')
           use_symanzik = .FALSE.
        END SELECT
     END IF
     ncolours = ncolours_for_action(use_symanzik)
     UUpdated = U
     !====================================
     ! MAIN UPDATE
     !====================================
     DO mu = 1, 4
        DO colour = 0, ncolours - 1
           CALL derive_stage_key(master_key, sweep_id, &
                stage_tag=3, mu=mu, colour=colour, key=key)
           DO CONCURRENT (k = 1:counts(colour+1, mu)) LOCAL(it,ix,iy,iz,Uloc)
              !------------------------------------
              ! Unit-stride coordinate loads (SoA)
              !------------------------------------
              it = sites_it(k,colour+1,mu)
              ix = sites_ix(k,colour+1,mu)
              iy = sites_iy(k,colour+1,mu)
              iz = sites_iz(k,colour+1,mu)
              !------------------------------------
              ! Compute update (fully scalar args)
              !------------------------------------
              Uloc = su3_updated_link( UUpdated, beta, it, ix, iy, iz, mu, key, dims4, use_symanzik )
              !------------------------------------
              ! Store result
              !------------------------------------
              UUpdated(:, :, mu, it, ix, iy, iz) = Uloc
           END DO
        END DO
     END DO
   END SUBROUTINE updateLinks

   PURE SUBROUTINE stapleWilson(U, V, coord, mu)
     !!
     !! Compute the Wilson staple for one SU(3) link.
     !!
     !!   V = sum_{nu != mu} [ path(nu,-mu,-nu) + path(-nu,-mu,nu) ]
     !! with the starting coordinate shifted to x + mu.
     !!
     COMPLEX(kind=WC), DIMENSION(:, :, :, :, :, :, :), INTENT(IN) :: U
     INTEGER, DIMENSION(4), INTENT(IN) :: coord
     INTEGER, INTENT(IN) :: mu
     COMPLEX(kind=WC), DIMENSION(3, 3), INTENT(OUT) :: V
     INTEGER, DIMENSION(4) :: thisCoord, step
     INTEGER :: nu
     step = 0
     step(mu) = 1
     thisCoord = coord + step
     thisCoord = periodCoord(thisCoord, shape(U))
     V = cmplx(0.0_WP, 0.0_WP, kind=WC)
     DO nu = 1, 4
        IF (nu == mu) CYCLE
        V = V + genericPath(U, thisCoord, [  nu, -mu, -nu ]) &
             + genericPath(U, thisCoord, [ -nu, -mu,  nu ])
     END DO
   END SUBROUTINE stapleWilson

   PURE SUBROUTINE stapleSymanzik(U, V, coord, mu)
     !!
     !! Compute the tree-level Symanzik staple for one SU(3) link.
     !!
     !! This consists of:
     !!   * the plaquette staple contribution (Wilson staple),
     !!   * minus the weighted rectangle contribution.
     !!
     COMPLEX(kind=WC), DIMENSION(:, :, :, :, :, :, :), INTENT(IN) :: U
     INTEGER, DIMENSION(4), INTENT(IN) :: coord
     INTEGER, INTENT(IN) :: mu
     COMPLEX(kind=WC), DIMENSION(3, 3), INTENT(OUT) :: V
     INTEGER, DIMENSION(4) :: thisCoord, step
     INTEGER, DIMENSION(5) :: r5
     COMPLEX(kind=WC), DIMENSION(3, 3) :: Vplaq, Vrect
     INTEGER :: nu
     CALL stapleWilson(U, Vplaq, coord, mu)
     Vrect = cmplx(0.0_WP, 0.0_WP, kind=WC)
     step = 0
     step(mu) = 1
     thisCoord = periodCoord(coord + step, shape(U))  ! x + mu
     DO nu = 1, 4
        IF (nu == mu) CYCLE
        ! Long in nu
        r5 = [  nu,  nu, -mu, -nu, -nu ]
        Vrect = Vrect + genericPath(U, thisCoord, r5)
        r5 = [ -nu, -nu, -mu,  nu,  nu ]
        Vrect = Vrect + genericPath(U, thisCoord, r5)
        ! Long in mu, link is the first mu
        r5 = [  mu,  nu, -mu, -mu, -nu ]
        Vrect = Vrect + genericPath(U, thisCoord, r5)
        r5 = [  mu, -nu, -mu, -mu,  nu ]
        Vrect = Vrect + genericPath(U, thisCoord, r5)
        ! Long in mu, link is the second mu
        r5 = [  nu, -mu, -mu, -nu,  mu ]
        Vrect = Vrect + genericPath(U, thisCoord, r5)
        r5 = [ -nu, -mu, -mu,  nu,  mu ]
        Vrect = Vrect + genericPath(U, thisCoord, r5)
     END DO
     V = (5.0_WP / 3.0_WP) * Vplaq - (1.0_WP / 12.0_WP) * Vrect
   END SUBROUTINE stapleSymanzik

END MODULE FLUE_heatbath
