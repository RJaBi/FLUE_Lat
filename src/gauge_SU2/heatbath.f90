module FLUE_SU2_heatbath
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
   use FLUE_constants, only: WP, WC, TWO_PI
   use FLUE_philox_helpers, only: derive_stage_key, site_linear_index, site_colour
   use FLUE_SU2_random, only: constructSU2Matrix
   use FLUE_SU2_wloops, only: SU2_genericPath
   use FLUE_wloops, only: periodCoord
   use Philox, only: philox_uniform_co, philox_uniform_oo, c64
   implicit none(type, external)
   private

   public :: SU2_updateLinks
   public :: constructXMatrix

contains

   pure function getLambda2(alpha, beta, key, counter0) result(lambda2)
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
      real(kind=WP), intent(IN) :: alpha, beta
      integer(C64), intent(IN) :: key(2), counter0(4)
      real(kind=WP) :: lambda2, rr
      real(kind=WP) :: ri(3)
      integer(C64) :: c(4), attempt
      real(kind=WP), parameter :: eps = TINY(1.0_WP)
      if (alpha <= eps .OR. beta <= eps) then
         lambda2 = 0.0_WP
         return
      end if
      attempt = 0_C64
      do
         c = counter0
         c(3) = attempt
         c(4) = 1_C64
         ! c = [linearised index, subgroup ID, attempt, 1]
         ri(1) = philox_uniform_oo(c, key, 1_C64, 0.0_WP, 1.0_WP)
         ! counter, key, idx, bounds
         ri(2) = philox_uniform_oo(c, key, 2_C64, 0.0_WP, 1.0_WP)
         ri(3) = philox_uniform_oo(c, key, 3_C64, 0.0_WP, 1.0_WP)
         rr = philox_uniform_oo(c, key, 4_C64, 0.0_WP, 1.0_WP)
         ! Got random numbers, construct lambda2
         lambda2 = -(LOG(ri(1)) + LOG(ri(3)) * COS(TWO_PI * ri(2))**2) / &
                   (2.0_WP * alpha * beta)
         ! if lambda2 isn't suitable, try again
         if (lambda2 < 0.0_WP) then
            attempt = attempt + 1_C64
            cycle
         end if
         if (lambda2 > 1.0_WP) then
            attempt = attempt + 1_C64
            cycle
         end if
         ! If it is suitable exit
         if (rr * rr <= 1.0_WP - lambda2) exit
         ! If it isn't, go again
         attempt = attempt + 1_C64
      end do
   end function getLambda2

   pure function getXVec(x0, key, counter0) result(xvec)
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
      real(kind=WP), intent(IN) :: x0
      integer(C64), intent(IN) :: key(2), counter0(4)
      real(kind=WP) :: xvec(3)
      real(kind=WP) :: xlen, requiredLen
      integer(C64) :: c(4), attempt
      real(kind=WP), parameter :: eps = TINY(1.0_WP)
      requiredLen = MAX(0.0_WP, 1.0_WP - x0 * x0)
      if (requiredLen <= eps) then
         xvec = 0.0_WP
         return
      end if
      attempt = 0_C64
      do
         c = counter0
         c(3) = attempt
         c(4) = 2_C64
         ! c = [linearised site index, subgroup ID, attempt, 2]
         xvec(1) = philox_uniform_co(c, key, 1_C64, -1.0_WP, 1.0_WP)
         xvec(2) = philox_uniform_co(c, key, 2_C64, -1.0_WP, 1.0_WP)
         xvec(3) = philox_uniform_co(c, key, 3_C64, -1.0_WP, 1.0_WP)
         xlen = SUM(xvec**2)
         ! If suitable, exit
         if (xlen <= 1.0_WP .AND. xlen > eps) exit
         ! else go again
         attempt = attempt + 1_C64
      end do
      xvec = xvec * (SQRT(requiredLen) / SQRT(xlen))
   end function getXVec

   pure function constructXMatrix(alpha, beta, key, counter0, subgroup_id) result(X)
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
      real(kind=WP), intent(IN) :: alpha, beta
      integer(C64), intent(IN) :: key(2), counter0(4)
      integer, intent(IN) :: subgroup_id
      complex(kind=WC), dimension(2, 2) :: X
      real(kind=WP) :: lambda2
      real(kind=WP), dimension(0:3) :: xAll
      integer(C64) :: c(4)
      c = counter0
      c(2) = INT(subgroup_id, C64)
      ! c = [linearised index, subgroupID, 0, 0]
      lambda2 = getLambda2(alpha, beta, key, c)
      xAll(0) = 1.0_WP - 2.0_WP * lambda2
      xAll(1:3) = getXVec(xAll(0), key, c)
      X = constructSU2Matrix(xAll)
   end function constructXMatrix

   pure function su2_updated_link(U, beta, coord, mu, key, dims) result(Unew)
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
      complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(IN) :: U
      real(kind=WP), intent(IN) :: beta
      integer, intent(IN) :: coord(4), mu, dims(4)
      integer(C64), intent(IN) :: key(2)
      complex(kind=WC), dimension(2, 2) :: Unew

      complex(kind=WC), dimension(2, 2) :: XMatrix, V, Vdag
      complex(kind=WC) :: detV
      real(kind=WP) :: alpha
      integer(C64) :: counter0(4)
      real(kind=WP), parameter :: eps = TINY(1.0_WP)
      ! Calculate the staple
      call stapleAt(U, V, coord, mu)
      ! get dterminant
      detV = V(1, 1) * V(2, 2) - V(1, 2) * V(2, 1)
      ! Calculate alpha
      alpha = SQRT(MAX(real(detV, kind=WP), 0.0_WP))
      ! Convert the site index into a linear index
      ! So that each site has it's own independent random number space
      counter0 = [site_linear_index(coord, dims), 0_C64, 0_C64, 0_C64]
      if (alpha < eps) then
         XMatrix = constructXMatrix(1.0_WP, beta, key, counter0, 1)
         Unew = XMatrix
      else
         Vdag = CONJG(TRANSPOSE(V / alpha))
         XMatrix = constructXMatrix(alpha, beta, key, counter0, 1)
         ! i.e. alpha, beta, random stream, linearised index, SU2 subgroup ID (needed for the 3 SU2 updates in SU3)
         Unew = MATMUL(XMatrix, Vdag)
      end if
   end function su2_updated_link

   subroutine SU2_updateLinks(U, beta, UUpdated, master_key, sweep_id)
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
      complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(IN) :: U
      real(kind=WP), intent(IN) :: beta
      complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(INOUT) :: UUpdated
      integer(kind=C64), dimension(2), intent(IN) :: master_key
      integer, intent(IN) :: sweep_id
      integer, dimension(7) :: dataShape
      integer :: nt, nx, ny, nz
      integer :: it, ix, iy, iz, mu, colour
      integer :: dims4(4)
      integer(C64) :: key(2)
      integer, dimension(4) :: coord
      real(kind=WP), parameter :: eps = TINY(1.0_WP)
      dataShape = SHAPE(U)
      nt = dataShape(4)
      nx = dataShape(5)
      ny = dataShape(6)
      nz = dataShape(7)
      dims4 = [nt, nx, ny, nz]
      UUpdated = U
      if (beta <= eps) then
         write (*, *) 'beta =', beta, ' must be positive. stopping'
         stop
      end if
      ! Loop over directions, checkerboard
      do mu = 1, 4
         do colour = 0, 1
            ! Get a key for the random number for this direction, checker, sweep and 'run' (master)
            call derive_stage_key(master_key, sweep_id, stage_tag=2, mu=mu, colour=colour, key=key)
            ! Can now parallelise over
            do concurrent(it=1:nt, ix=1:nx, iy=1:ny, iz=1:nz) &
               DEFAULT(none) SHARED(UUpdated, beta) LOCAL_INIT(mu, key, dims4, colour, coord)
               coord = (/it, ix, iy, iz/)
               if (site_colour(coord, mu, .FALSE.) /= colour) then
                  ! Skip this execution cause it's on the other 'colour'
                  cycle
               end if
               UUpdated(:, :, mu, it, ix, iy, iz) = su2_updated_link( &
                                                    UUpdated, beta, coord, mu, key, dims4)
            end do
         end do
      end do
   end subroutine SU2_updateLinks

   pure subroutine stapleAt(U, V, coord, mu)
      !!
      !! Compute the standard Wilson-type staple for an SU(2) link.
      !!
      !! The staple is evaluated at the link U_mu(coord), and is the sum of the
      !! forward and backward staples in all transverse directions nu != mu.
      !!
      complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(IN) :: U
      integer, dimension(4), intent(IN) :: coord
      integer, intent(IN) :: mu
      complex(kind=WC), dimension(2, 2), intent(OUT) :: V

      integer, dimension(4) :: thisCoord, step
      integer :: nu
      integer, dimension(7) :: dataShape
      dataShape = SHAPE(U)
      step = 0
      step(mu) = 1
      thisCoord = coord + step
      thisCoord = periodCoord(thisCoord, dataShape)
      V = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      do nu = 1, 4
         if (nu == mu) cycle
         V = V + SU2_genericPath(U, thisCoord, [nu, -mu, -nu]) &
             + SU2_genericPath(U, thisCoord, [-nu, -mu, nu])
      end do
   end subroutine stapleAt

end module FLUE_SU2_heatbath
