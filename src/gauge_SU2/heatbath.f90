module FLUE_SU2_heatbath
   !<
   !< SU(2) heatbath update
   !<
   !< Key points:
   !<   * No stdlib random routines are used.
   !<   * getLambda2() uses philox_uniform_oo on (0,1),
   !<   * getXVec() uses philox_uniform_co on [-1,1)
   !<   * The update is staged by direction and colour so the inner site loop can be
   !<     expressed as do concurrent.
   !<
   use FLUE_constants, only: WP, WC, TWO_PI
   use FLUE_philox_helpers, only: derive_stage_key, site_linear_index, site_colour
   use FLUE_SU2_random, only: constructSU2Matrix
   use FLUE_SU2_wloops, only: SU2_genericPath
   use FLUE_wloops, only: periodCoord
   use Philox, only: philox_uniform_co, philox_uniform_oo, c64, philox_fill_uniform_oo, philox_fill_uniform_co
   implicit none(type, external)
   private

   public :: SU2_updateLinks
   public :: constructXMatrix

contains

  pure function getLambda2(alpha, beta, key, counter0) result(lambda2)
    !$omp declare target
     !<
    !< Sample lambda^2 for the SU(2) heatbath accept/reject step.
    !<
    !< uses Philox draws directly.
    !<
    !< Random-number conventions:
    !<   * ri(1:3) and rr are strictly in (0,1), so we use philox_fill_uniform_oo.
    !<
    !< Counter usage:
    !<   counter(1) = site id        (passed in through counter0)
    !<   counter(2) = subgroup id    (set in constructXMatrix)
    !<   counter(3) = attempt number (incremented here)
    !<   counter(4) = phase id = 1   (lambda2 phase)
    !<
     real(kind=WP), intent(IN) :: alpha, beta
     !< alpha and beta parameters (beta is fixed for given sim.)
     integer(kind=C64), intent(IN), dimension(4) :: counter0
     !< Philox keys
      integer(kind=C64), intent(IN), dimension(2) :: key
      !< Philox keys
      real(kind=WP) :: lambda2
      !< Output lambda 2
      real(kind=WP), dimension(4) :: rvec
      !< 4 random numbers
      integer(kind=C64), dimension(4) :: c
      !< The full counter for the philox draw
      integer(kind=C64) :: attempt
      !< A counter/philox key to handle the case where the random number draw does not meet requirements
      real(kind=WP), parameter :: eps = TINY(1.0_WP)
      !< Smallest real number in WP precision
      if (alpha <= eps .OR. beta <= eps) then
         lambda2 = 0.0_WP
         return
      end if
      attempt = 0_C64
      PhiloxAttempt: do
         c = counter0
         c(3) = attempt
         c(4) = 1_C64
         call philox_fill_uniform_oo(c, key, 1_C64, 0.0_WP, 1.0_WP, rvec)
         ! Got random numbers, construct lambda2
         lambda2 = -(LOG(rvec(1)) + LOG(rvec(3)) * COS(TWO_PI * rvec(2))**2) / &
                   (2.0_WP * alpha * beta)
         ! if lambda2 isn't suitable, try again
         if (lambda2 < 0.0_WP) then
            attempt = attempt + 1_C64
            cycle PhiloxAttempt
         end if
         if (lambda2 > 1.0_WP) then
            attempt = attempt + 1_C64
            cycle PhiloxAttempt
         end if
         ! If it is suitable exit
         if (rvec(4) * rvec(4) <= 1.0_WP - lambda2) exit PhiloxAttempt
         ! If it isn't, go again
         attempt = attempt + 1_C64
      end do PhiloxAttempt
   end function getLambda2

   pure function getXVec(x0, key, counter0) result(xvec)
    !<
    !< Sample a 3-vector uniformly from the unit ball, then rescale it so that
    !< ||xvec|| = sqrt(1 - x0^2).
    !<
    !< Random-number convention:
    !<   * Each component is drawn from [-1,1), so we use philox_fill_uniform_co.
    !<
    !< Counter usage:
    !<   counter(1) = site id
    !<   counter(2) = subgroup id
    !<   counter(3) = attempt number
    !<   counter(4) = phase id = 2   (x-vector phase)
    !<
     real(kind=WP), intent(IN) :: x0
     !< The zeroth random number that sets the condition on the 3-vec of random numbers
     integer(kind=C64), intent(IN), dimension(4) :: counter0
     !< Philox keys
      integer(kind=C64), intent(IN), dimension(2) :: key
      !< Philox key
      real(kind=WP), dimension(3) :: xvec
      !< 3-vector of random numbers
      real(kind=WP) :: xlen, requiredLen
      !< Length before rescaling, rescaled length
      integer(kind=C64), dimension(4) :: c
      !< The full counter for the philox draw
      integer(kind=C64) :: attempt
      !< A counter/philox key to handle the case where the random number draw does not meet the requirements
      real(kind=WP), parameter :: eps = TINY(1.0_WP)
      !< Smallest real number in WP precision
      requiredLen = MAX(0.0_WP, 1.0_WP - x0 * x0)
      if (requiredLen <= eps) then
         xvec = 0.0_WP
         return
      end if
      attempt = 0_C64
      PhiloxAttempt: do
         c = counter0
         c(3) = attempt
         c(4) = 2_C64
         ! c = [linearised site index, subgroup ID, attempt, 2]
         call philox_fill_uniform_co(c, key, 1_C64, -1.0_WP, 1.0_WP, xvec)
         !xvec(1) = philox_uniform_co(c, key, 1_C64, -1.0_WP, 1.0_WP)
         !xvec(2) = philox_uniform_co(c, key, 2_C64, -1.0_WP, 1.0_WP)
         !xvec(3) = philox_uniform_co(c, key, 3_C64, -1.0_WP, 1.0_WP)
         xlen = SUM(xvec**2)
         ! If suitable, exit
         if (xlen <= 1.0_WP .AND. xlen > eps) exit PhiloxAttempt
         ! else go again
         attempt = attempt + 1_C64
      end do PhiloxAttempt
      xvec = xvec * (SQRT(requiredLen) / SQRT(xlen))
   end function getXVec

   pure function constructXMatrix(alpha, beta, key, counter0, subgroup_id) result(X)
     !$omp declare target
     !<
     !< Construct the SU(2) heatbath matrix X used in the update.
     !<
     !< This wraps together:
     !<   1) sampling lambda^2,
     !<   2) constructing x0 = 1 - 2 lambda^2,
     !<   3) sampling the spatial part x(1:3),
     !<   4) building the corresponding SU(2) matrix.
     !<
     !< subgroup_id is carried in counter(2) so the three SU(2) sub-updates inside
     !< an SU(3) heatbath update use disjoint Philox substreams.
     !<
     real(kind=WP), intent(IN) :: alpha, beta
     !< alpha and beta parameters (beta fixed for given sim.)
     integer(kind=C64), intent(IN), dimension(4) :: counter0
     !< Philox Keys
      integer(kind=C64), intent(IN), dimension(2) :: key
      !< Philox keys
      integer, intent(IN) :: subgroup_id
      !< Describes the SU2 subgroup
      complex(kind=WC), dimension(2, 2) :: X
      !< SU2 matrix X
      real(kind=WP) :: lambda2
      !< The appropriately distritubted lambda2
      real(kind=WP), dimension(0:3) :: xAll
      !< the 3-vec of random numbers of given length
      integer(kind=C64), dimension(4) :: c
      !< The full counter for the philox draw [linearised index, subgroupID, 0, 0]
      c = counter0
      c(2) = INT(subgroup_id, C64)
      ! c = [linearised index, subgroupID, 0, 0]
      lambda2 = getLambda2(alpha, beta, key, c)
      xAll(0) = 1.0_WP - 2.0_WP * lambda2
      xAll(1:3) = getXVec(xAll(0), key, c)
      X = constructSU2Matrix(xAll)
   end function constructXMatrix

   pure function su2_updated_link(U, beta, coord, mu, key, dims) result(Unew)
      !<
      !< Compute the new SU(2) link at one site and one direction.
      !<
      !< Steps:
      !<   1) Compute the staple V.
      !<   2) Compute alpha from det(V).
      !<   3) Build the heatbath matrix X using Philox.
      !<   4) Return the updated link:
      !<        U' = X * (V / alpha)^dagger
      !<
      !< If alpha is too small, the code falls back to a mild update with alpha=1.
      !<
     complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(IN) :: U
     !< Input gaugefield [2,2,4,nt,nx,ny,nz]
      real(kind=WP), intent(IN) :: beta
      !< beta-value
      integer, intent(IN), dimension(4) :: coord
      !< Base coordinate of link to update [it,ix,iy,zi]
      integer, intent(in) :: mu
      !< Direction to update
      integer, intent(in), dimension(4) :: dims
      !< Lattice dimensions [nt,nx,ny,nz]
      integer(kind=C64), intent(IN), dimension(2) :: key
      !< Philox key
      complex(kind=WC), dimension(2, 2) :: Unew
      !< Updated SU2 link
      complex(kind=WC), dimension(2, 2) :: XMatrix, V, Vdag
      !< Staple, StapleDagger, X-matrix
      complex(kind=WC) :: detV
      !< determinant of staple
      real(kind=WP) :: alpha
      !< sqrt(det(staple))
      integer(kind=C64), dimension(4) :: counter0
      !< Philox counter
      real(kind=WP), parameter :: eps = TINY(1.0_WP)
      !< Smallest real number in WP precision
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
      !<
      !< Update all SU(2) links on the lattice.
      !<
      !< The update is staged by:
      !<   * direction mu
      !<   * checkerboard parity colour
      !<
      !< so that each do concurrent block updates an independent set of links.
      !<
      !< Inputs:
      !<   U          : current lattice
      !<   beta       : gauge coupling
      !<   master_key : user-provided Philox master key
      !<   sweep_id   : identifies the sweep so keys vary deterministically
      !<
      !< Output:
      !<   UUpdated   : updated lattice
      !<
     complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(IN) :: U
     !< Gaugefield [2,2,4,nt,nx,ny,nz]
      real(kind=WP), intent(IN) :: beta
      !< beta-value
      complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(INOUT) :: UUpdated
      !< Output updated gaugefield [2,2,4,nt,nx,ny,nz]
      integer(kind=C64), dimension(2), intent(IN) :: master_key
      !< Philox Key
      integer, intent(IN) :: sweep_id
      !< A counter for Philox
      integer, dimension(7) :: dataShape
      !< Gaugefield diemnsions [2,2,4,nt,nx,ny,nz]
      integer :: nt, nx, ny, nz
      !< Lattice dimensions
      integer :: it, ix, iy, iz, mu, colour
      !< Loop counters
      integer, dimension(4) :: dims4
      !< Lattice dimensions [nt,nx,ny,nz]
      integer(kind=C64), dimension(2) :: key
      !< Philox key
      integer, dimension(4) :: coord
      !< Single site coordinate [it,ix,iy,zi]
      real(kind=WP), parameter :: eps = TINY(1.0_WP)
      !< Smallest real number in WP precision
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
            site: do concurrent(it=1:nt, ix=1:nx, iy=1:ny, iz=1:nz) &
               DEFAULT(none) SHARED(UUpdated, beta) LOCAL_INIT(mu, key, dims4, colour, coord)
               coord = (/it, ix, iy, iz/)
               if (site_colour(coord, mu, .FALSE.) /= colour) then
                  ! Skip this execution cause it's on the other 'colour'
                  cycle site
               end if
               UUpdated(:, :, mu, it, ix, iy, iz) = su2_updated_link( &
                                                    UUpdated, beta, coord, mu, key, dims4)
            end do site
         end do
      end do
   end subroutine SU2_updateLinks

   pure subroutine stapleAt(U, V, coord, mu)
      !<
      !< Compute the standard Wilson-type staple for an SU(2) link.
      !<
      !< The staple is evaluated at the link U_mu(coord), and is the sum of the
      !< forward and backward staples in all transverse directions nu != mu.
      !<
     complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(IN) :: U
     !< Input gaugefield [2,2,4,nt,nx,ny,nz]
      integer, dimension(4), intent(IN) :: coord
      !< Staple is at this coordinate
      integer, intent(IN) :: mu
      !< Staple is in this direction
      complex(kind=WC), dimension(2, 2), intent(OUT) :: V
      !< Output Wilson-type staple
      integer, dimension(4) :: thisCoord
      !< single coordinate [it,ix,iy,iz]
      integer, dimension(4) :: step
      !< Used to navigate around gaugefield
      integer :: nu
      !< counter
      integer, dimension(7) :: dataShape
      !< Shape of gaugefield [2,2,4,it,ix,iy,iz]
      dataShape = SHAPE(U)
      step = 0
      step(mu) = 1
      thisCoord = coord + step
      thisCoord = periodCoord(thisCoord, dataShape)
      V = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      direction: do nu = 1, 4
         if (nu == mu) cycle direction
         V = V + SU2_genericPath(U, thisCoord, [nu, -mu, -nu]) &
             + SU2_genericPath(U, thisCoord, [-nu, -mu, nu])
      end do direction
   end subroutine stapleAt

end module FLUE_SU2_heatbath
