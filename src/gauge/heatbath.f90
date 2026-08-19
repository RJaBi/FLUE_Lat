module FLUE_heatbath
  !< SU(3) heatbath update using Cabibbo–Marinari SU(2) subgroups.
  !<
  !< This module implements a lattice gauge heatbath update using
  !< three embedded SU(2) subgroup updates.
  !<
  !< Algorithm:
  !< - Compute staple
  !< - Form W = U * staple
  !< - Apply SU(2) updates in (1,2), (1,3), (2,3)
  !< - Reunitarize
  !<
  !< Update strategy:
  !< - Loop over direction mu
  !< - Loop over colour classes
  !< - Parallelise with do concurrent
  !<
  !< Supported actions:
  !< - Wilson (2 colours)
  !< - Symanzik (4 colours)
  !< - Iwasaki (4 colours)
   use FLUE_constants, only: WP, WC
   use FLUE_matrixConstants, only: Ident3x3
   use FLUE_philox_helpers, only: derive_stage_key, site_linear_index, &
                                  site_colour, ncolours_for_action
   use FLUE_SU2_heatbath, only: constructXMatrix
   use FLUE_SU2_random, only: constructSU2Matrix
   use FLUE_SU3MatrixOps, only: FixSU3Matrix, MultiplyMatMat
   use FLUE_wloops, only: genericPath, periodCoord
   use philox, only: c64
   use stdlib_ascii, only: to_lower
   implicit none(type, external)
   private
   public :: updateLinks
   public :: build_colour_sites
   abstract interface
      !< for staple computation kernels
      !<
      !< Implementations must compute the staple `V` for a given link.
      !<
      !< Used to switch between Wilson / Symanzik / Iwasaki actions.
      pure subroutine stapleInterface(U, V, coord, mu, xi)
        !< StapleInterface
        import :: WP, WC
        implicit none(type, external)
        !$omp declare target
         complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(IN) :: U
         !< Input gaugefield [3,3,4,nt,nx,ny,nz]
         integer, dimension(4), intent(IN) :: coord
         !< Coordinate of the staple [it,ix,iy,iz]
         integer, intent(IN) :: mu
         !< Staple is along mu direction
         real(kind=WP), intent(in) :: xi
         !< The auge Anisotropy
         complex(kind=WC), dimension(3, 3), intent(OUT) :: V
         !< Output staple at this coordinate
      end subroutine stapleInterface
   end interface

contains

  pure function su3_updated_link(U, beta, coord, mu, key, dims4, stapleKernel, xi) result(Uout)
    !$omp declare target
    !<
    !< Compute the updated SU(3) link at one site and one direction.
    !<
    !< Steps:
    !<   1) Extract the current link U_mu(coord).
    !<   2) Compute either the Wilson or Symanzik staple.
    !<   3) Form W = U * staple.
    !<   4) Apply the three Cabibbo-Marinari SU(2) subgroup updates:
    !<        (1,2), (1,3), (2,3)
    !<   5) Reunitarize the final matrix.
    !<
    complex(kind=WC), intent(IN), dimension(:, :, :, :, :, :, :) :: U
    !< Input gaugefield [3,3,4,nt,nx,ny,nz]
      real(kind=WP), intent(IN) :: beta
      !< beta value
      integer, intent(IN) :: mu
      !< Direction of link
      integer, dimension(4), intent(IN) :: dims4
      !< Lattice dimensions [nt,nx,ny,nz]
      integer, dimension(4), intent(in) :: coord
      !< Coordinate of link to be updated [it,ix,iy,iz]
      real(kind=WP), intent(in) :: xi
      !< The gauge anisotropy
      procedure(stapleInterface), pointer, intent(in) :: stapleKernel
      !< Pointer used to switch which action is used
      integer(kind=C64), dimension(2), intent(IN) :: key
      !< Philox variable
      complex(kind=WC), dimension(3, 3) :: Uout
      !< Updated SU(3) link
      complex(kind=WC), dimension(3, 3) :: staple, W
      !< The staple and W= U * staple
      integer(kind=C64), dimension(4) :: counter0
      !< Philox variable
      Uout = U(:, :, mu, coord(1), coord(2), coord(3), coord(4))

      call stapleKernel(U, staple, coord, mu, xi)
      call MultiplyMatMat(W, Uout, staple)
      counter0 = [site_linear_index(coord, dims4), 0_C64, 0_C64, 0_C64]
      call apply_su2_subgroup_update(Uout, W, 1, 2, 1, beta, key, counter0)
      ! current link, product * staple being updated in line with ULink,
      ! i1, i2 to identify which SU3 subgroup,
      ! subgroup ID used for Philox substream
      ! beta, Philox key
      ! counter with linearised site index
      call apply_su2_subgroup_update(Uout, W, 1, 3, 2, beta, key, counter0)
      call apply_su2_subgroup_update(Uout, W, 2, 3, 3, beta, key, counter0)
      call FixSU3Matrix(Uout)
    end function su3_updated_link

    pure subroutine apply_su2_subgroup_update(ULink, W, i1, i2, subgroup_id, beta, key, counter0)
      !$omp declare target
     !<
     !< Apply one embedded SU(2) heatbath update inside SU(3).
     !<
     !< Inputs:
     !<   ULink       : current SU(3) link being updated
     !<   W           : current product U * staple, updated in tandem with ULink
     !<   i1, i2      : the two SU(3) row/column indices defining the subgroup
     !<   subgroup_id : 1, 2, or 3, used for the Philox substream
     !<   beta, key   : heatbath parameters and Philox key
     !<   counter0    : counter with site id in counter(1)
     !<
     !< Method:
     !<   1) Extract the 2x2 block from W.
     !<   2) Project it to the SU(2) quaternion representation.
     !<   3) Build the SU(2) heatbath matrix.
     !<   4) Embed the result back into SU(3).
     !<   5) Left-multiply both ULink and W by the embedded matrix.
     !<
     complex(kind=WC), dimension(3, 3), intent(INOUT) :: ULink, W
     !< Input link, output link
      integer, intent(IN) :: i1, i2
      !< Select which 2x2 block
      integer, intent(IN) :: subgroup_id
      !< Philox variable
      real(kind=WP), intent(IN) :: beta
      !< beta variables
      integer(kind=C64), intent(IN), dimension(4) :: counter0
      !< Philox keys
      integer(kind=C64), intent(IN), dimension(2) :: key
      !< Philox key
      complex(kind=WC), dimension(2, 2) :: MfromW, SU2_M, SU2_X, SU2_U
      !< 2x2 (SU2) sub-matrices
      complex(kind=WC), dimension(3, 3) :: embed, Wtemp
      !< SU3 working matrices
      real(kind=WP), dimension(0:3) :: aQuart
      !< Quartenion representation
      real(kind=WP) :: alpha
      !< sqrt(sum(aQuart))
      real(kind=WP), parameter :: eps = TINY(1.0_WP)
      !< Smallest real number in WP precision
      MfromW(1, 1) = W(i1, i1)
      MfromW(2, 2) = W(i2, i2)
      MfromW(1, 2) = W(i1, i2)
      MfromW(2, 1) = W(i2, i1)
      aQuart(0) = real(MfromW(1, 1) + MfromW(2, 2), kind=WP)
      aQuart(1) = AIMAG(MfromW(1, 2) + MfromW(2, 1))
      aQuart(2) = real(MfromW(1, 2) - MfromW(2, 1), kind=WP)
      aQuart(3) = AIMAG(MfromW(1, 1) - MfromW(2, 2))
      aQuart = 0.5_WP * aQuart
      alpha = SQRT(SUM(aQuart**2))
      if (alpha <= eps) then
         SU2_U = constructXMatrix(1.0_WP, 2.0_WP * beta / 3.0_WP, key, counter0, subgroup_id)
      else
         aQuart = aQuart / alpha
         SU2_M = constructSU2Matrix(aQuart)
         SU2_X = constructXMatrix(alpha, 2.0_WP * beta / 3.0_WP, key, counter0, subgroup_id)
         SU2_M = CONJG(TRANSPOSE(SU2_M))
         SU2_U = MATMUL(SU2_X, SU2_M)
      end if
      embed = Ident3x3
      embed(i1, i1) = SU2_U(1, 1)
      embed(i2, i2) = SU2_U(2, 2)
      embed(i1, i2) = SU2_U(1, 2)
      embed(i2, i1) = SU2_U(2, 1)
      call MultiplyMatMat(ULink, embed, ULink)
      WTemp = W
      call MultiplyMatMat(W, embed, WTemp)
   end subroutine apply_su2_subgroup_update


   subroutine build_colour_sites(nt, nx, ny, nz, use_symanzik, &
        sites_it, sites_ix, sites_iy, sites_iz, counts)
     !< Construct checkerboard / multi-colour site lists
     !< The output layout is **structure-of-arrays (SoA)** for optimal
     !< memory access in `do concurrent` loops.
     !< Separates sites out into checkerboard
     !< Separate arrays for each of t, x, y, z
     integer, intent(IN) :: nt, nx, ny, nz
     !< Lattice dimensions
     logical, intent(IN) :: use_symanzik
     !< Selects checkerboard scheme
     integer, allocatable, intent(OUT), dimension(:, :, :) :: sites_it
     !< T: Output coordinate arrays (SoA layout)
     integer, allocatable, intent(OUT), dimension(:, :, :) :: sites_ix
     !< X: Output coordinate arrays (SoA layout)
     integer, allocatable, intent(OUT), dimension(:, :, :) :: sites_iy
     !< Y: Output coordinate arrays (SoA layout)
     integer, allocatable, intent(OUT), dimension(:, :, :) :: sites_iz
     !< Z: Output coordinate arrays (SoA layout)
     integer, allocatable, intent(OUT), dimension(:, :) :: counts
     !< Number of sites per colour/checkerboard and direction
     integer :: mu, colour, it, ix, iy, iz
     !< counters
     integer :: idx, ncolours, max_sites
     !< linearised index, checkerboard counter, lattice volume
      integer, dimension(4) :: coord
      !< lattice coordinate [it,ix,iy,iz]
      ncolours = ncolours_for_action(use_symanzik)
      max_sites = nt * nx * ny * nz
      allocate (sites_it(max_sites, ncolours, 4))
      allocate (sites_ix(max_sites, ncolours, 4))
      allocate (sites_iy(max_sites, ncolours, 4))
      allocate (sites_iz(max_sites, ncolours, 4))
      allocate (counts(ncolours, 4))
      counts = 0
      do mu = 1, 4
         do it = 1, nt
            do ix = 1, nx
               do iy = 1, ny
                  do iz = 1, nz
                     coord = (/it, ix, iy, iz/)
                     colour = site_colour(coord, mu, use_symanzik)
                     counts(colour + 1, mu) = counts(colour + 1, mu) + 1
                     idx = counts(colour + 1, mu)
                     sites_it(idx, colour + 1, mu) = it
                     sites_ix(idx, colour + 1, mu) = ix
                     sites_iy(idx, colour + 1, mu) = iy
                     sites_iz(idx, colour + 1, mu) = iz
                  end do
               end do
            end do
         end do
      end do
   end subroutine build_colour_sites

   subroutine updateLinks(U, beta, UUpdated, master_key, sweep_id, sites_it, sites_ix, sites_iy, sites_iz, counts, actionTag, xi)
     !< Perform one full heatbath sweep of the lattice
     !< ## Parallelism
     !< Uses `do concurrent` over site lists for each `(mu, colour)` stage.
     !<
     !< ## RNG strategy
     !< Each stage derives an independent Philox substream via
     !< `derive_stage_key`.
     complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(IN) :: U
     !< Input gaugefield [3,3,4,nt,nx,ny,nz]
      real(kind=WP), intent(IN) :: beta
      !< The beta value
      complex(kind=WC), dimension(:,:,:,:,:,:,:), intent(INOUT) :: UUpdated
      !< Updated gaugefield [3,3,4,nt,nx,ny,nz]
      integer(kind=C64), intent(IN), dimension(2) :: master_key
      !< Master key for philox
      integer, intent(IN) :: sweep_id
      !< which sweep/trajectory it is. Used for Philox
      integer, dimension(:,:,:), intent(IN) :: sites_it
      !< T: Output coordinates array in Structure of Arrays for better memomry access
      integer, dimension(:,:,:),intent(IN) :: sites_ix
      !< X: Output coordinates array in Structure of Arrays for better memomry access
      integer, dimension(:,:,:),intent(IN) :: sites_iy
      !< Y: Output coordinates array in Structure of Arrays for better memomry access
      integer, dimension(:,:,:),intent(IN) :: sites_iz
      !< Z: Output coordinates array in Structure of Arrays for better memomry access
      integer, dimension(:,:),intent(IN) :: counts
      !< Number of sites per colour and direction
      character(len=*), intent(IN), optional :: actionTag
      !< Specifies the action. Wilson|Symanzik|Iwasaki
      real(kind=WP), intent(in), optional :: xi
      !< gauge anisotropy
      !------------------------------------
      ! Locals
      !------------------------------------
      integer :: mu, colour, k
      !< counters
      integer :: it, ix, iy, iz
      !< counters
      integer :: ncolours
      !< Setup the checkerboarding layout. 2 for Wilson, 4 for Symanzik
      integer, dimension(4) :: dims4
      !< The lattice dimensions [nt,nx,ny,nz]
      logical :: use_symanzik
      !< Specifies if it is a symanzik-type action (needed for checkerboarding)
      integer(kind=C64), dimension(2) :: key
      !< Per checkerboard/direction Philox key
      complex(kind=WC), dimension(3,3) :: Uloc
      !< Local updated link
      integer, dimension(7) :: shapeU
      !< Shape of the gaugefield[3,3,4,nt,nx,ny,nz]
      real(kind=WP) :: xig
      !< The used gauge anisotropy - set to 1 if xi not present, else xi
      integer, dimension(4) :: coord
      !< current coordinate [it,ix,iy,iz]
      procedure(stapleInterface), pointer :: stapleKernel
      !< Pointer used to switch which action is used
      !------------------------------------
      ! Setup
      !------------------------------------
      shapeU = SHAPE(U)
      dims4 = shapeU(4:7)
      use_symanzik = .FALSE.
      if (PRESENT(actionTag)) then
         select case (to_lower(TRIM(actionTag)))
         case ('symanzik')
            use_symanzik = .TRUE.
            stapleKernel => stapleSymanzik
         case ('wilson')
            use_symanzik = .FALSE.
            stapleKernel => stapleWilson
         case ('iwasaki')
            use_symanzik = .true.
            stapleKernel => stapleIwasaki
         end select
      end if
      ncolours = ncolours_for_action(use_symanzik)
      UUpdated = U
      if (PRESENT(xi)) then
         xig = xi
      else
         xig = 1.0_WP
      end if
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! MAIN UPDATE
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do mu = 1, 4
         do colour = 0, ncolours - 1
            call derive_stage_key(master_key, sweep_id, &
                                  stage_tag=3, mu=mu, colour=colour, key=key)
            do concurrent(k=1:counts(colour + 1, mu)) LOCAL(it, ix, iy, iz, Uloc, coord)
               !------------------------------------
               ! Unit-stride coordinate loads (SoA)
               !------------------------------------
               it = sites_it(k, colour + 1, mu)
               ix = sites_ix(k, colour + 1, mu)
               iy = sites_iy(k, colour + 1, mu)
               iz = sites_iz(k, colour + 1, mu)
               coord = (/it, ix, iy, iz/)
               !------------------------------------
               ! Compute update (fully scalar args)
               !------------------------------------
               Uloc = su3_updated_link(UUpdated, beta, coord, mu, key, dims4, stapleKernel, xig)
               !Uloc = su3_updated_link(UUpdated, beta, coord, mu, key, dims4, xig)
               !Uloc = Ident3x3
               !------------------------------------
               ! Store result
               !------------------------------------
               UUpdated(:, :, mu, it, ix, iy, iz) = Uloc
            end do
         end do
      end do
   end subroutine updateLinks

   pure subroutine stapleWilson(U, V, coord, mu, xi)
     !$omp declare target
     !<
     !< Compute the Wilson staple for one site
     !<
     !<   V = sum_{nu != mu} [ path(nu,-mu,-nu) + path(-nu,-mu,nu) ]
     !< with the starting coordinate shifted to x + mu.
     !<
     complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(IN) :: U
     !< Input gaugefield [3,3,4,nt,nx,ny,nz]
      integer, dimension(4), intent(IN) :: coord
      !< Coordinate of the staple [it,ix,iy,iz]
      integer, intent(IN) :: mu
      !< Staple is along mu direction
      real(kind=WP), intent(in) :: xi
      !< The gauge anisotropy
      complex(kind=WC), dimension(3, 3), intent(OUT) :: V
      !< Output Wilson staple at this coordinate
      integer, dimension(4) :: thisCoord, step
      !< coordinate and step [it,ix,iy,iz], [...,mu,...] where otherwise mu is 0
      integer :: nu
      !< counter
      real(kind=WP) :: aniFac
      !< Used to set anisotropy factor for either space-space or space-time staples
      integer, dimension(7) :: dataShape
      !< Shape of the gaugefield [3,3,4,nt,nx,ny,nz]
      dataShape = SHAPE(U)
      step = 0
      step(mu) = 1
      thisCoord = coord + step
      thisCoord = periodCoord(thisCoord, dataShape)
      V = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      direction: do nu = 1, 4
         if (nu == mu) cycle direction
         if (mu == 1 .OR. nu == 1) then
            aniFac = xi
         else
            aniFac = 1.0_WP / xi
         end if
         V = V + aniFac * (genericPath(U, thisCoord, [nu, -mu, -nu]) &
                           + genericPath(U, thisCoord, [-nu, -mu, nu]))
      end do direction
   end subroutine stapleWilson

   pure subroutine stapleRectangle(U, V, coord, mu, xi)
     !$omp declare target
     !< Computes the rectangle staple (1x2) as in the Symanzik
     !< improved action for one site
     complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(IN) :: U
     !< Input gaugefield [3,3,4,nt,nx,ny,nz]
     integer, dimension(4), intent(IN) :: coord
     !< Coordinate of the staple [it,ix,iy,iz]
     integer, intent(IN) :: mu
     !< Staple is along mu direction
     real(kind=WP), intent(in) :: xi
     !< The gauge anisotropy
     complex(kind=WC), dimension(3, 3), intent(OUT) :: V
     !< Output rectangle staple at this coordinate
     integer, dimension(7) :: dataShape
     !< Shape of the gaugefield [3,3,4,nt,nx,ny,nz]
     integer, dimension(4) :: thisCoord, step
     !< coordinate and step [it,ix,iy,iz], [...,mu,..] where otherwise step is 0
     integer, dimension(5) :: r5
     !< Path for staples
     integer :: nu
     !< counter
     real(kind=WP) :: aniFac
     !< Used to set anisotropy factor for either space-space or space-time staples
     dataShape = SHAPE(U)
     V = CMPLX(0.0_WP, 0.0_WP, kind=WC)
     step = 0
     step(mu) = 1
     ! Reusing step to avoid array temporary's
     step = step + coord
      thisCoord = periodCoord(step, dataShape)  ! x + mu
      direction: do nu = 1, 4
         if (nu == mu) cycle direction
         if (mu == 1 .OR. nu == 1) then
            aniFac = xi
         else
            aniFac = 1.0_WP / xi
         end if
         ! Long in nu
         r5 = [nu, nu, -mu, -nu, -nu]
         V = V + genericPath(U, thisCoord, r5)
         r5 = [-nu, -nu, -mu, nu, nu]
         V = V + aniFac * genericPath(U, thisCoord, r5)
         ! Long in mu, link is the first mu
         r5 = [mu, nu, -mu, -mu, -nu]
         V = V + aniFac * genericPath(U, thisCoord, r5)
         r5 = [mu, -nu, -mu, -mu, nu]
         V = V + aniFac * genericPath(U, thisCoord, r5)
         ! Long in mu, link is the second mu
         r5 = [nu, -mu, -mu, -nu, mu]
         V = V + aniFac * genericPath(U, thisCoord, r5)
         r5 = [-nu, -mu, -mu, nu, mu]
         V = V + aniFac * genericPath(U, thisCoord, r5)
      end do direction
    end subroutine stapleRectangle

    pure subroutine stapleSymanzik(U, V, coord, mu, xi)
      !$omp declare target
     !<
     !< Compute the tree-level Symanzik staple for one SU(3) link.
     !<
     !< This consists of:
     !<   * the plaquette staple contribution (Wilson staple),
     !,   * minus the weighted rectangle contribution.
     !,
     !, With c0= 5/3, c1= -1/12
     !,
     complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(IN) :: U
     !< Input gaugefield [3,3,4,nt,nx,ny,nz]
      integer, dimension(4), intent(IN) :: coord
      !< Coordinate of the staple [it,ix,iy,iz]
      integer, intent(IN) :: mu
      !< Staple is along mu direction
      real(kind=WP), intent(in) :: xi
      !< The gauge anisotropy
      complex(kind=WC), dimension(3, 3), intent(OUT) :: V
      !< Output tree-level Symanzik staple at this coordinate
      complex(kind=WC), dimension(3, 3) :: Vplaq, Vrect
      !< Holders for plaquette, rectangle staples
      integer, dimension(7) :: dataShape
      !< Shape of the gaugefield [3,3,4,nt,nx,ny,nz]
      dataShape = SHAPE(U)
      call stapleWilson(U, Vplaq, coord, mu, xi)
      call stapleRectangle(U, Vrect, coord, mu, xi)
      V = (5.0_WP / 3.0_WP) * Vplaq - (1.0_WP / 12.0_WP) * Vrect
   end subroutine stapleSymanzik


   pure subroutine stapleIwasaki(U, V, coord, mu, xi)
     !$omp declare target
     !<
     !< Compute the Iwasaki staple for one SU(3) link.
     !<
     !< This consists of:
     !<   * the plaquette staple contribution (Wilson staple),
     !<   * minus the weighted rectangle contribution.
     !<
     !<
     !< With c0= 3.648, c1 = -0.331
     !<
     complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(IN) :: U
     !< Input gaugefield [3,3,4,nt,nx,ny,nz]
      integer, dimension(4), intent(IN) :: coord
      !< Coordinate of the staple [it,ix,iy,iz]
      integer, intent(IN) :: mu
      !< Staple is along mu direction
      real(kind=WP), intent(in) :: xi
      !< The gauge anisotropy
      complex(kind=WC), dimension(3, 3), intent(OUT) :: V
      !< Output Iwasaki staple at this coordinate
      complex(kind=WC), dimension(3, 3) :: Vplaq, Vrect
      !< Holders for plaquette, rectangle staples
      integer, dimension(7) :: dataShape
      !< Shape of the gaugefield [3,3,4,nt,nx,ny,nz]
      dataShape = SHAPE(U)
      call stapleWilson(U, Vplaq, coord, mu, xi)
      call stapleRectangle(U, Vrect, coord, mu, xi)
      V = 3.648_WP * Vplaq - 0.331_WP * Vrect
   end subroutine stapleIwasaki



end module FLUE_heatbath
