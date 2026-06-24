module FLUE_heatbath
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
      pure subroutine stapleInterface(U, V, coord, mu, xi)
         import :: WP, WC
         implicit none(type, external)
         complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(IN) :: U
         integer, dimension(4), intent(IN) :: coord
         integer, intent(IN) :: mu
         real(kind=WP), intent(in) :: xi
         complex(kind=WC), dimension(3, 3), intent(OUT) :: V
      end subroutine stapleInterface
   end interface

contains

   ! The force inline seems to help speed sometimes?
  !!!d ir$ forceinline
   pure function su3_updated_link(U, beta, coord, mu, key, dims4, use_symanzik, xi) result(Uout)
      complex(kind=WC), intent(IN) :: U(:, :, :, :, :, :, :)
      real(kind=WP), intent(IN) :: beta
      integer, intent(IN) :: mu
      integer, intent(IN) :: dims4(4)
      logical, intent(IN) :: use_symanzik
      integer, dimension(4), intent(in) :: coord
      real(kind=WP), intent(in) :: xi
      integer(C64), intent(IN) :: key(2)
      complex(kind=WC) :: Uout(3, 3)

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
      complex(kind=WC), dimension(3, 3) :: staple, W
      integer(C64) :: counter0(4)
      Uout = U(:, :, mu, coord(1), coord(2), coord(3), coord(4))

      if (use_symanzik) then
         call stapleSymanzik(U, staple, coord, mu, xi)
      else
         call stapleWilson(U, staple, coord, mu, xi)
      end if
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
      complex(kind=WC), dimension(3, 3), intent(INOUT) :: ULink, W
      integer, intent(IN) :: i1, i2, subgroup_id
      real(kind=WP), intent(IN) :: beta
      integer(C64), intent(IN) :: key(2), counter0(4)
      complex(kind=WC), dimension(2, 2) :: MfromW, SU2_M, SU2_X, SU2_U
      complex(kind=WC), dimension(3, 3) :: embed, Wtemp
      real(kind=WP), dimension(0:3) :: aQuart
      real(kind=WP) :: alpha
      real(kind=WP), parameter :: eps = TINY(1.0_WP)
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
     !! Separates sites out into checkerboard
     !! Separate arrays for each of t, x, y, z
      integer, intent(IN) :: nt, nx, ny, nz
      logical, intent(IN) :: use_symanzik
      integer, allocatable, intent(OUT) :: sites_it(:, :, :)
      integer, allocatable, intent(OUT) :: sites_ix(:, :, :)
      integer, allocatable, intent(OUT) :: sites_iy(:, :, :)
      integer, allocatable, intent(OUT) :: sites_iz(:, :, :)
      integer, allocatable, intent(OUT) :: counts(:, :)
      integer :: mu, colour, it, ix, iy, iz
      integer :: idx, ncolours, max_sites
      integer, dimension(4) :: coord
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
      !------------------------------------
      ! Arguments
      !------------------------------------
      complex(kind=WC), intent(IN) :: U(:, :, :, :, :, :, :)
      real(kind=WP), intent(IN) :: beta
      complex(kind=WC), intent(INOUT) :: UUpdated(:, :, :, :, :, :, :)
      integer(C64), intent(IN) :: master_key(2)
      integer, intent(IN) :: sweep_id
      integer, intent(IN) :: sites_it(:, :, :)
      integer, intent(IN) :: sites_ix(:, :, :)
      integer, intent(IN) :: sites_iy(:, :, :)
      integer, intent(IN) :: sites_iz(:, :, :)
      integer, intent(IN) :: counts(:, :)
      character(len=*), intent(IN), optional :: actionTag
      real(kind=WP), intent(in), optional :: xi
      !------------------------------------
      ! Locals
      !------------------------------------
      integer :: mu, colour, k
      integer :: it, ix, iy, iz
      integer :: ncolours
      integer :: dims4(4)
      logical :: use_symanzik
      integer(C64) :: key(2)
      complex(kind=WC) :: Uloc(3, 3)
      integer :: shapeU(7)
      real(kind=WP) :: xig
      integer, dimension(4) :: coord
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
         case ('wilson')
            use_symanzik = .FALSE.
         end select
      end if
      ncolours = ncolours_for_action(use_symanzik)
      UUpdated = U
      if (PRESENT(xi)) then
         xig = xi
      else
         xig = 1.0_WP
      end if
      !====================================
      ! MAIN UPDATE
      !====================================
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
               Uloc = su3_updated_link(UUpdated, beta, coord, mu, key, dims4, use_symanzik, xig)
               !------------------------------------
               ! Store result
               !------------------------------------
               UUpdated(:, :, mu, it, ix, iy, iz) = Uloc
            end do
         end do
      end do
   end subroutine updateLinks

   pure subroutine stapleWilson(U, V, coord, mu, xi)
     !!
     !! Compute the Wilson staple for one SU(3) link.
     !!
     !!   V = sum_{nu != mu} [ path(nu,-mu,-nu) + path(-nu,-mu,nu) ]
     !! with the starting coordinate shifted to x + mu.
     !!
      complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(IN) :: U
      integer, dimension(4), intent(IN) :: coord
      integer, intent(IN) :: mu
      real(kind=WP), intent(in) :: xi
      complex(kind=WC), dimension(3, 3), intent(OUT) :: V
      integer, dimension(4) :: thisCoord, step
      integer :: nu
      real(kind=WP) :: aniFac
      integer, dimension(7) :: dataShape
      dataShape = SHAPE(U)
      step = 0
      step(mu) = 1
      thisCoord = coord + step
      thisCoord = periodCoord(thisCoord, dataShape)
      V = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      do nu = 1, 4
         if (nu == mu) cycle
         if (mu == 1 .OR. nu == 1) then
            aniFac = xi
         else
            aniFac = 1.0_WP / xi
         end if
         V = V + aniFac * (genericPath(U, thisCoord, [nu, -mu, -nu]) &
                           + genericPath(U, thisCoord, [-nu, -mu, nu]))
      end do
   end subroutine stapleWilson

   pure subroutine stapleSymanzik(U, V, coord, mu, xi)
     !!
     !! Compute the tree-level Symanzik staple for one SU(3) link.
     !!
     !! This consists of:
     !!   * the plaquette staple contribution (Wilson staple),
     !!   * minus the weighted rectangle contribution.
     !!
      complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(IN) :: U
      integer, dimension(4), intent(IN) :: coord
      integer, intent(IN) :: mu
      real(kind=WP), intent(in) :: xi
      complex(kind=WC), dimension(3, 3), intent(OUT) :: V
      integer, dimension(4) :: thisCoord, step
      integer, dimension(5) :: r5
      complex(kind=WC), dimension(3, 3) :: Vplaq, Vrect
      integer :: nu
      real(kind=WP) :: aniFac
      integer, dimension(7) :: dataShape
      dataShape = SHAPE(U)
      call stapleWilson(U, Vplaq, coord, mu, xi)
      Vrect = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      step = 0
      step(mu) = 1
      ! Reusing step to avoid array temporary's
      step = step + coord
      thisCoord = periodCoord(step, dataShape)  ! x + mu
      do nu = 1, 4
         if (nu == mu) cycle
         if (mu == 1 .OR. nu == 1) then
            aniFac = xi
         else
            aniFac = 1.0_WP / xi
         end if
         ! Long in nu
         r5 = [nu, nu, -mu, -nu, -nu]
         Vrect = Vrect + genericPath(U, thisCoord, r5)
         r5 = [-nu, -nu, -mu, nu, nu]
         Vrect = Vrect + aniFac * genericPath(U, thisCoord, r5)
         ! Long in mu, link is the first mu
         r5 = [mu, nu, -mu, -mu, -nu]
         Vrect = Vrect + aniFac * genericPath(U, thisCoord, r5)
         r5 = [mu, -nu, -mu, -mu, nu]
         Vrect = Vrect + aniFac * genericPath(U, thisCoord, r5)
         ! Long in mu, link is the second mu
         r5 = [nu, -mu, -mu, -nu, mu]
         Vrect = Vrect + aniFac * genericPath(U, thisCoord, r5)
         r5 = [-nu, -mu, -mu, nu, mu]
         Vrect = Vrect + aniFac * genericPath(U, thisCoord, r5)
      end do
      V = (5.0_WP / 3.0_WP) * Vplaq - (1.0_WP / 12.0_WP) * Vrect
   end subroutine stapleSymanzik

end module FLUE_heatbath
