module FLUE_stoutSmearing
  !< Module to do stout-link smearing on SU(3) gaugelinks
  !< i.e. Morningstar & Peardon: hep-lat/0311018
   use FLUE_constants, only: WP, WC
   use FLUE_matrixConstants, only: Ident3x3
   use flue_wloops, only: genericPath
   use FLUE_SU3MatrixOps, only: MultiplyMatMatDag, ExpIQ, MultiplyMatMat, FixSU3Matrix, TraceMat
   implicit none(type, external)
   private

   public :: StoutSmearLinks

contains

  pure subroutine StoutSmearLinks(U, rho, nSweeps, USmeared)
    !< Stout-Smear the gaugelinks in U
    !< with a smearing weight of rho
    !< Do it nSweeps times
    !< output is in USmeared
    !< i.e. MorningStar & Peardon
    complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(in) :: U
    !< Input gaugefield [3,3,4,nt,nx,ny,nz]
    complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(out) :: USmeared
    !< Output smeared gaugefield [3,3,4,nt,nx,ny,nz]
    real(kind=WP), intent(in) :: rho
    !< Smearing strength
    integer, intent(in) :: nSweeps
    !< Number of sweeps of smearing
      ! working arrays
    complex(kind=WC), dimension(:, :, :, :, :, :, :), allocatable :: UTemp
    !< Temporary storage of smeared links
      ! Geometry
    integer, dimension(7) :: dataShape
    !< Shape of the gaugefield [3,3,4,nt,nx,ny,nz]
    integer :: nt, nx, ny, nz
    !< lattice dimensions
    integer :: ISweep
    !< counter
      dataShape = SHAPE(U)
      nt = dataShape(4)
      nx = dataShape(5)
      ny = dataShape(6)
      nz = dataShape(7)
      ! working arrays
      allocate (UTemp(3, 3, 4, nt, nx, ny, nz))
      ! Initialise with current field
      USmeared = U
      ! iterate over
      do iSweep = 1, nSweeps
         call StoutSmearOnce(USmeared, rho, UTemp)
         USmeared = UTemp
      end do
      deallocate (UTemp)
   end subroutine StoutSmearLinks

   pure subroutine StoutSmearOnce(U_old, rho, U_new)
     !< Does a single sweep of stout-link smearing
     complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(in) :: U_old
     !< Input gaugefield [3,3,4,nt,nx,ny,nz]
     complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(out) :: U_new
     !< Output gaugefield after 1 sweep of smearing [3,3,4,nt,nx,ny,nz]
     real(kind=WP), intent(in) :: rho
     !< Smering strength
      ! lattice geometry
     integer, dimension(7) :: dataShape
     !< Shape of the gaugefield [3,3,4,nt,nx,ny,nz]
     integer :: nx, ny, nz, nt
     !< Lattice dimensions
     integer :: it, ix, iy, iz, mu
     !< counters
     complex(kind=WC), dimension(:, :, :, :, :, :, :), allocatable :: staples
     !< Holds staples
     complex(kind=WC), dimension(3, 3) :: U_link, C, Q, V, U_updated
      !< for updating the links
      ! Get geometry
      dataShape = SHAPE(U_old)
      nt = dataShape(4)
      nx = dataShape(5)
      ny = dataShape(6)
      nz = dataShape(7)
      ! initialise output with input
      ! this is a quick way get the temporal links
      U_new = U_old
      ! Allocate and compute SPATIAL staples
      allocate (staples(3, 3, 3, nt, nx, ny, nz))
      call ComputeSpatialStaples(U_old, rho, Staples)

#ifdef LOCALITYSUPPORT
      do concurrent(mu=2:4, ix=1:nx, iy=1:ny, iz=1:nz, it=1:nt) &
         local(U_link, C, Q, V, U_updated) shared(U_old, U_new, Staples)
#else
         do mu = 2, 4
            do it = 1, nt
               do ix = 1, nx
                  do iy = 1, ny
                     do iz = 1, nz
#endif
                        ! Get the current link
                        U_link = U_old(:, :, mu, it, ix, iy, iz)
                        ! Get staple sum
                        C = Staples(:, :, mu - 1, it, ix, iy, iz)
                        ! Get the update term
                        Q = computeQMatrix(U_link, C)
                        V = ExpIQ(Q)
                        ! Update the link
                        call MultiplyMatMat(U_updated, V, U_link)
                        ! Reunitarise
                        call FixSU3Matrix(U_updated)
                        U_new(:, :, mu, it, ix, iy, iz) = U_updated
#ifdef LOCALITYSUPPORT
                     end do
#else
                  end do
               end do
            end do
         end do
      end do
#endif
      deallocate (staples)

   end subroutine StoutSmearOnce

   pure subroutine ComputeSpatialStaples(data, rho, staples)
     !< Computes *all* of the spatial staples for the gaugefield
     !< with factor of smearing strength rho included
     !< Computes sum of backward/forward staples at given site/direction
     complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(in) :: data
     !< Input gaugefield [3,3,4,nt,nx,ny,nz]
     real(kind=WP), intent(in) :: rho
     !< Smearing strength rho
      ! staples is only mu=1,2,3, i.e. spatial only
     complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(out) :: staples
     !< Output staples [3,3,3,nt,nx,ny,nz]
     integer, dimension(7) :: dataShape
     !< The shape of the gaugefield [3,3,4,nt,nx,ny,nz]
     integer, dimension(4) :: coord
     !< Single site coordiinate [it,ix,iy,iz]
     integer :: nt, nx, ny, nz
     !< Lattice dimensions
     integer :: it, iz, iy, ix, mu, nu
     !< counters
     complex(kind=WC), dimension(3, 3) :: stapleFwd, stapleBwd
     !< Per site/direction forward and backward staples
      dataShape = SHAPE(data)
      nt = dataShape(4)
      nx = dataShape(5)
      ny = dataShape(6)
      nz = dataShape(7)
      staples = CMPLX(0.0_WP, 0.0_WP, kind=WC)
#ifdef LOCALITYSUPPORT
      do concurrent(ix=1:nx, iy=1:ny, iz=1:nz, it=1:nt) &
         local(coord, mu, nu, stapleFwd, stapleBwd) shared(data, staples)
#else
         do it = 1, nt
            do ix = 1, nx
               do iy = 1, ny
                  do iz = 1, nz
#endif
                     coord = (/it, ix, iy, iz/)
                     ! do spatial links only
                     do mu = 2, 4
                        ! do staples from all other directions
                        direction: do nu = 2, 4
                           if (nu == mu) cycle direction
                           ! Calculate forward staple
                           stapleFwd = genericPath(data, coord, (/nu, mu, -nu/))
                           staples(:, :, mu - 1, it, ix, iy, iz) = staples(:, :, mu - 1, it, ix, iy, iz) &
                                                                   + stapleFwd
                           ! and the backward staple
                           stapleBwd = genericPath(data, coord, (/-nu, mu, nu/))
                           staples(:, :, mu - 1, it, ix, iy, iz) = staples(:, :, mu - 1, it, ix, iy, iz) &
                                                                   + stapleBwd
                        end do direction
                     end do
#ifdef LOCALITYSUPPORT
                  end do
#else
               end do
            end do
         end do
      end do
#endif
      staples = rho * staples
   end subroutine ComputeSpatialStaples

   pure function ComputeQMatrix(U, C) result(Q)
     !< Computes the Q Matrix required
     !< Eqn2 of Morningstar & Peardon
     complex(kind=WC), dimension(3, 3), intent(in) :: U, C
     !< Input gauge link and sum of staples
     complex(kind=WC), dimension(3, 3) :: Q
     !< Output Q Matrix
     complex(kind=WC), dimension(3, 3) :: Omega, Diff
     !< Omega = CU^\dagger. Diff is helper variable
      ! C is the raw staple sum.
      ! i.e. Eqn1 of Morningstar & Peadron
      call MultiplyMatMatDag(Omega, C, U)
      Diff = CONJG(TRANSPOSE(Omega)) - Omega
      Q = CMPLX(0.0_WP, 0.5_WP, kind=WC) * (Diff - Ident3x3 * (1.0_WP / 3.0_WP) * TraceMat(Diff))
   end function computeQMatrix

end module FLUE_stoutSmearing
