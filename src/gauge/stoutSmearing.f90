module FLUE_stoutSmearing
use FLUE_constants, only: WP, WC
use FLUE_matrixConstants, only: Ident3x3
use flue_wloops, only: genericPath
use FLUE_SU3MatrixOps, only: MultiplyMatMatDag, ExpIQ, MultiplyMatMat, FixSU3Matrix, TraceMat
use stdlib_linalg, only: is_hermitian
implicit none (type, external)
private

public :: StoutSmearLinks

contains

pure subroutine StoutSmearLinks(U, rho, nSweeps, USmeared)
    complex(kind=WC), dimension(:,:,:,:,:,:,:), intent(in) :: U
    complex(kind=WC), dimension(:,:,:,:,:,:,:), intent(out) :: USmeared
    real(kind=WP), intent(in) :: rho
    integer, intent(in) :: nSweeps
    ! working arrays
    complex(kind=WC), dimension(:,:,:,:,:,:,:), allocatable :: UTemp
    ! Geometry
    integer, dimension(7) :: dataShape
    integer :: nt, nx, ny, nz
    ! counters
    integer :: ISweep
    dataShape = shape(U)
    nt = dataShape(4)
    nx = dataShape(5)
    ny = dataShape(6)
    nz = dataShape(7)
    ! working arrays
    allocate(UTemp(3, 3, 4, nt, nx, ny, nz))
    ! Initialise with current field
    USmeared = U
    ! iterate over
    do iSweep =1, nSweeps
        call StoutSmearOnce(USmeared, rho, UTemp)
        USmeared = UTemp
    end do
    deallocate(UTemp)
end subroutine StoutSmearLinks

pure subroutine StoutSmearOnce(U_old, rho, U_new)
    complex(kind=WC), dimension(:,:,:,:,:,:,:), intent(in) :: U_old
    complex(kind=WC), dimension(:,:,:,:,:,:,:), intent(out) :: U_new
    real(kind=WP), intent(in) :: rho
    ! lattice geometry
    integer, dimension(7) :: dataShape
    integer :: nx, ny, nz, nt
    ! counters
    integer :: it, ix, iy, iz, mu
    ! staples holder
    complex(kind=WC), dimension(:,:,:,:,:,:,:), allocatable :: staples
    ! for updating the links
    complex(kind=WC), dimension(3,3) :: U_link, C, Q, V, U_updated

    ! Get geometry
    dataShape = shape(U_old)
    nt = dataShape(4)
    nx = dataShape(5)
    ny = dataShape(6)
    nz = dataShape(7)
    ! initialise output with input
    ! this is a quick way get the temporal links
    U_new = U_old
    ! Allocate and compute SPATIAL staples
    allocate(staples(3,3,3,nt,nx,ny,nz))
    call ComputeSpatialStaples(U_old, rho, Staples)

#ifdef LOCALITYSUPPORT
    do concurrent(mu=2:4, ix=1:nx, iy=1:ny, iz=1:nz, it=1:nt) &
        local(U_link, C, Q, V, U_updated) shared(U_old, U_new, Staples)
#else
    do mu=2, 4
        do it=1, nt
            do ix=1, nx
                do iy=1, ny
                    do iz=1, nz
#endif
    ! Get the current link
    U_link = U_old(:, :, mu, it, ix, iy, iz)
    ! Get staple sum
    C = Staples(:, :, mu-1, it, ix, iy, iz)
    ! Get the update term
    Q = computeQMatrix(U_link, C)
    V = ExpIQ(Q)
    ! Update the link
    call MultiplyMatMat(U_updated, V, U_link)
    ! Reunitarise
    call FixSU3Matrix(U_updated)
    U_new(:, :, mu, it,ix,iy,iz) = U_updated
#ifdef LOCALITYSUPPORT
    end do
#else
                    end do
                end do
            end do
        end do
    end do
#endif
deallocate(staples)

end subroutine StoutSmearOnce

pure subroutine ComputeSpatialStaples(data, rho, staples)
    complex(kind=WC), dimension(:,:,:,:,:,:,:), intent(in)  :: data
    real(kind=WP), intent(in) :: rho
    ! staples is only mu=1,2,3, i.e. spatial only
    complex(kind=WC), dimension(:,:,:,:,:,:,:), intent(out) :: staples
    ! Lattice dimensions
    integer, dimension(7) :: dataShape
    integer, dimension(4) :: coord
    integer :: nt, nx, ny, nz
    ! counters
    integer :: it, iz, iy, ix, mu, nu
    ! holders
    complex(kind=WC), dimension(3,3) :: stapleFwd, stapleBwd
    dataShape = shape(data)
    nt = dataShape(4)
    nx = dataShape(5)
    ny = dataShape(6)
    nz = dataShape(7)

    staples = cmplx(0.0_WP, 0.0_WP, kind=WC)
#ifdef LOCALITYSUPPORT
    do concurrent(ix=1: nx, iy=1:ny, iz=1:nz, it=1:nt) &
        local(coord, mu, nu, stapleFwd, stapleBwd) shared(data, staples)
#else
    do it=1, nt
        do ix=1, nx
            do iy=1, ny
                do iz=1, nz
#endif
    coord = (/it, ix, iy, iz/)
    ! do spatial links only
    do mu= 2, 4
        ! do staples from all other directions
        do nu=2, 4
            if (nu == mu) cycle
            ! Calculate forward staple
            stapleFwd = genericPath(data, coord, (/nu, mu, -nu/))
            staples(:,:,mu-1,it,ix,iy,iz) = staples(:,:,mu-1,it,ix,iy,iz) &
            + stapleFwd
            ! and the backward staple
            stapleBwd = genericPath(data, coord, (/-nu, mu, nu/))
            staples(:,:,mu-1,it,ix,iy,iz) = staples(:,:,mu-1,it,ix,iy,iz) &
                + stapleBwd
        end do
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
    complex(kind=WC), dimension(3,3), intent(in) :: U, C
    complex(kind=WC), dimension(3,3) :: Q
    ! working matricesb
    complex(kind=WC), dimension(3,3) :: Omega, Diff
    ! C is the raw staple sum.
    ! i.e. Eqn1 of Morningstar & Peadron
    call MultiplyMatMatDag(Omega, C, U)
    Diff = conjg(transpose(Omega)) - Omega
    Q = cmplx(0.0_WP, 0.5_WP, kind=WC) * (Diff - Ident3x3 * (1.0_WP/3.0_WP) * TraceMat(Diff))
end function computeQMatrix

end module FLUE_stoutSmearing
