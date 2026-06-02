module FLUE_SU2_heatbath
   use FLUE_constants, only: WP, WC, TWO_PI
   use FLUE_matrixConstants, only: Ident2x2
   use FLUE_SU2_random, only: constructSU2Matrix
   use FLUE_SU2_wloops, only: SU2_genericPath
   use FLUE_wloops, only: periodCoord
   use stdlib_stats_distribution_uniform, only: rvs_uniform
   implicit none(type, external)
   private
   public :: SU2_updateLinks
contains
   function getLambda2(alpha, beta) result(lambda2)
      real(kind=WP), intent(in) :: alpha, beta
      real(kind=WP) :: lambda2, rr
      real(kind=WP) :: ri(3)
      real(kind=WP), parameter :: eps = TINY(1.0_WP)
      ! Defensive: beta or alpha too small -> distribution becomes Haar-like.
      ! Caller should handle beta~0 separately, but we also guard here.
      if (alpha <= eps .OR. beta <= eps) then
         lambda2 = 0.0_WP
         return
      end if
      do
         ri = rvs_uniform(loc=0.0_WP, scale=1.0_WP, array_size=3)
         if (ANY(ri <= eps)) cycle  ! avoid log(0)
         ! Eq. (4.45): proposal for lambda^2
         lambda2 = -(LOG(ri(1)) + LOG(ri(3)) * COS(TWO_PI * ri(2))**2) &
                   / (2.0_WP * alpha * beta)
         ! Make sure inside bounds
         if (lambda2 < 0.0_WP) cycle
         if (lambda2 > 1.0_WP) cycle
         ! Eq. (4.46): accept with prob sqrt(1 - lambda^2)
         rr = rvs_uniform(scale=1.0_WP)
         if (rr <= eps) cycle
         ! done?
         if (rr * rr <= 1.0_WP - lambda2) exit
      end do
   end function getLambda2

   function getXVec(x0) result(xvec)
      real(kind=WP), intent(in) :: x0
      real(kind=WP) :: xvec(3)
      real(kind=WP) :: xlen, requiredLen
      real(kind=WP), parameter :: eps = TINY(1.0_WP)
      requiredLen = MAX(0.0_WP, 1.0_WP - x0 * x0)          ! requiredLen = 1 - x0^2
      if (requiredLen <= eps) then
         xvec = 0.0_WP
         return
      end if
      do
         xvec = rvs_uniform(loc=-1.0_WP, scale=2.0_WP, array_size=3)
         ! Reject exactly +1 only if you insist on [-1,1) — but using eps is safer
         if (ANY(xvec >= 1.0_WP - eps)) cycle
         ! check conditions
         xlen = SUM(xvec**2.0_WP)
         if (xlen <= 1.0_WP .AND. xlen > eps) exit
      end do
      ! Rescale to |xvec| = sqrt(1 - x0^2)
      xvec = xvec * (SQRT(requiredLen) / SQRT(xlen))
   end function getXVec

   function constructXMatrix(alpha, beta) result(X)
    !! section 4.3.1
      real(kind=WP), intent(in) :: alpha, beta
      complex(kind=WC), dimension(2, 2) :: X
      real(kind=WP) :: lambda2
      real(kind=WP), dimension(0:3) :: xAll
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
      end function constructXMatrix

      subroutine SU2_updateLinks(U, beta, UUpdated)
         complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(in) :: U
         real(kind=WP), intent(in) :: beta
         complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(inout) :: UUpdated
         !complex(kind=WC), allocatable, dimension(:,:,:,:,:,:,:) :: staples
         complex(kind=WC), dimension(2, 2) :: XMatrix, V, Vdag
         complex(kind=WC) :: detV
         real(kind=WP) :: alpha
    !! lattice geometry
         integer, dimension(7) :: dataShape
         integer :: nt, nx, ny, nz
    !! counters
         integer :: ix, iy, iz, it, mu
         integer, dimension(4) :: coord
         real(kind=WP), parameter :: eps = TINY(1.0_WP)
         dataShape = SHAPE(U)
         nt = dataShape(4)
         nx = dataShape(5)
         ny = dataShape(6)
         nz = dataShape(7)
         !allocate(staples(2, 2, 4, nt, nx, ny, nz))
         !call computeStaples(U, staples)
         UUpdated = U
         if (beta .LE. eps) then
            write (*, *) 'beta', beta, ' must not be zero. stopping'
            stop
         end if
         do it = 1, nt
            do ix = 1, nx
               do iy = 1, ny
                  do iz = 1, nz
                     coord = (/it, ix, iy, iz/)
                     do mu = 1, 4
                        call stapleAt(UUpdated, V, coord, mu)
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
                        if (alpha < eps) then
                           ! Degenerate staple: fall back to Haar-ish update (or identity)
                           XMatrix = constructXMatrix(1.0_WP, beta)    ! mild fallback
                           UUpdated(:, :, mu, it, ix, iy, iz) = XMatrix
                        else
                           ! Vtilde = V / alpha, and we need Vtilde^dagger on the right
                           Vdag = CONJG(TRANSPOSE(V / alpha))
                           !Vdag = V / alpha
                           XMatrix = constructXMatrix(alpha, beta)
                           ! Update: U' = X * Vtilde^dagger  (Montvay–Münster / Gattringer–Lang)
                           UUpdated(:, :, mu, it, ix, iy, iz) = MATMUL(XMatrix, Vdag)
                        end if
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
                     end do
                  end do
               end do
            end do
         end do

      end subroutine SU2_updateLinks

      pure subroutine stapleAt(U, V, coord, mu)
         complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(in) :: U
         integer, dimension(4), intent(in) :: coord
         integer, intent(in) :: mu
         complex(kind=WC), dimension(2, 2), intent(out) :: V
         integer, dimension(4) :: thisCoord, step
         integer :: nu
         step = 0
         step(mu) = 1
         thisCoord = coord + step
         thisCoord = periodCoord(thisCoord, SHAPE(U))
         V = CMPLX(0.0_WP, 0.0_WP, kind=WC)
         do nu = 1, 4
            if (nu == mu) cycle
            V = V + SU2_genericPath(U, thisCoord, (/nu, -mu, -nu/)) &
                + SU2_genericPath(U, thisCoord, (/-nu, -mu, nu/))
         end do
      end subroutine stapleAt

   end module FLUE_SU2_heatbath
