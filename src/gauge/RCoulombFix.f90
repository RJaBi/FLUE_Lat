module FLUE_gaugeFix_Coulomb
  !< Module FLUE_gaugeFix_Coulomb
  !< Contains the functional, residual & gauge fixing observable theta
  !< Current implementation uses a standard CG to apply the inverse Laplacian,
  !< in a similar way to as mentioned in arXiv:hep-lat/9711047
  use FLUE_constants, only: WP, WC, assertSTDLIBKindCompatibility
  use FLUE_SU3MatrixOps, only: RealTraceMat, TraceLessConjgSubtract, realTraceMultMatMatDag, &
       colourDecompResidual, adjointToResidual, expIQ
  use FLUE_wloops, only: periodCoord
  use stdlib_kinds, only: DP
  use stdlib_linalg_iterative_solvers, only: &
       stdlib_linop_dp_type,                    &
       stdlib_solver_workspace_dp_type,         &
       stdlib_solve_cg_kernel,                  &
       stdlib_size_wksp_cg

  implicit none(type, external)
  private
  public :: Coulomb_functional, Coulomb_residual, Coulomb_theta, Coulomb_constructUpdateField
contains
  function Coulomb_functional(U) result(funcVal)
    !< The value of the coulomb functional
    !< $\sum_{x}\sum_{\mu=2}^{4} ReTr[U_\mu(x)]$
    !< Note this is a sum over spatial directions
    !< and that we just take the real part because it should be real
    !< Functional normallised by (3*3*nt*nx*ny*nz)
    complex(kind=WC), dimension(:,:,:,:,:,:,:), intent(in) :: U
    !< The gaugefield of shape 3,3,4,nt,nx,ny,nz
    real(kind=WP) :: funcVal
    !< Output value of functional
    integer :: mu, ix, iy, iz, it
    !< counters
    integer :: nt, nx, ny , nz
    !< lattice dimensions
    complex(kind=WC), dimension(3,3) :: UHold
    !< Holds each $U_\mu(x) + U_\mu(x)^\dagger$
    real(kind=WP) :: contribution
    !< Per lattice site value
    nt = size(U, 4)
    nx = size(U, 5)
    ny = size(U, 6)
    nz = size(U, 7)
    funcVal = 0.0_WP
    do concurrent(it=1:nt, ix=1:nx, iy=1:ny, iz=1:nz, mu=2:4) &
         default(none) reduce(+:funcVal) local(UHold, contribution) shared(U)
       UHold = U(:,:,mu,it,ix,iy,iz)
       contribution = RealTraceMat(UHold)
       funcVal = funcVal + contribution
    end do
    ! Normallise functional
    funcVal = funcVal / real(9*nt*nx*ny*nz, kind=WP)
  end function Coulomb_functional

   subroutine Coulomb_residual(U, deltaX)
     !< Construct anti-Hermitian traceless Coulomb Fix Residual
     !< Construct at each site x
     !< $delta^G(x) = \sum_{\mu} (A^G_\mu(x) - A^G_\mu(x-\hat{\mu}))$
     !< where
     !< $A_\mu(x) = \frac{1}{2i}(U_\mu(x) - U_\mu^\dagger(x)) traceless$
     !< but actually it doesn't have factor {1}{2i} so that it is actually anti-hermitian
     !< $A_\mu(x) = \frac{1}{2i}(U_\mu(x) - U_\mu^\dagger(x)) traceless$
     complex(kind=WC), dimension(:,:,:,:,:,:,:), intent(in) :: U
     !< The gauge-fixed gaugefield of shape 3,3,4,nt,nx,ny,nz
     complex(kind=WC), dimension(:,:,:,:,:,:), intent(out) :: deltaX
     !< The residual of shape 3,3 nt,nx,ny,nz
     integer :: mu, it, ix, iy, iz
     !< counters
     integer :: nt, nx, ny, nz
     !< lattice dimensions
     integer, dimension(7) :: dataShape
     !< gaugefield dimensions
     complex(kind=WC), dimension(3,3) :: UHold
     !< Holder for SU3 matrix
     integer, dimension(4) :: coord, pathCoord
     !< Lattice coordinate, step coordinate
     deltaX = 0.0_WP
     dataShape = shape(U)
     nt = size(U, 4)
     nx = size(U, 5)
     ny = size(U, 6)
     nz = size(U, 7)
     UHold = cmplx(8.0, 4.0)
     do mu=2, 4
        pathCoord = 0
        pathCoord(mu) = -1
        do concurrent(it=1:nt, ix=1:nx, iy=1:ny, iz=1:nz) default(none)&
             local(UHold, coord) shared(mu, pathCoord, U, deltaX, dataShape)
           call TraceLessConjgSubtract(UHold, U(:,:,mu,it,ix,iy,iz), &
                U(:,:,mu,it,ix,iy,iz))
           deltaX(:,:, it, ix, iy, iz) = deltaX(:,:, it, ix, iy, iz) + UHold
           coord = (/it, ix, iy, iz/) + pathCoord
           coord = periodCoord(coord, dataShape)
           call TraceLessConjgSubtract(UHold, &
                U(:,:,mu,coord(1), coord(2), coord(3), coord(4)), &
                U(:,:,mu,coord(1), coord(2), coord(3), coord(4)))
           deltaX(:,:, it, ix, iy, iz) = deltaX(:,:, it, ix, iy, iz) - UHold
           end do
     end do
     deltaX = 0.5_WP * deltaX
   end subroutine Coulomb_residual

   subroutine Coulomb_theta(deltaX, thetaTimeslice)
     !< Measures theta
     !< Eqn 25
     !< $theta = \frac{1}{V} \sum_{x} Tr(\delta^G(x)\delta^{\dagger,G}(x))$
     complex(kind=WC), dimension(:,:,:,:,:,:), intent(in) :: deltaX
     !< Residual quantity of Eqn 20
     !< Theta, Eqn 25
     real(kind=WP), dimension(:), intent(out) :: thetaTimeslice
     !< Theta, Eqn 25 per timeslice
     integer :: nt, nx, ny, nz
     !< lattice dimensions
     real(kind=WP) :: perVal
     !< inside loop value
     integer :: it
     !< counter
     nt = size(deltaX, 3)
     nx = size(deltaX, 4)
     ny = size(deltaX, 5)
     nz = size(deltaX, 6)
     thetaTimeSlice = 0.0_WP
     do concurrent(it=1:nt) default(none) local(perVal) shared(deltaX, thetaTimeSlice, nx, ny, nz)
        perVal = sum(real(deltaX(:,:, it, :, :, :), kind=WP)**2.0_WP &
             + aimag(deltaX(:,:, it, :, :, :))**2.0_WP)
        perVal = perVal / real(3 * nx * ny * nz, kind=WP)
        thetaTimeslice(it) = perVal
     end do
   end subroutine Coulomb_theta

   subroutine Coulomb_constructUpdateField(deltaX, innerRelTol, innerAbsTol, maxInnerIterations, innerIterations, updateField)
     !< Computes the gauge-update field $\omega(x)$ for eventual use in
     !< $g(x) = \text{exp}(a\,\Omega(x))$
     !< where $\Omega$ is the 3x3 colour form of the adjoint-colour $\omega(x)$
     !< Update field is obtained from
     !< $-\nabla^2\omega(x) = \delta(x)$
     !< after $\delta(x)$ is put into adjoint-colour space
     !< This equation is solved via CG
     complex(kind=WC), dimension(:,:,:,:,:,:), intent(in) :: deltaX
     !< $\delta(x)$. Gauge fixing residual 3,3,nt,nx,ny,nz
     real(kind=WP), intent(in) :: innerRelTol
     !< Inner relative tolerance (CG)
     real(kind=WP), intent(in) :: innerAbsTol
     !< Inner absolute tolerance (CG)
     integer, intent(in) :: maxInnerIterations
     !< Maximum number of CG iterations (single step)
     integer, intent(inout) :: innerIterations
     !< total number of CG iterations
     complex(kind=WC), dimension(:,:,:,:,:,:), intent(out) :: updateField
     !< output field $\Omega(x)$
     ! TODO MAKE WORKSPACE
     real(kind=WP), dimension(:, :, :, :, :), allocatable :: rhs
     !< workspace 8,nt,nx,ny,nz in adjoint colour
     real(kind=WP), dimension(:, :, :, :, :), allocatable :: solution
     !< workspace 8,nt,nx,ny,nz in adjoint colour
     real(kind=DP), dimension(:), allocatable :: rhsLinear, solutionLinear
     !< linearised workspace 8*nx*ny*nz
     integer :: timeSliceIterations
     !< number of iterations at single time slice
     integer :: nt, nx, ny, nz
     !< Lattice dimensions
     integer :: it
     !< counters
     nt = size(deltaX, 3)
     nx = size(deltaX, 4)
     ny = size(deltaX, 5)
     nz = size(deltaX, 6)
     allocate(rhs(8, nt, nx, ny, nz))
     allocate(solution(8, nt, nx, ny, nz))
     allocate(rhsLinear(size(rhs)/nt), solutionLinear(size(solution)/nt))
     !
     ! The Coulomb residual deltaX is anti-Hermitian. Convert
     !
     !    Q_residual = -i deltaX
     !
     ! into eight real coefficients in the T_a = lambda_a/2 basis.
     call residualFieldToAdjoint(deltaX, rhs)
     !
     ! The periodic spatial Laplacian has one constant zero mode for each algebra component.
     ! The Poisson problem is solved on the zero-mean subspace
     call projectCoulombZeroModes(rhs)
     solution = 0.0_WP
     do it = 1, nt
        ! Start CG from zero. This also ensures that the initial guess lies
        ! in the zero-mean subspace.
        !
        solutionLinear = 0.0_DP
        rhsLinear = real(reshape(rhs(:,it,:,:,:), (/size(rhsLinear)/)), kind=DP)
        call solveLaplacianTimeslice( &
             rhsLinear, &
             solutionLinear, &
             real(innerRelTol, kind=DP), &
             real(innerAbsTol, kind=DP), &
             maxInnerIterations, &
             nx, ny, nz, &
             timeSliceIterations)
        solution(:, it, :, :, :) = real(reshape(solutionLinear, (/8, nx, ny, nz/)), kind=WP)
        innerIterations = &
             innerIterations + timesliceIterations
     end do
     call projectCoulombZeroModes(solution)
     ! Take it from 8 real value coefficients to 3x3 matrix
     call adjointFieldToResidual(solution, updateField)
   end subroutine Coulomb_constructUpdateField

   pure subroutine residualFieldToAdjoint(deltaX, omega)
     !< Transforms the anti-hermitian traceless deltaX
     !< Into the adjoint rep for each site
     complex(kind=WC), dimension(:,:,:,:,:,:), intent(in) :: deltaX
     !< Residual field: (3,3,nt,nx,ny,nz)
     real(kind=WP), dimension(:,:,:,:,:), intent(out) :: omega
     !< Adjoint field: (8,nt,nx,ny,nz)
     integer :: it, ix, iy, iz
     !< counters
     do concurrent( it=1:size(deltaX,3), &
          ix=1:size(deltaX,4), &
          iy=1:size(deltaX,5), &
          iz=1:size(deltaX,6)) &
          default(none) shared(omega, deltaX)
        call colourDecompResidual( &
             omega(:,it,ix,iy,iz), &
             deltaX(:,:,it,ix,iy,iz) )
     end do
   end subroutine residualFieldToAdjoint

   pure subroutine adjointFieldToResidual(omega, deltaX)
     !< Does inverse of the below at each site
     !< Decompose traceless anti-Hermitian matrix
     !<
     !<    A = i * sum_a omega(a) * lambda_a
     !<
     !< into the 8 real coefficients omega(a).
     complex(kind=WC), dimension(:,:,:,:,:,:), intent(out) :: deltaX
     !< output traceless anti-hermitian matrices
     real(kind=WP),    dimension(:,:,:,:,:),   intent(in)  :: omega
     !< input real coefficents of adjoin
     integer :: it, ix, iy, iz
     !< counters
     do concurrent(it=1:size(omega,2), &
          ix=1:size(omega,3), &
          iy=1:size(omega,4), &
          iz=1:size(omega,5)) &
          default(none) shared(omega, deltaX)
        call adjointToResidual( &
             deltaX(:,:,it,ix,iy,iz), &
             omega(:,it,ix,iy,iz) )
     end do
   end subroutine adjointFieldToResidual

   pure subroutine projectCoulombZeroModes(omega)
     !< Remove the constant gauge zero modes on every time slice
     real(kind=WP), dimension(:,:,:,:,:), intent(inout) :: omega
     !< 8, nt,nx,ny,nz
     integer :: it, ix, iy, iz
     !< counters
     integer :: nx, ny, nz
     !< lattice dimensions
     real(kind=WP), dimension(size(omega,1)) :: spatialMean
     !< The mean across space sites
     nx = size(omega,3)
     ny = size(omega,4)
     nz = size(omega,5)
     do it = 1, size(omega,2)
        spatialMean = 0.0_WP
        do concurrent(ix=1:nx, iy=1:ny, iz=1:nz) &
             default(none) shared(omega,it) reduce(+:spatialMean)
           spatialMean = spatialMean + omega(:,it,ix,iy,iz)
        end do
        spatialMean = spatialMean / real(nx*ny*nz, kind=WP)
        do concurrent(ix=1:nx, iy=1:ny, iz=1:nz) &
             default(none) shared(omega,it,spatialMean)
           omega(:,it,ix,iy,iz) = omega(:,it,ix,iy,iz) - spatialMean
        end do
     end do
   end subroutine projectCoulombZeroModes

   subroutine solveLaplacianTimeSlice( rhsLinear, solutionLinear, rtol, atol, maxIter, nx, ny, nz, iterations)
     !< Solves the projected spatial Poisson equatioon
     !< $ \nabla^2 \phi = f$
     !< Where $\nabla^2$ is the Laplace operator, $\phi$ is the input residual field
     !< and $f$ is the solution
     !< The solve is on a single timeslice
     !< The projection handles the constant gauge mode
     !< No matrix is explicitly assembled
     !< Uses the stdlib CG
     real(kind=DP), dimension(:), intent(in) :: rhsLinear
     !< 8*nx*ny*nz. Eight adjoint colour components on all spatial sites
     real(kind=DP), dimension(:), intent(inout) :: solutionLinear
     !< 8*nx*ny*nz
     real(kind=WP), intent(in) :: rtol, atol
     !< relative/absolute tolerance for convergence
     integer, intent(in) :: maxIter
     !< max number of iterations
     integer, intent(in) :: nx, ny, nz
     !< lattice spatial dims
     integer, intent(out) :: iterations
     !< number of iterations of CG
     type(stdlib_linop_dp_type) :: lapOp
     !< the linear operator -- i.e. laplacianMatVec
     type(stdlib_solver_workspace_dp_type) :: workspace
     !< workspace for the solver
     integer :: callBackIteration
     !< iterations used - separate variable required for stdlib CG
     lapOp%matvec => laplacianMatVec
     allocate(workspace%tmp(size(rhsLinear), stdlib_size_wksp_cg), &
          source = 0.0_DP)
     callBackIteration = 0
     workspace%callback => recordIteration
     ! Solve the CG
     call stdlib_solve_cg_kernel(&
          lapOp, &
          rhsLinear, &
          solutionLinear, &
          rtol, &
          atol, &
          maxIter, &
          workspace)
     iterations = callBackIteration
     ! Deallocate workspace in type...
     deallocate(workspace%tmp)
   contains
     subroutine recordIteration(x, norm_sq, iter)
       !< stdlib CG callback routine
       !< Records the most recent iteration number reported by solver
       real(kind=DP), dimension(:), intent(in) :: x
       !< solution vector at this stage of CG
       real(kind=DP), intent(in) :: norm_sq
       !< measure of convergence
       integer, intent(in) :: iter
       !< current iteration of CG
       callbackIteration = iter
     end subroutine recordIteration
     subroutine laplacianMatVec(x, y, alpha, beta, op)
       !< Applies the projected 3Dim periodic lattice Laplacian
       !< In adjoint rep
       !< For each colour component a, and spatial lattice site x
       !< $(-\nabla^2\,f)(x) = 6f(x) - \sum_{i=1}^{3}\left[f(x+\hat{i}) + f(x-\hat{i})\right]$
       !< Additionally, it forms
       !< $L_p = -nabla^2 + P_0$
       !< where $P_0\,f = \langle f \rangle$
       !< So that $L_p$ is strictly positive definite.
       real(kind=DP), dimension(:), intent(in) :: x
       !< Input f, 8*nx*ny*nz
       real(kind=DP), dimension(:), intent(inout) :: y
       !< output $ y = \alpha * L_p + \beta * y$, 8*nx*ny*nz
       real(kind=DP), intent(in) :: alpha
       !< $\alpha$
       real(kind=DP), intent(in) :: beta
       !< $\beta$
       character(1), intent(in) :: op
       !< operation type flag. Not really used
       real(kind=DP), dimension(8, nx, ny, nz) :: xfield
       !< reshape of x
       real(kind=DP), dimension(8, nx, ny, nz) :: axfield
       !< $L_p$
       real(kind=DP), dimension(8) :: colourMean
       !< $\langle f \rangle$ for single lattice site
       integer :: ix, iy, iz
       !< lattice site counters
       integer :: ixm, iym, izm
       !< lattice site - 1
       integer :: ixp, iyp, izp
       !< lattice site + 1
       !! TODO replace with spares_op_none/transpose/hermitian from stdlib_sparse
       if (op /= 'N' .and. op /= 'T' .and. op /= 'C') then
          error stop "laplacianMatvec: invalid operator flag: "//op
       end if
       xfield = reshape(x, shape(xfield))
       colourMean = sum(sum(sum(xfield, dim=4), dim=3), dim=2)  ! sum over space
       colourMean = colourMean / real(nx * ny * nz, kind=DP)
       do concurrent (ix = 1:nx, iy = 1:ny, iz = 1:nz) &
            default(none) shared(axfield, xfield, colourMean, nx, ny, nz) &
            local(ixm, ixp, iym, iyp, izm, izp)
          ixm = periodic_minus(ix, nx)
          ixp = periodic_plus(ix, nx)
          iym = periodic_minus(iy, ny)
          iyp = periodic_plus(iy, ny)
          izm = periodic_minus(iz, nz)
          izp = periodic_plus(iz, nz)
          axfield(:, ix, iy, iz) = &
               6.0_DP * xfield(:, ix, iy, iz)       &
               - xfield(:, ixm, iy, iz)             &
               - xfield(:, ixp, iy, iz)             &
               - xfield(:, ix, iym, iz)             &
               - xfield(:, ix, iyp, iz)             &
               - xfield(:, ix, iy, izm)             &
               - xfield(:, ix, iy, izp)             &
               + colourMean
       end do
       y = alpha * reshape(axfield, [size(y)]) + beta * y
     end subroutine laplacianMatVec

     pure function periodic_plus(index, extent) result(next)
       !< Gets the next positive (plus) index periodically
       integer, intent(in) :: index
       !< The current index
       integer, intent(in) :: extent
       !< The total extent
       integer :: next
       !< the next index
       if (index == extent) then
          next = 1
       else
          next = index + 1
       end if
     end function periodic_plus

     pure function periodic_minus(index, extent) result(next)
       !< Gets the next negative (minus) index periodically
       integer, intent(in) :: index
       !< The current index
       integer, intent(in) :: extent
       !< The total extent
       integer :: next
       !< the next index
       if (index == 1) then
          next = extent
       else
          next = index - 1
       end if
     end function periodic_minus

   end subroutine solveLaplacianTimeSlice

end module flue_gaugeFix_Coulomb
