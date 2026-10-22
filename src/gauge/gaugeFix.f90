module FLUE_gaugeFix
  !< Main module for gauge fixing $SU(3)$
  !< Contains interfaces & main subroutines that must be called
  use FLUE_constants, only: WP, WC, assertSTDLIBKindCompatibility
  use FlUE_gaugeFix_Coulomb, only: Coulomb_functional, Coulomb_residual, Coulomb_theta, Coulomb_constructUpdateField
  use FLUE_matrixConstants, only: Ident3x3
  use FLUE_SU3MatrixOps, only: adjointToResidual, expIQ, multiplymatmat, multiplymatmatdag
  use FLUE_wloops, only: periodCoord
  use stdlib_ascii, only: to_lower
  implicit none(type, external)
  private

  abstract interface
     !< Abstract interface for the gauge fixing functional evaluation
     function gaugeFix_functional(U) result(funcVal)
       !< The function interface for the gauge fixing functional evaluation
       import :: WC, WP
       !< complex, real kinds
       implicit none(type, external)
       complex(kind=WC), dimension(:,:,:,:,:,:,:), intent(in) :: U
       !< The gaugefield 3,3,4,nt,nx,ny,nz
       real(kind=WP) :: funcVal
       !< the value of the gauge fixing functional
     end function gaugeFix_functional
  end interface

  abstract interface
     !< Abstract interface for the gauge fixing residual evaluation
     subroutine gaugeFix_residual(U, deltaX)
       !< The subroutine interface for the gauge fixing residual evaulation
       import :: WC
       !< complex kind
       implicit none(type, external)
       complex(kind=WC), dimension(:,:,:,:,:,:,:), intent(in) :: U
       !< The gaugefield 3,3,4,nt,nx,ny,nz
       complex(kind=WC), dimension(:,:,:,:,:,:), intent(out) :: deltaX
       !< The residual 3,3,nt,nx,ny,nz
     end subroutine gaugeFix_residual
  end interface

  abstract interface
     !< Abstract interface for gauge fixing condition/tolerance evaluation
     subroutine gaugeFix_theta(deltaX, thetaParts)
       !< The subroutine interface for gauge fixing condition/tolerance evaluation
       import :: WC, WP
       !< complex, real kinds
       implicit none(type, external)
       complex(kind=WC), dimension(:,:,:,:,:,:), intent(in) :: deltaX
       !< The gauge fixing residual
       real(kind=WP), dimension(:), intent(out) :: thetaParts
       !< the gauge fixing condition, split across timeslices if applicable
     end subroutine gaugeFix_theta
  end interface

  abstract interface
     !< Abstract interface for constructing the update field $\Omega(x)$
     subroutine gaugeFix_constructUpdateField(deltaX, innerRelTol, innerAbsTol, maxInnerIterations, innerIterations, updateField)
       !< The subroutine interface for constructing the update field $\Omega(x)$
       import :: WC, WP
       !< complex, real kinds
       implicit none(type, external)
       complex(kind=WC), dimension(:,:,:,:,:,:), intent(in) :: deltaX
       !< gauge fixing residual 3,3,nt,nx,ny,nz
       real(kind=WP), intent(in) :: innerRelTol
       !< inner iteration relative tolerance
       real(kind=WP), intent(in) :: innerAbsTol
       !< inner iteration absolute tolerance
       integer, intent(in) :: maxInnerIterations
       !< max number of inner iterations (per timeslice if relevant elsse total)
       integer, intent(inout) :: innerIterations
       !< number of inner iterations
       complex(kind=WC), dimension(:,:,:,:,:,:), intent(out) :: updateField
       !< update field $\Omega(x)$ 3,3,nt,nx,ny,nz
     end subroutine gaugeFix_constructUpdateField
  end interface

  type :: gaugeFix_info
     !< Type for information about how the gauge fixing went
     logical :: converged = .false.
     !< convergence status
     integer :: OuterIterations = 0
     !< Number of outer iterations
     integer :: InnerIterations = 0
     !< Number of inner iterations
     real(kind=WP) :: initialTheta = 0.0_WP
     !< Initial value of gauge fixing condition theta
     real(kind=WP) :: finalTheta = 0.0_WP
     !< Final value of gauge fixing condition theta
     real(kind=WP) :: maxFinalThetaPart = 0.0_WP
     !< Final value of maximum of gauge fixing condition theta parts, i.e. over timeslices
     real(kind=WP) :: initialFunctional = 0.0_WP
     !< initial value of gauge fixing functional
     real(kind=WP) :: finalFunctional = 0.0_WP
     !< final value of gauge fixing functional
  end type gaugeFix_info

  type :: gaugeFix_options
     !< Type for gauge fixing options
     real(WP) :: outerTolerance = 1.0e-6_WP
     !< Tolerance of gauge fixing condition theta
     real(WP) :: innerReltol = 1.0e-6_WP
     !< Relative tolerance for inner updates (i.e. CG)
     real(WP) :: innerAbstol = 1.0e-14_WP
     !< Absolute tolerance for inner updates (i.e. CG)
     real(WP) :: initialStep = 1.0_WP
     !< initial $a$ for $g(x) = exp(a * \Omega(x)$
     real(WP) :: minimumStep = 1.0e-8_WP
     !< minimum value of $a$
     integer :: maxOuterIterations = 5000
     !< max number of outer iterations total
     integer :: maxInnerIterations = 1000
     !< max number of inner iteration at each outer step
     integer :: maxLinesearch = 28
     !< max number of steps for $a$ linesearch
     logical :: verbose = .false.
     !< control extra print statements (not currently used)
     character(len=40) :: fixType = 'none'
     !< what kind of gauge fixing to do
  end type gaugeFix_options

  public :: gaugeFix
  public :: applyFiniteGaugeTransformField
  public :: gaugeFix_info, gaugeFix_options

contains

  subroutine gaugeFix(U, options, info, GTrans)
    !< main public gauge fixing subroutine to call
    complex(kind=WC), dimension(:,:,:,:,:,:,:), intent(in) :: U
    !< gaugefield 3,3,4,nt,nx,ny,nz
    type(gaugeFix_options), intent(in) :: options
    !< options for the gauge fixing
    type(gaugeFix_info), intent(out) :: info
    !< info on how the gauge fixing went
    complex(kind=WC), dimension(:,:,:,:,:,:), intent(out) :: GTrans
    !< output gauge transform
    complex(kind=WC), dimension(:,:,:,:,:,:), allocatable :: antiHermitianUpdateField
    !< the anti-hermitian update field $\Omega(x)$
    complex(kind=WC), dimension(:,:,:,:,:,:), allocatable :: GInc
    !< For holding the incremental update to GTrans
    procedure(gaugeFix_functional), pointer :: functional
    !< The gauge fixing functional function
    procedure(gaugeFix_residual), pointer :: residual
    !< the gauge fixing residual subroutine
    procedure(gaugeFix_theta), pointer :: theta
    !< the gauge fixing condition theta subroutine
    procedure(gaugeFix_constructUpdateField), pointer :: constructUpdateField
    !< The subroutine to construct the update field $\Omega(x)$
    integer :: nt, nx, ny, nz
    !< lattice dimensions
    integer :: it, ix, iy, iz, outIter, lineIteration
    !< counters
    real(kind=WP) :: oldFuncVal, newFuncVal, trialFuncVal
    !< the value of the functional before / after
    real(kind=WP), dimension(:), allocatable :: thetaParts
    !< the theta observable
    complex(kind=WC), dimension(:,:,:,:,:,:), allocatable :: deltaX
    !< the residual
    complex(kind=WC), dimension(:,:,:,:,:,:,:), allocatable :: UFixed, UTrial
    !< A copy of the gaugefield where we apply the gauge fix..
    complex(kind=WC), dimension(3,3) :: product
    !< holder for assembling gtrans from incremental
    logical :: accepted
    !< whether the linesearch was accepted
    real(kind=WP) :: alpha
    !< The step in $g(x) = exp(a*\Omega(x)$
    ! check that stdlib DP is WP
    call assertSTDLIBKindCompatibility()
    ! Get lattice dimensions
    nt = size(U, 4)
    nx = size(U, 5)
    ny = size(U, 6)
    nz = size(U, 7)
    select case(trim(to_lower(options%fixType)))
    case('none')
       error stop "You forgot to set the fixType."
    case('coulomb')
       functional => Coulomb_functional
       residual => Coulomb_residual
       theta => Coulomb_theta
       constructUpdateField => Coulomb_constructUpdateField
       ! allocate here cause thetaParts is different for 3d/4d functionals
       allocate(thetaParts(nt))
    case default
       error stop "Bad fixType: "//trim(to_lower(options%fixType))
    end select

    allocate(deltaX(3,3,nt,nx,ny,nz))
    allocate(GInc, mold=GTrans)
    allocate(antiHermitianUpdateField, mold=deltaX)
    allocate(UFixed, source=U)
    allocate(UTrial, source=U)
    ! Get initial values
    oldFuncVal = functional(U)
    call residual(U, deltaX)
    call theta(deltaX, thetaParts)
    info%initialFunctional = oldFuncVal
    info%initialTheta = sum(thetaParts) / real(size(thetaParts), kind=WP)
    ! Set the gauge transform to identity at start
    do concurrent (it = 1:nt, ix = 1:nx, iy = 1:ny, iz = 1:nz)
       GTrans(:, :, it, ix, iy, iz) = Ident3x3
    end do
    ! Short circuit if already gauge fixed
    if (maxval(thetaParts) <= options%outerTolerance) then
       info%converged = .true.
       info%finalFunctional = info%initialFunctional
       info%finalTheta = info%initialTheta
       info%maxFinalThetaPart = maxval(thetaParts)
       return
    end if
    outerLoop: do outIter =1, options%maxOuterIterations
       UTrial = UFixed
       ! This is a repeat calc of residual on iter 1
       call residual(UFixed, deltaX)
       call constructUpdateField(deltaX, &
            options%innerRelTol, &
            options%innerAbsTol, &
            options%maxInnerIterations, &
            info%innerIterations, &
            antiHermitianUpdateField)
       accepted = .false.
       alpha = options%initialStep
       trialFuncVal = oldFuncVal
       linesearchLoop: do lineIteration =1, options%maxLinesearch
          ! write(*,*) 'lineIter', lineIteration, oldFuncVal
          call makeIncrementalTransformation(cmplx(0.0_WP, -1.0_WP) * antiHermitianUpdateField,&
               -alpha, GInc)
          call applyFiniteGaugeTransformField(UTrial, GInc, UFixed)

          trialFuncVal = functional(UFixed)
          if (trialFuncVal > oldFuncVal) then
             accepted = .true.
             exit linesearchLoop
          end if
          alpha = 0.5_WP * alpha
          if (alpha < options%minimumStep) exit linesearchLoop
       end do linesearchLoop
       if (.not. accepted ) then
          exit outerLoop
       end if
       UTrial = UFixed
       do concurrent (it = 1:nt, ix = 1:nx, iy = 1:ny, iz = 1:nz) &
            local(product) shared(gInc, gTrans)
          call multiplymatmat( &
               product,                              &
               gInc(:, :, it, ix, iy, iz),         &
               gTrans(:, :, it, ix, iy, iz))
          gTrans(:, :, it, ix, iy, iz) = product
       end do
       oldFuncVal = trialFuncVal
       call residual(UFixed, deltaX)
       call theta(deltaX, thetaParts)
       info%innerIterations = info%innerIterations + 1
       write(*,*) 'outIter, oldFuncVal, maxval(thetaParts)', outIter, oldFuncVal, maxval(thetaParts)
       if (maxval(thetaParts) <= options%outerTolerance) then
          info%finalTheta = sum(thetaParts) / real(size(thetaParts), kind=WP)
          info%maxFinalThetaPart = maxval(thetaParts)
          info%finalFunctional = oldFuncVal
          info%OuterIterations = outIter
          info%converged = .true.
          exit outerLoop
       end if


    end do outerLoop

  end subroutine gaugeFix


   subroutine applyFiniteGaugeTransform(U, omega, Upert)
     !< Apply the finite gauge transformation
     !<
     !<   Upert_mu(x) = g(x) U_mu(x) g(x+mu)^dagger
     !<
     !< where
     !<
     !<   Omega(x) = i sum_a omega_a(x) lambda_a
     !<   Q(x)     = -i Omega(x)
     !<   g(x)     = exp(i Q(x)) = exp(Omega(x)).
     !<
     !< omega has shape:
     !<
     !<   (8,nt,nx,ny,nz)
     !<
     !< and U/Upert have shape:
     !<
     !<   (3,3,4,nt,nx,ny,nz)
     complex(kind=WC), dimension(:,:,:,:,:,:,:), intent(in) :: U
     !< gaugefield 3,3,4,nt,nx,ny,nz
     real(kind=WP), dimension(:,:,:,:,:), intent(in) :: omega
     !< gauge transform g(x), 8,nt,nx,ny,nz in adjoint
     complex(kind=WC), dimension(:,:,:,:,:,:,:), intent(out) :: Upert
     !< gauge fixed gaugefield 3,3,4,nt,nx,ny,nz
     integer :: mu, it, ix, iy, iz
     !< counters
     integer :: nt, nx, ny, nz
     !< lattice dimensions
     integer, dimension(4) :: coord
     !< current lattice coordinate it,ix,iy,iz
     integer, dimension(4) :: pathCoord
     !< Used to step around the lattice
     integer, dimension(7) :: dataShape
     !< shape of U, 3,3,4,nt,nx,ny,nz

     complex(kind=WC), dimension(3,3) :: OmegaHere
     !< $\omega(x)$ in 3x3 rather than adjoint
     complex(kind=WC), dimension(3,3) :: QHere
     !< Make $\omega(x)$ hermitian
     complex(kind=WC), dimension(3,3) :: GHere
     !< $g(x) = expIQ(QHere)$
     complex(kind=WC), dimension(3,3) :: GNext
     !< $G(x + mu)$
     complex(kind=WC), dimension(3,3) :: work
     !< working matrix for matmul

     nt = size(U,4)
     nx = size(U,5)
     ny = size(U,6)
     nz = size(U,7)
     dataShape = shape(U)
     ! Transform all four link directions. Although the Coulomb
     ! functional and residual involve only mu=2:4, a gauge
     ! transformation acts on the temporal links as well.
     do mu = 1, 4
        pathCoord = 0
        pathCoord(mu) = 1
        do concurrent( &
             it=1:nt, &
             ix=1:nx, &
             iy=1:ny, &
             iz=1:nz) &
             default(none) &
             shared(U,omega,Upert,mu,pathCoord,dataShape) &
             local(coord,OmegaHere,QHere, &
             GHere,GNext,work)
           !-----------------------------------------------------
           ! Construct Omega(x)
           !-----------------------------------------------------
           call adjointToResidual( &
                OmegaHere, &
                omega(:,it,ix,iy,iz))
           ! adjointToResidual returns anti-Hermitian Omega.
           ! ExpIQ requires Hermitian Q, with Omega = i Q.
           QHere = cmplx(0.0_WP, -1.0_WP, kind=WC) * OmegaHere
           GHere = ExpIQ(QHere)
           !-----------------------------------------------------
           ! Locate x + mu using periodic boundary conditions
           !-----------------------------------------------------
           coord = (/it,ix,iy,iz/) + pathCoord
           coord = periodCoord(coord, dataShape)
           !-----------------------------------------------------
           ! Construct Omega(x+mu)
           !-----------------------------------------------------
           call adjointToResidual( &
                OmegaHere, &
                omega(:,coord(1),coord(2),coord(3),coord(4)))
           QHere = cmplx(0.0_WP, -1.0_WP, kind=WC) * OmegaHere
           GNext = ExpIQ(QHere)
           !-----------------------------------------------------
           ! Upert_mu(x) = g(x) U_mu(x) g(x+mu)^dagger
           !-----------------------------------------------------
           call MultiplyMatMat(work, GHere, U(:,:,mu,it,ix,iy,iz))
           call MultiplyMatMatDag(Upert(:,:,mu,it,ix,iy,iz), work, GNext)
        end do
     end do
   end subroutine applyFiniteGaugeTransform

   subroutine applyFiniteGaugeTransformField(U, Gtrans, Upert)
     !< Apply the finite gauge transformation
     !<
     !<   Upert_mu(x) = g(x) U_mu(x) g(x+mu)^dagger
     !< and U/Upert have shape:
     !<
     !<   (3,3,4,nt,nx,ny,nz)
     complex(kind=WC), dimension(:,:,:,:,:,:,:), intent(in) :: U
     !< gaugefield 3,3,4,nt,nx,ny,nz
     complex(kind=WC), dimension(:,:,:,:,:,:), intent(in) :: Gtrans
     !< gauge transform g(x), 3,3,nt,nx,ny,nz
     complex(kind=WC), dimension(:,:,:,:,:,:,:), intent(out) :: Upert
     !< gauge fixed gaugefield 3,3,4,nt,nx,ny,nz
     integer :: mu, it, ix, iy, iz
     !< counters
     integer :: nt, nx, ny, nz
     !< lattice dimensions
     integer, dimension(4) :: coord
     !< current lattice coordinate it,ix,iy,iz
     integer, dimension(4) :: pathCoord
     !< used to step around the lattice
     integer, dimension(7) :: dataShape
     !< shape of U, 3,3,4,nt,nx,ny,nz
     complex(kind=WC), dimension(3,3) :: GHere
     !< g(x)
     complex(kind=WC), dimension(3,3) :: GNext
     !< g(x+mu)
     complex(kind=WC), dimension(3,3) :: work
     !< working matrix for matmul
     nt = size(U,4)
     nx = size(U,5)
     ny = size(U,6)
     nz = size(U,7)
     dataShape = shape(U)
     ! Transform all four link directions. Although the Coulomb
     ! functional and residual involve only mu=2:4, a gauge
     ! transformation acts on the temporal links as well.
     do mu = 1, 4
        pathCoord = 0
        pathCoord(mu) = 1
        do concurrent( &
             it=1:nt, &
             ix=1:nx, &
             iy=1:ny, &
             iz=1:nz) &
             default(none) &
             shared(U,Upert,mu,pathCoord,dataShape,GTrans) &
             local(coord,GHere,GNext,work)
           GHere = GTrans(:,:, it, ix, iy, iz)
           !-----------------------------------------------------
           ! Locate x + mu using periodic boundary conditions
           !-----------------------------------------------------
           coord = (/it,ix,iy,iz/) + pathCoord
           coord = periodCoord(coord, dataShape)
           GNext = GTrans(:,:,coord(1), coord(2), coord(3), coord(4))
           !-----------------------------------------------------
           ! Upert_mu(x) = g(x) U_mu(x) g(x+mu)^dagger
           !-----------------------------------------------------
           call MultiplyMatMat(work, GHere, U(:,:,mu,it,ix,iy,iz))
           call MultiplyMatMatDag(Upert(:,:,mu,it,ix,iy,iz), work, GNext)
        end do
     end do
   end subroutine applyFiniteGaugeTransformField

   subroutine makeIncrementalTransformation(qfield, alpha, gTrans)
     !< Build g(x) = ExpIQ(alpha Q(x))
     !< where Q(x) is Hermititan traceless
     complex(kind=WC), dimension(:,:,:,:,:,:), intent(in) :: qfield
     !< hermitian traceless3,3,nt,nx,ny,nz
     real(kind=WP), intent(in) :: alpha
     !< amount of incremental
     complex(wc), dimension(:, :, :, :, :, :), intent(out) :: gTrans
     !< g(x)3,3, nt,nx,ny,nz
     integer :: nt, nx, ny, nz
     !< lattice dimensions
     integer :: it, ix, iy, iz
     !< counters
     nt = size(qfield, 3)
     nx = size(qfield, 4)
     ny = size(qfield, 5)
     nz = size(qfield, 6)
     if (any(shape(gTrans) /= shape(qfield))) then
        error stop "make_incremental_transformation: shape mismatch"
     end if
     do concurrent (it = 1:nt, ix = 1:nx, iy = 1:ny, iz = 1:nz)
        gTrans(:, :, it, ix, iy, iz) = &
             expIQ(alpha * qfield(:, :, it, ix, iy, iz))
     end do
   end subroutine makeIncrementalTransformation



end module FLUE_gaugeFix
