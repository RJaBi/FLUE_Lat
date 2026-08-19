
module flue_wloops
  !< Module to do calculations for wilson loops and generic paths of gauge links in SU3
  use flue_constants, only: wp, wc, sp
  use flue_su3matrixops, only: multiplymatmat, multiplymatdagmatdag, &
       realtracemultmatmat, tracelessconjgsubtract, colourdecomp, realtracemat
  use flue_matrixconstants, only: ident3x3
  use m_stopwatch, only: watchtype, create_watch, start_watch, stop_watch, destroy_watch, read_watch
  implicit none(external)
  public


contains
  !$omp declare target
  pure function genericpath(data, coordbase, path) result(u_xd)
    !< Returns the links starting from coordbase multiplied together along path
    !< Used to construct wilson loops, staples, etc
    complex(kind=wc), dimension(:, :, :, :, :, :, :), intent(in) :: data
    !< The gaugefield of shape 3,3,4,nt,nx,ny,nz
    integer, dimension(4), intent(in) :: coordbase
    !< The starting coordinate [it, ix, iy, iz]
    integer, dimension(:), intent(in) :: path
    !< The path to follow in terms of +ve/-ve mu/nu/etc
    complex(kind=wc), dimension(3, 3) :: u_xd
    !< The resulting path ordered multiplication of links
    ! internal
    complex(kind=wc), dimension(3, 3) :: amat, bmat
    !< Working SU3 variables
    integer, dimension(4) :: coord, pcoord
    !< coordinate holders for moving around [it, ix, iy, iz]
    integer, dimension(7) :: datashape
    !< The shape of the data [3,3,4,nt,nx,ny,nz]
    integer :: pp
    !< counter for step along path
    integer :: mu
    !< direction to move in
    coord = coordbase
    u_xd = ident3x3
    ! Do this to prevent the array temporary in periodCoord
    datashape = SHAPE(data)
    do pp = 1, SIZE(path)
       bmat = u_xd
       mu = path(pp)
       pcoord = 0
       pcoord(ABS(mu)) = 1
       if (mu < 0) then
          ! the step is backwards so subtract off from coord first
          coord = coord - pcoord
          coord = periodcoord(coord, datashape)
          ! get the link here going forward in mu
          ! and dagger it
          amat = CONJG(TRANSPOSE(data(:, :, ABS(mu), coord(1), coord(2), coord(3), coord(4))))
       else
          amat = data(:, :, mu, coord(1), coord(2), coord(3), coord(4))
          coord = coord + pcoord
          coord = periodcoord(coord, datashape)
       end if
       ! now multiply it into U_xd from the right
       call multiplymatmat(u_xd, bmat, amat)
    end do
  end function genericpath

  subroutine genplaquette(data, nt, nx, ny, nz, mustart, muend, nuend, sumtrp, np, time)
    !< Calculates the plaquette using genplaquette
    !< Outputs the sum over all real trace plaquettes, the number of plaquettes and the time taken
    !< Do all plaquettes using mustart=1, muend=4, nuend=4
    !< Do spatial plaquettes using mustart=2, muend=4, nuend=4
    !< Do temporal plaquettes using mustart=1, muend=1, nuend=4
    complex(kind=wc), dimension(3, 3, 4, nt, nx, ny, nz), intent(in) :: data
    !< The gaugefield
    integer, intent(in) :: mustart, muend, nuend
    !< Which dimensions to loop over
    integer, intent(in) :: nt, nx, ny, nz
    !< The size of the lattice
    real(kind=wp), intent(out) :: sumtrp
    !< Sum of ReTr(P) of all plaquettes
    real(kind=WP), intent(out) :: time
    !< Time taken for calculation
    integer, intent(out) :: np
    !< Number of plaquettes
    integer :: mu, nu, nnx, nny, nnz, nnt
    !< counters
    ! other variables
    complex(kind=wc), dimension(3, 3) :: plaq
    !< plaquette variable holder
    integer, dimension(7) :: datashape
    !< Size of the gaugefield
    integer, dimension(4) :: plaqpath, coordbase
    !< Variables for constructing all plaquettes using genericPath
    real(kind=wp) :: p
    !< real trace of single plaquette
    ! Timers
    type(watchtype) :: watch
    !< holds multiple types of timer (cpu, wall, etc)
    real(kind=sp) :: watchtime
    !< output wall time
    call create_watch(watch)
    call start_watch(watch)
    datashape = (/3, 3, 4, nt, nx, ny, nz/)
    !# hold the sum
    sumtrp = 0.0_WP
    !# hold the number measured
    np = 0
    do mu = mustart, muend
       do nu = mu + 1, nuend
          plaqpath = (/mu, nu, -mu, -nu/)
#ifdef  LOCALITYSUPPORT
          do concurrent(nnx=1:nx, nny=1:ny, nnz=1:nz, nnt=1:nt) &
               local(coordbase, plaq, p) reduce(+:sumtrp, np) shared(data)
#elif OMP
             !$omp parallel do collapse(4) reduction(+:sumTrp,nP) shared(data) private(coordBase, plaq,P)
             do nnx = 1, nx
                do nny = 1, ny
                   do nnz = 1, nz
                      do nnt = 1, nt
#else
                         do nnx = 1, nx
                            do nny = 1, ny
                               do nnz = 1, nz
                                  do nnt = 1, nt
#endif
                                     coordbase = (/nnt, nnx, nny, nnz/)
                                     plaq = genericpath(data, coordbase, plaqpath)

                                     call realtracemultmatmat(p, ident3x3, plaq)
                                     !sumTrP = sumTrP + 1
                                     sumtrp = sumtrp + p
                                     np = np + 1
#ifdef OMP
                                  end do
                               end do
                            end do
                         end do
#elif LOCALITYSUPPORT
                      end do
#else
                   end do
                end do
             end do
          end do
#endif
       end do
    end do
    call stop_watch(watch)
    call read_watch(watchtime, watch, 'wall')
    time = 0.0_WP
    time = watchtime
    !time = end - start
  end subroutine genplaquette

  pure function cloverloopmunucoord(data, nt, nx, ny, nz, mu, nu) result(u_xd)
    !< Calculates the standard clover leaf definition of f_munu
    complex(kind=wc), dimension(3, 3, 4, nt, nx, ny, nz), intent(in) :: data
    !< The gaugefield
    integer, intent(in) :: nt, nx, ny, nz
    !< Dimensions of gaugefield
    integer, intent(in) :: mu, nu
    !< planes to calculate f_munu for
    integer, dimension(4) :: coord
    !< Single coordinate [it, ix, iy, iz]
    complex(kind=wc), dimension(3, 3, nt, nx, ny, nz) :: u_xd
    !< Holds the clover term for lattice site
    complex(kind=wc), dimension(3, 3) :: clovleaf, thisu
    !< working SU3 matrices
    integer, dimension(4, 4) :: plaqpath
    !< The 4 leafs of the clover
    integer, dimension(4) :: coordbase
    !< Working coordinate storage
    complex(kind=wc) :: c12, c13, c23
    !< for projection
    integer :: nnx, nny, nnz, nnt
    !< counters
    ! top left
    plaqpath(1, :) = (/nu, -mu, -nu, mu/)
    ! top right
    plaqpath(2, :) = (/mu, nu, -mu, -nu/)
    ! bottom left
    plaqpath(3, :) = (/-mu, -nu, mu, nu/)
    ! bottom right
    plaqpath(4, :) = (/-nu, mu, nu, -mu/)
#ifdef LOCALITYSUPPORT
    do concurrent(nnx=1:nx, nny=1:ny, nnz=1:nz, nnt=1:nt) &
         default(none) local_init(plaqpath) &
         local(coordbase, clovleaf, thisu, c12, c13, c23) shared(u_xd, data)
#else
       do nnx = 1, nx
          do nny = 1, ny
             do nnz = 1, nz
                do nnt = 1, nt
#endif
                   ! Calculate the clover for this coordinate
                   coordbase = (/nnt, nnx, nny, nnz/)
                   ! 1
                   clovleaf = genericpath(data, coordbase, plaqpath(1, :))
                   thisu = +clovleaf
                   ! 2
                   clovleaf = genericpath(data, coordbase, plaqpath(2, :))
                   thisu = thisu + clovleaf
                   ! 3
                   clovleaf = genericpath(data, coordbase, plaqpath(3, :))
                   thisu = thisu + clovleaf
                   ! 4
                   clovleaf = genericpath(data, coordbase, plaqpath(4, :))
                   thisu = thisu + clovleaf
                   ! Projection & factors
                   ! projection to anti-hermition but not traceless matrix
                   ! as in Borsanyi wilson_flow.c
                   ! a_chm = 0.5*( a - conj(b))
                   c12 = 0.5_WP * (thisu(1, 2) - CONJG(thisu(2, 1)))
                   c13 = 0.5_WP * (thisu(1, 3) - CONJG(thisu(3, 1)))
                   c23 = 0.5_WP * (thisu(2, 3) - CONJG(thisu(3, 2)))
                   thisu(1, 2) = c12
                   thisu(1, 3) = c13
                   thisu(2, 3) = c23
                   thisu(2, 1) = -CONJG(c12)
                   thisu(3, 1) = -CONJG(c13)
                   thisu(3, 2) = -CONJG(c23)
                   thisu(1, 1) = CMPLX(0.0_WP, AIMAG(thisu(1, 1)))
                   thisu(2, 2) = CMPLX(0.0_WP, AIMAG(thisu(2, 2)))
                   thisu(3, 3) = CMPLX(0.0_WP, AIMAG(thisu(3, 3)))
                   u_xd(:, :, nnt, nnx, nny, nnz) = thisu * 0.25_WP
#ifdef LOCALITYSUPPORT
                end do
#else
             end do
          end do
       end do
    end do
#endif
  end function cloverloopmunucoord
#ifdef OMP
  function loop5munucoord(data, nt, nx, ny, nz, mu, nu) result(u_xd)
#else
    pure function loop5munucoord(data, nt, nx, ny, nz, mu, nu) result(u_xd)
#endif
      !< Calculates f_munu according to the 5-loop improved definition. See i.e. hep-lat/0203008
      complex(kind=wc), dimension(3, 3, 4, nt, nx, ny, nz), intent(in) :: data
      !< The gaugefield
      integer, intent(in) :: nt, nx, ny, nz
      !< Dimensions of gaugefield
      integer, intent(in) :: mu, nu
      !< plane to calculate f_munu for
      integer, dimension(4) :: coord
      !< Single coordinate [it, ix, iy, iz]
      complex(kind=wc), dimension(3, 3, nt, nx, ny, nz) :: u_xd
      !< Holds the clover term for each lattice site
      complex(kind=wc), dimension(3, 3) :: clovleaf, thisu, worku
      !< working SU3 Matrices
      integer, dimension(4, 4) :: path1x1
      !< Standard clover leaf paths
      integer, dimension(4, 8) :: path2x2
      !< 2x2 clover leaf paths
      integer, dimension(4, 12) :: path3x3
      !< 3x3 clover leaf paths
      integer, dimension(8, 6) :: path1x2
      !< 1x2 rectangle paths
      integer, dimension(8, 8) :: path1x3
      !< 1x3 rectangle paths
      integer, dimension(4) :: coordbase
      !< Working coordinate storage
      complex(kind=wc) :: c12, c13, c23
      !< for projection
      real(kind=wp), parameter :: k5 = 1.0_WP / 180.0_WP  ! 5 loop
      !< coefficient on 3x3 loops
      !real(kind=WP), parameter :: k5 = 1.0_WP / 90.0_WP   ! 3 loop
      real(kind=wp), parameter :: k1 = (19.0_WP / 9.0_WP) - 55.0_WP * k5
      !< coefficient on 1x1 loops
      real(kind=wp), parameter :: k2 = (1.0_WP / 36.0_WP) - 16.0_WP * k5
      !< coefficient on 2x2 loops
      real(kind=wp), parameter :: k3 = 64.0_WP * k5 - (32.0_WP / 45.0_WP)
      !< coefficient on 1x2 rectangle loops
      real(kind=wp), parameter :: k4 = (1.0_WP / 15.0_WP) - 6.0_WP * k5
      !< coefficient on 1x3 rectangle loops
      !            real(kind=WP), parameter :: k5 = 0.0_WP
      !           real(kind=WP), parameter :: k1 = 0.0_WP
      !          real(kind=WP), parameter :: k2 = 0.0_WP
      !         real(kind=WP), parameter :: k3 = 0.0_WP
      !        real(kind=WP), parameter :: k4 = 1.0_WP
      ! counters
      integer :: nnx, nny, nnz, nnt
      !< counters

      ! Clover 1x1
      ! top left
      path1x1(1, :) = (/nu, -mu, -nu, mu/)
      ! top right
      path1x1(2, :) = (/mu, nu, -mu, -nu/)
      ! bottom left
      path1x1(3, :) = (/-mu, -nu, mu, nu/)
      ! bottom right
      path1x1(4, :) = (/-nu, mu, nu, -mu/)

      ! Clover 2x2
      ! top left
      path2x2(1, :) = (/nu, nu, -mu, -mu, -nu, -nu, mu, mu/)
      ! top right
      path2x2(2, :) = (/mu, mu, nu, nu, -mu, -mu, -nu, -nu/)
      ! bottom left
      path2x2(3, :) = (/-mu, -mu, -nu, -nu, mu, mu, nu, nu/)
      ! bottom right
      path2x2(4, :) = (/-nu, -nu, mu, mu, nu, nu, -mu, -mu/)

      ! Clover 3x3
      ! top left
      path3x3(1, :) = (/nu, nu, nu, -mu, -mu, -mu, -nu, -nu, -nu, mu, mu, mu/)
      ! top right
      path3x3(2, :) = (/mu, mu, mu, nu, nu, nu, -mu, -mu, -mu, -nu, -nu, -nu/)
      ! bottom left
      path3x3(3, :) = (/-mu, -mu, -mu, -nu, -nu, -nu, mu, mu, mu, nu, nu, nu/)
      ! bottom right
      path3x3(4, :) = (/-nu, -nu, -nu, mu, mu, mu, nu, nu, nu, -mu, -mu, -mu/)

      ! Clover 1x2
      ! top left
      path1x2(1, :) = (/nu, nu, -mu, -nu, -nu, mu/)
      path1x2(5, :) = (/nu, -mu, -mu, -nu, mu, mu/)
      ! top right
      path1x2(2, :) = (/mu, nu, nu, -mu, -nu, -nu/)
      path1x2(6, :) = (/mu, mu, nu, -mu, -mu, -nu/)
      ! bottom left
      path1x2(3, :) = (/-mu, -nu, -nu, mu, nu, nu/)
      path1x2(7, :) = (/-mu, -mu, -nu, mu, mu, nu/)
      ! bottom right
      path1x2(4, :) = (/-nu, -nu, mu, nu, nu, -mu/)
      path1x2(8, :) = (/-nu, mu, mu, nu, -mu, -mu/)

      ! Clover 1x3
      ! top left
      path1x3(1, :) = (/nu, nu, nu, -mu, -nu, -nu, -nu, mu/)
      path1x3(5, :) = (/nu, -mu, -mu, -mu, -nu, mu, mu, mu/)
      ! top right
      path1x3(2, :) = (/mu, nu, nu, nu, -mu, -nu, -nu, -nu/)
      path1x3(6, :) = (/mu, mu, mu, nu, -mu, -mu, -mu, -nu/)
      ! bottom left
      path1x3(3, :) = (/-mu, -nu, -nu, -nu, mu, nu, nu, nu/)
      path1x3(7, :) = (/-mu, -mu, -mu, -nu, mu, mu, mu, nu/)
      ! bottom right
      path1x3(4, :) = (/-nu, -nu, -nu, mu, nu, nu, nu, -mu/)
      path1x3(8, :) = (/-nu, mu, mu, mu, nu, -mu, -mu, -mu/)
      u_xd = 0.0
#ifdef LOCALITYSUPPORT
      do concurrent(nnx=1:nx, nny=1:ny, nnz=1:nz, nnt=1:nt) &
           default(none) &
           local_init(path1x1, path1x2, path2x2, path1x3, path3x3) &
           local(coordbase, clovleaf, thisu, c12, c13, c23, worku) shared(u_xd, data)
#elif OMP
         !$omp parallel do collapse(4) shared(data,U_xd) default(none)&
         !$ompprivate(coordBase, clovLeaf, thisU, c12, c13, c23, workU) &
         !$omp firstprivate(path1x1, path1x2, path2x2, path1x3, path3x3, &
         !$omp nt, nx, ny, nz)
         do nnx = 1, nx
            do nny = 1, ny
               do nnz = 1, nz
                  do nnt = 1, nt
#else
                     do nnx = 1, nx
                        do nny = 1, ny
                           do nnz = 1, nz
                              do nnt = 1, nt
#endif

                                 ! Calculate the clover 1x1 for this coordinate
                                 coordbase = (/nnt, nnx, nny, nnz/)
                                 ! 1
                                 clovleaf = genericpath(data, coordbase, path1x1(1, :))
                                 thisu = +clovleaf
                                 ! 2
                                 clovleaf = genericpath(data, coordbase, path1x1(2, :))
                                 thisu = thisu + clovleaf
                                 ! 3
                                 clovleaf = genericpath(data, coordbase, path1x1(3, :))
                                 thisu = thisu + clovleaf
                                 ! 4
                                 clovleaf = genericpath(data, coordbase, path1x1(4, :))
                                 thisu = thisu + clovleaf
                                 thisu = thisu * k1
                                 ! 2x2
                                 ! 1
                                 clovleaf = genericpath(data, coordbase, path2x2(1, :))
                                 worku = +clovleaf
                                 ! 2
                                 clovleaf = genericpath(data, coordbase, path2x2(2, :))
                                 worku = worku + clovleaf
                                 ! 3
                                 clovleaf = genericpath(data, coordbase, path2x2(3, :))
                                 worku = worku + clovleaf
                                 ! 4
                                 clovleaf = genericpath(data, coordbase, path2x2(4, :))
                                 worku = worku + clovleaf
                                 thisu = thisu + worku * k2

!!!!!
!!!!! 3x3
!!!!!
                                 ! 1
                                 clovleaf = genericpath(data, coordbase, path3x3(1, :))
                                 worku = +clovleaf
                                 ! 2
                                 clovleaf = genericpath(data, coordbase, path3x3(2, :))
                                 worku = worku + clovleaf
                                 ! 3
                                 clovleaf = genericpath(data, coordbase, path3x3(3, :))
                                 worku = worku + clovleaf
                                 ! 4
                                 clovleaf = genericpath(data, coordbase, path3x3(4, :))
                                 worku = worku + clovleaf
                                 thisu = thisu + worku * k5

                                 ! 1x2
                                 ! 1
                                 clovleaf = genericpath(data, coordbase, path1x2(1, :))
                                 worku = +clovleaf
                                 clovleaf = genericpath(data, coordbase, path1x2(5, :))
                                 worku = worku + clovleaf
                                 ! 2
                                 clovleaf = genericpath(data, coordbase, path1x2(2, :))
                                 worku = worku + clovleaf
                                 clovleaf = genericpath(data, coordbase, path1x2(6, :))
                                 worku = worku + clovleaf
                                 ! 3
                                 clovleaf = genericpath(data, coordbase, path1x2(3, :))
                                 worku = worku + clovleaf
                                 clovleaf = genericpath(data, coordbase, path1x2(7, :))
                                 worku = worku + clovleaf
                                 ! 4
                                 clovleaf = genericpath(data, coordbase, path1x2(4, :))
                                 worku = worku + clovleaf
                                 clovleaf = genericpath(data, coordbase, path1x2(8, :))

                                 worku = worku + clovleaf
                                 thisu = thisu + worku * k3 * 0.5_WP

                                 ! 1x3
                                 ! 1
                                 clovleaf = genericpath(data, coordbase, path1x3(1, :))
                                 worku = +clovleaf
                                 clovleaf = genericpath(data, coordbase, path1x3(5, :))
                                 worku = worku + clovleaf
                                 ! 2
                                 clovleaf = genericpath(data, coordbase, path1x3(2, :))
                                 worku = worku + clovleaf
                                 clovleaf = genericpath(data, coordbase, path1x3(6, :))
                                 worku = worku + clovleaf
                                 ! 3
                                 clovleaf = genericpath(data, coordbase, path1x3(3, :))
                                 worku = worku + clovleaf
                                 clovleaf = genericpath(data, coordbase, path1x3(7, :))
                                 worku = worku + clovleaf
                                 ! 4
                                 clovleaf = genericpath(data, coordbase, path1x3(4, :))
                                 worku = worku + clovleaf
                                 clovleaf = genericpath(data, coordbase, path1x3(8, :))
                                 worku = worku + clovleaf
                                 thisu = thisu + worku * k4 * 0.5_WP

                                 ! Projection & factors
                                 ! projection to anti-hermition but not traceless matrix
                                 ! as in Borsanyi wilson_flow.c
                                 ! a_chm = 0.5*( a - conj(b))
                                 c12 = 0.5_WP * (thisu(1, 2) - CONJG(thisu(2, 1)))
                                 c13 = 0.5_WP * (thisu(1, 3) - CONJG(thisu(3, 1)))
                                 c23 = 0.5_WP * (thisu(2, 3) - CONJG(thisu(3, 2)))
                                 thisu(1, 2) = c12
                                 thisu(1, 3) = c13
                                 thisu(2, 3) = c23
                                 thisu(2, 1) = -CONJG(c12)
                                 thisu(3, 1) = -CONJG(c13)
                                 thisu(3, 2) = -CONJG(c23)
                                 thisu(1, 1) = CMPLX(0.0_WP, AIMAG(thisu(1, 1)))
                                 thisu(2, 2) = CMPLX(0.0_WP, AIMAG(thisu(2, 2)))
                                 thisu(3, 3) = CMPLX(0.0_WP, AIMAG(thisu(3, 3)))
                                 u_xd(:, :, nnt, nnx, nny, nnz) = thisu * 0.25_WP
#ifdef LOCALITYSUPPORT
                              end do
#else
                           end do
                        end do
                     end do
                  end do
#endif
                end function loop5munucoord

                function plaquettemunucoord(data, nt, nx, ny, nz, mu, nu) result(u_xd)
                  !< Calculates f_munu from the plaquette only
                  complex(kind=wc), dimension(3, 3, 4, nt, nx, ny, nz), intent(in) :: data
                  !< The gaugefield
                  integer, intent(in) :: nt, nx, ny, nz
                  !< Dimensions of gaugefield
                  integer, intent(in) :: mu, nu
                  !< plane to calculate f_munu for
                  integer, dimension(4) :: coord
                  !< Single coordiinate [it, ix, iy, iz]
                  complex(kind=wc), dimension(3, 3, nt, nx, ny, nz) :: u_xd
                  !< Holds the clover term for each lattice site
                  complex(kind=wc), dimension(3, 3) :: clovleaf
                  !< Working SU3 matrix
                  integer, dimension(4) :: plaqpath
                  !< Path for plaquette
                  integer, dimension(4) :: coordbase
                  !< Working coordinate storage
                  complex(kind=wc) :: c12, c13, c23
                  !< for projection
                  integer :: nnx, nny, nnz, nnt
                  !< counters
                  ! top left
                  ! Untested
                  plaqpath = (/nu, -mu, -nu, mu/)
#ifdef LOCALITYSUPPORT
                  do concurrent(nnx=1:nx, nny=1:ny, nnz=1:nz, nnt=1:nt) &
                       default(none) local_init(plaqpath) &
                       local(coordbase, clovleaf, c12, c13, c23) shared(data, u_xd)
#else
                     do concurrent(nnx=1:nx, nny=1:ny, nnz=1:nz, nnt=1:nt)
#endif
                        ! Calculate the clover for this coordinate
                        coordbase = (/nnt, nnx, nny, nnz/)
                        clovLeaf = genericPath(data, coordBase, plaqPath)
                        ! Projection & factors
                        ! projection to anti-hermition but not traceless matrix
                        ! as in Borsanyi wilson_flow.c
                        ! a_chm = 0.5*( a - conj(b))
                        c12 = 0.5_WP * (clovleaf(1, 2) - CONJG(clovleaf(2, 1)))
                        c13 = 0.5_WP * (clovleaf(1, 3) - CONJG(clovleaf(3, 1)))
                        c23 = 0.5_WP * (clovleaf(2, 3) - CONJG(clovleaf(3, 2)))
                        clovleaf(1, 2) = c12
                        clovleaf(1, 3) = c13
                        clovleaf(2, 3) = c23
                        clovleaf(2, 1) = -CONJG(c12)
                        clovleaf(3, 1) = -CONJG(c13)
                        clovleaf(3, 2) = -CONJG(c23)
                        clovleaf(1, 1) = CMPLX(0.0_WP, AIMAG(clovleaf(1, 1)), kind=WC)
                        clovleaf(2, 2) = CMPLX(0.0_WP, AIMAG(clovleaf(2, 2)), kind=WC)
                        clovleaf(3, 3) = CMPLX(0.0_WP, AIMAG(clovleaf(3, 3)), kind=WC)
                        u_xd(:, :, nnt, nnx, nny, nnz) = clovleaf * 0.25_WP
                     end do
                   end function plaquettemunucoord

                   function magnetic(data, nt, nx, ny, nz) result(b)
                     !< Computes the magnetic contribution to the gauge-field action density.
                     !<
                     !< This routine evaluates the magnetic component of the field-strength
                     !< tensor \(F_{\mu\nu}\) using a five-loop clover discretization. The
                     !< three independent spatial components,
                     !<
                     !<
                     !< $B^2 = -\sum_{i=1}^{3} \mathrm{Tr}\left(F_{jk}F_{jk}\right)$,
                     !<
                     !<
                     !< are constructed from the clover operators returned by
                     !< `loop5munucoord`. The result is averaged over all lattice sites.
                     !< The magnetic field components are therefore constructed from
                     !< the spatial field-strength tensors:
                     !<
                     !< - \(F_{23}\)  → `bfield(:,:,:,:,:,:,1)`
                     !< - \(F_{31}\)  → `bfield(:,:,:,:,:,:,2)`
                     !< - \(F_{12}\)  → `bfield(:,:,:,:,:,:,3)`
                     !<
                     !< For each lattice site and magnetic component, the quantity
                     !<
                     !<
                     !<$ -\mathrm{Tr}(F_{ij}F_{ij})$
                     !<
                     !<
                     !< is accumulated and averaged over the lattice volume.
                     complex(kind=wc), dimension(3, 3, 4, nt, nx, ny, nz), intent(in) :: data
                     !< The gaugefield
                     integer, intent(in) :: nt, nx, ny, nz
                     !< Dimensions of gaugefield
                     real(kind=wp) :: b
                     !< magnetic contribution
                     complex(kind=wc), dimension(3, 3, nt, nx, ny, nz, 3) :: bfield
                     !< spatial field-strength tensors in each of 3 planes at each site
                     complex(kind=wc), dimension(3, 3) :: tempeval
                     !< working SU3 matrix
                     real(kind=wp) :: cloverplaq, temppl
                     !< Working variables for after tracing
                     integer :: nnx, nny, nnz, nnt, ii
                     !< counters

                     ! The indices are so that they match the wilson_flow.c values exactly
                     ! there t,x,y,z = 0,3,2,1
                     ! here t,x,y,z = 1,2,3,4
                     ! match 23
                     bfield(:, :, :, :, :, :, 1) = loop5munucoord(data, nt, nx, ny, nz, 3, 2)
                     ! match 31
                     bfield(:, :, :, :, :, :, 2) = loop5munucoord(data, nt, nx, ny, nz, 2, 4)
                     ! match 12
                     bfield(:, :, :, :, :, :, 3) = loop5munucoord(data, nt, nx, ny, nz, 4, 3)
                     cloverplaq = 0.0_WP
#ifdef LOCALITYSUPPORT
                     do concurrent(nny=1:ny, nnz=1:nz, nnt=1:nt, nnx=1:nx, ii=1:3) &
                          default(none) shared(bfield) &
                          local(tempeval, temppl) &
                          reduce(+:cloverplaq)
#elif OMP
                        !$omp parallel do collapse(5) default(none) &
                        !$omp reduction(+:cloverPlaq) shared(BField) &
                        !$omp firstprivate(nx,ny,nz,nt) private(tempPL)
                        do nnx = 1, nx
                           do nny = 1, ny
                              do nnz = 1, nz
                                 do nnt = 1, nt
                                    do ii = 1, 3
#else
                                       do nnx = 1, nx
                                          do nny = 1, ny
                                             do nnz = 1, nz
                                                do nnt = 1, nt
                                                   do ii = 1, 3

#endif
                                                      call realtracemultmatmat(temppl, bfield(:,:, nnt, nnx, nny, nnz, ii), &
                                                           bfield(:, :, nnt, nnx, nny, &
                                                           nnz, ii))
                                                      cloverplaq = cloverplaq - temppl
#ifdef LOCALITYSUPPORT
                                                   end do
#else
                                                end do
                                             end do
                                          end do
                                       end do
                                    end do
#endif
                                    b = cloverplaq / real(nt * nx * ny * nz, kind=wp)
                                  end function magnetic

                                  pure function periodcoord(coord, datashape)
                                    !< Handles periodic boundary conditions
                                    !< Only handles steps of 1!
                                    integer, dimension(7), intent(in) :: datashape
                                    !< Dimensions of gaugefield [colour, colour, mu, NT, NX, NY, NZ]
                                    integer, dimension(4), intent(in) :: coord
                                    !< Coord to check [it, ix, iy, iz]
                                    integer, dimension(4) :: periodcoord
                                    !< Periodically wrapped coordinate
                                    integer :: cc
                                    !< counter
                                    ! A lazy function to handle the periodic boundary conditions
                                    ! dataShape is (colour, colour, mu, NT, nx, ny, nz)
                                    ! coord is (nt, nx, ny, nz)
                                    ! checks if the value in coord is greater than corresponding N in datashape
                                    ! if so sets it to 1
                                    ! i.e. only handles steps of 1
                                    periodcoord = coord
                                    do cc = 1, SIZE(coord)
                                       if (coord(cc) > datashape(cc + 3)) then
                                          periodcoord(cc) = 1
                                       else if (coord(cc) == 0) then
                                          periodcoord(cc) = datashape(cc + 3)
                                       end if
                                    end do
                                  end function periodcoord

                                  pure function wmunutau(data, tau, mu, nu, nmu, nnu) result(loopval)
                                    !< Old unused function to calculate
                                    !< The nMuxNu loop at t=tau in plane mu nu
                                    !< Returns average value of realtrace of loop
                                    !< DOES NOT USE GENERIC PATH
                                    complex(kind=wc), dimension(:, :, :, :, :, :, :), intent(in) :: data
                                    !< The gaugefield [colour, colour, mu, NT, NX, NY, NZ]
                                    integer, intent(in) :: tau
                                    !< time slice to calculate at
                                    integer, intent(in) :: mu, nu
                                    !< plane to calcluate in
                                    integer, intent(in) :: nmu, nnu
                                    !< Dimensions of loop to calculate
                                    real(kind=wp) :: loopval
                                    !< real trace result
                                    !logical(kind=c_bool),                                intent(in) :: verbose
                                    !"""
                                    !Calculates the nMuxnNu loop at t = tau
                                    ! Nmu link in mu direction
                                    ! nNu links in nu direction
                                    !data is [colour, colour, mu, nt, nx, ny, nz] complex
                                    ! Returns the average value of the loop
                                    !"""
                                    integer, dimension(7) :: datashape
                                    !< Dimensions of gaugefield [colour, colour, mu, nt, nx, nz, nz]
                                    integer, dimension(4) :: mucoord, nucoord, coordbase, coord
                                    !< working coordinates [it, ix, iy, iz]
                                    integer :: nx, ny, nz, nn
                                    !< counters
                                    complex(kind=wc), dimension(3, 3) :: umu_x, unu_xpamuh, umu_xpanuh, unu_x
                                    !< Intermediate SU3 matrix
                                    complex(kind=wc), dimension(3, 3) :: utemp, utemp2
                                    !< Intermediate SU3 matrix
                                    complex(kind=wc), dimension(3, 3) :: umuunu, umudagunudag
                                    !< Intermediate SU3 matrix
                                    real(kind=wp) :: p
                                    !< Working sum
                                    datashape = SHAPE(data)
                                    !# hold the sum
                                    loopval = 0.0_WP

                                    mucoord(:) = 0
                                    nucoord(:) = 0
                                    !# This is a single shift in mu
                                    mucoord(mu) = 1
                                    !# This is a single shift in nu
                                    nucoord(nu) = 1
                                    !# loop over all sites
                                    do nx = 1, datashape(5)
                                       do ny = 1, datashape(6)
                                          do nz = 1, datashape(7)
                                             ! The starting coordinate
                                             coordbase = (/tau, nx, ny, nz/)
                                             ! Do U_nMu(x)
                                             coord = coordbase
                                             umu_x = ident3x3
                                             utemp2 = ident3x3
                                             do nn = 1, nmu
                                                ! Get data
                                                utemp = data(:, :, mu, coord(1), coord(2), coord(3), coord(4))
                                                ! Multiply it in
                                                call multiplymatmat(utemp2, umu_x, utemp)
                                                ! and re-assign
                                                umu_x = utemp2
                                                ! Update the coordinate
                                                coord = periodcoord(coord + mucoord, datashape)
                                             end do
                                             ! Do U_nNu(x+nMu*amu)
                                             coord = periodcoord(coordbase + mucoord, datashape)
                                             do nn = 2, nmu
                                                coord = periodcoord(coord + mucoord, datashape)
                                             end do
                                             unu_xpamuh = ident3x3
                                             utemp2 = ident3x3
                                             do nn = 1, nnu
                                                ! Get data
                                                utemp = data(:, :, nu, coord(1), coord(2), coord(3), coord(4))
                                                ! Multiply it in
                                                call multiplymatmat(utemp2, unu_xpamuh, utemp)
                                                ! and re-assign
                                                unu_xpamuh = utemp2
                                                ! Update the coordinate
                                                coord = periodcoord(coord + nucoord, datashape)
                                             end do
                                             ! U_nMu(x+nNu)
                                             coord = periodcoord(coordbase + nucoord, datashape)
                                             do nn = 2, nnu
                                                coord = periodcoord(coord + nucoord, datashape)
                                             end do
                                             umu_xpanuh = ident3x3
                                             utemp2 = ident3x3
                                             do nn = 1, nmu
                                                ! Get data
                                                utemp = data(:, :, mu, coord(1), coord(2), coord(3), coord(4))
                                                ! Multiply it in
                                                call multiplymatmat(utemp2, umu_xpanuh, utemp)
                                                ! and re-assign
                                                umu_xpanuh = utemp2
                                                ! Update the coordinate
                                                coord = periodcoord(coord + mucoord, datashape)
                                             end do
                                             ! U_nNu(x)
                                             coord = coordbase
                                             unu_x = ident3x3
                                             utemp2 = ident3x3
                                             do nn = 1, nnu
                                                ! Get data
                                                utemp = data(:, :, nu, coord(1), coord(2), coord(3), coord(4))
                                                ! Multiply it in
                                                call multiplymatmat(utemp2, unu_x, utemp)
                                                ! and re-assign
                                                unu_x = utemp2
                                                ! Update the coordinate
                                                coord = periodcoord(coord + nucoord, datashape)
                                             end do
                                             !# Multiply bottom, right together
                                             call multiplymatmat(umuunu, umu_x, unu_xpamuh)
                                             !# Multiply left, top together, take dagger
                                             call multiplymatdagmatdag(umudagunudag, umu_xpanuh, unu_x)
                                             !# multiply two halves together, take trace
                                             call realtracemultmatmat(p, umuunu, umudagunudag)
                                             loopval = loopval + p
                                          end do
                                       end do
                                    end do
                                    loopval = loopval / real(datashape(5) * datashape(6) * datashape(7), kind=wp)
                                  end function wmunutau

                                  pure function wmunu(data, mu, nu, nmu, nnu) result(loopval)
                                    !< Calculates the nMuxnNu loop
                                    !< averages over all sites in the mu-nu plane
                                    !< Does not use genericPath
                                    complex(kind=wc), dimension(:, :, :, :, :, :, :), intent(in) :: data
                                    !< The Gaugefield [colour, colour, mu, NT, NX, NY, NZ]
                                    integer, intent(in) :: mu, nu
                                    !< Plane to calculate in
                                    integer, intent(in) :: nmu, nnu
                                    !< Size of loop in mu, nu directions
                                    real(kind=wp) :: loopval
                                    !< Average real trace of loop
                                    !logical(kind=c_bool),                                intent(in) :: verbose
                                    !"""
                                    !Calculates the nMuxnNu loop
                                    ! Nmu link in mu direction
                                    ! nNu links in nu direction
                                    !data is [colour, colour, mu, nt, nx, ny, nz,] complex
                                    ! Returns the average value of the loop
                                    !"""
                                    integer, dimension(7) :: datashape
                                    !< dimensions of gaugefield [colour, colour, mu, NT, NX, NY, NZ]
                                    integer, dimension(4) :: mucoord, nucoord, coordbase, coord
                                    !< Working coordinates [it, ix, iy, iz]
                                    integer :: nx, ny, nz, nt, nn
                                    !< counters
                                    complex(kind=wc), dimension(3, 3) :: umu_x, unu_xpamuh, umu_xpanuh, unu_x
                                    !< Intermediate SU3 matrices
                                    complex(kind=wc), dimension(3, 3) :: utemp, utemp2
                                    !< Intermediate SU3 matrices
                                    complex(kind=wc), dimension(3, 3) :: umuunu, umudagunudag
                                    !< Intermediate SU3 matrices
                                    real(kind=wp) :: p
                                    !< Per site value
                                    datashape = SHAPE(data)
                                    !# hold the sum
                                    loopval = 0.0_WP

                                    mucoord(:) = 0
                                    nucoord(:) = 0
                                    !# This is a single shift in mu
                                    mucoord(mu) = 1
                                    !# This is a single shift in nu
                                    nucoord(nu) = 1
                                    !# loop over all sites
                                    do nx = 1, datashape(5)
                                       do ny = 1, datashape(6)
                                          do nz = 1, datashape(7)
                                             do nt = 1, datashape(4)
                                                ! The starting coordinate
                                                coordbase = (/nt, nx, ny, nz/)
                                                ! Do U_nMu(x)
                                                coord = coordbase
                                                umu_x = ident3x3
                                                utemp2 = ident3x3
                                                do nn = 1, nmu
                                                   ! Get data
                                                   utemp = data(:, :, mu, coord(1), coord(2), coord(3), coord(4))
                                                   ! Multiply it in
                                                   call multiplymatmat(utemp2, umu_x, utemp)
                                                   ! and re-assign
                                                   umu_x = utemp2
                                                   ! Update the coordinate
                                                   coord = periodcoord(coord + mucoord, datashape)
                                                end do
                                                ! Do U_nNu(x+nMu*amu)
                                                coord = periodcoord(coordbase + mucoord, datashape)
                                                do nn = 2, nmu
                                                   coord = periodcoord(coord + mucoord, datashape)
                                                end do
                                                unu_xpamuh = ident3x3
                                                utemp2 = ident3x3
                                                do nn = 1, nnu
                                                   ! Get data
                                                   utemp = data(:, :, nu, coord(1), coord(2), coord(3), coord(4))
                                                   ! Multiply it in
                                                   call multiplymatmat(utemp2, unu_xpamuh, utemp)
                                                   ! and re-assign
                                                   unu_xpamuh = utemp2
                                                   ! Update the coordinate
                                                   coord = periodcoord(coord + nucoord, datashape)
                                                end do
                                                ! U_nMu(x+nNu)
                                                coord = periodcoord(coordbase + nucoord, datashape)
                                                do nn = 2, nnu
                                                   coord = periodcoord(coord + nucoord, datashape)
                                                end do
                                                umu_xpanuh = ident3x3
                                                utemp2 = ident3x3
                                                do nn = 1, nmu
                                                   ! Get data
                                                   utemp = data(:, :, mu, coord(1), coord(2), coord(3), coord(4))
                                                   ! Multiply it in
                                                   call multiplymatmat(utemp2, umu_xpanuh, utemp)
                                                   ! and re-assign
                                                   umu_xpanuh = utemp2
                                                   ! Update the coordinate
                                                   coord = periodcoord(coord + mucoord, datashape)
                                                end do
                                                ! U_nNu(x)
                                                coord = coordbase
                                                unu_x = ident3x3
                                                utemp2 = ident3x3
                                                do nn = 1, nnu
                                                   ! Get data
                                                   utemp = data(:, :, nu, coord(1), coord(2), coord(3), coord(4))
                                                   ! Multiply it in
                                                   call multiplymatmat(utemp2, unu_x, utemp)
                                                   ! and re-assign
                                                   unu_x = utemp2
                                                   ! Update the coordinate
                                                   coord = periodcoord(coord + nucoord, datashape)
                                                end do
                                                !# Multiply bottom, right together
                                                call multiplymatmat(umuunu, umu_x, unu_xpamuh)
                                                !# Multiply left, top together, take dagger
                                                call multiplymatdagmatdag(umudagunudag, umu_xpanuh, unu_x)
                                                !# multiply two halves together, take trace
                                                call realtracemultmatmat(p, umuunu, umudagunudag)
                                                loopval = loopval + p
                                             end do
                                          end do
                                       end do
                                    end do
                                    loopval = loopval / real(SIZE(data) / 36, kind=wp)
                                  end function wmunu

                                  subroutine polyakov(data, sumtrp, np, time)
                                    !< Calculates the Polyakov loop
                                    !< i.e. the real trace of temporal links
                                    !< Does not use genericPath
                                    complex(kind=wc), dimension(:, :, :, :, :, :, :), intent(in) :: data
                                    !< The Gaugefield [colour, colour, mu, NT, NX, NY, NZ]
                                    real(kind=wp), intent(out) :: sumtrp
                                    !< The sum of all polyakov loops across lattice sites
                                    real(kind=WP), intent(out) :: time
                                    !< time taken to calculate
                                    integer, intent(out) :: np
                                    !< Number of loops measured
                                    !"""
                                    !Calculates the polyakov loop
                                    !data is [colour, colour, mu, nt, nx, ny, nz] complex
                                    !returns the sum of polyakov, number of loops measured,
                                    !the average polyakov and the time taken to calculate it
                                    !"""
                                    integer, dimension(7) :: datashape
                                    !< Dimensions of the gaugefield [colour, colour, mu, NT, NX, NY, NZ]
                                    integer, dimension(4) :: coord
                                    !< Working coordinate holder [it, ix, iy, iz]
                                    integer :: nx, ny, nz, nt
                                    !< counters
                                    complex(kind=wc), dimension(3, 3) :: unu_x, unu_xt, unu_temp
                                    !< Working SU3 matrices
                                    real(kind=wp) :: p
                                    !< single site value of polyakov loop
                                    ! Timers
                                    real(kind=wp) :: start, end
                                    !< start and end CPU times
                                    call CPU_TIME(start)
                                    datashape = SHAPE(data)
                                    !# hold the sum
                                    sumtrp = 0.0_WP
                                    !# hold the number measured
                                    np = 0
                                    !# loop over all sites
                                    do nx = 1, datashape(5)
                                       do ny = 1, datashape(6)
                                          do nz = 1, datashape(7)
                                             ! Get first link U_t(x, 1)
                                             coord = (/1, nx, ny, nz/)
                                             unu_x = data(:, :, 1, coord(1), coord(2), coord(3), coord(4))
                                             do nt = 2, datashape(4) - 1
                                                ! Get middle links U_t(x, t)
                                                coord = (/nt, nx, ny, nz/)
                                                unu_xt = data(:, :, 1, coord(1), coord(2), coord(3), coord(4))
                                                call multiplymatmat(unu_temp, unu_x, unu_xt)
                                                unu_x = unu_temp
                                             end do
                                             ! get final link U_t(x, NT)
                                             coord = (/datashape(4), nx, ny, nz/)
                                             unu_x = data(:, :, 1, coord(1), coord(2), coord(3), coord(4))
                                             ! Multiply and trace
                                             call realtracemultmatmat(p, unu_temp, unu_x)
                                             call multiplymatmat(unu_xt, unu_temp, unu_x)
                                             sumtrp = sumtrp + p
                                             np = np + 1
                                          end do
                                       end do
                                    end do
                                    call CPU_TIME(end)
                                    time = end - start
                                  end subroutine polyakov

                                end module flue_wloops
