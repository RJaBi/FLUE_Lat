!! Calculate, spatial, temporal plaquettes
!! Compiles for use with python using
!! 'f2py -c fortPlaq.f90 -m fortPlaq'
!!
!! Which can then be imported into python as
!! import fortPlaq
!! fPlaq = fortPlaq.plaq.plaquette
!! which would give you the plaquette function

module flue_wloops
   use flue_constants, only: wp, wc, sp
   use flue_su3matrixops, only: multiplymatmat, multiplymatdagmatdag, &
                                realtracemultmatmat, tracelessconjgsubtract, colourdecomp, realtracemat
   use flue_matrixconstants, only: ident3x3
   use m_stopwatch, only: watchtype, create_watch, start_watch, stop_watch, destroy_watch, read_watch
   implicit none(external)
   public

contains

   pure function genericpath(data, coordbase, path) result(u_xd)
      complex(kind=wc), dimension(:, :, :, :, :, :, :), intent(in) :: data
      integer, dimension(4), intent(in) :: coordbase
      integer, dimension(:), intent(in) :: path
      complex(kind=wc), dimension(3, 3) :: u_xd
      ! internal
      complex(kind=wc), dimension(3, 3) :: amat, bmat
      integer, dimension(4) :: coord, pcoord
      integer, dimension(7) :: datashape
      ! counters
      integer :: pp
      integer :: mu
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
      complex(kind=wc), dimension(3, 3, 4, nt, nx, ny, nz), intent(in) :: data
      integer, intent(in) :: mustart, muend, nuend
      integer, intent(in) :: nt, nx, ny, nz
      real(kind=wp), intent(out) :: sumtrp, time
      integer, intent(out) :: np
      ! Counters
      integer :: mu, nu, nnx, nny, nnz, nnt
      ! other variables
      complex(kind=wc), dimension(3, 3) :: plaq
      integer, dimension(7) :: datashape
      integer, dimension(4) :: plaqpath, coordbase
      real(kind=wp) :: p
      ! Timers
      type(watchtype) :: watch
      real(kind=sp) :: watchtime
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
      complex(kind=wc), dimension(3, 3, 4, nt, nx, ny, nz), intent(in) :: data
      integer, intent(in) :: nt, nx, ny, nz, mu, nu
      integer, dimension(4) :: coord
      complex(kind=wc), dimension(3, 3, nt, nx, ny, nz) :: u_xd
      complex(kind=wc), dimension(3, 3) :: clovleaf, thisu
      integer, dimension(4, 4) :: plaqpath
      integer, dimension(4) :: coordbase
      complex(kind=wc) :: c12, c13, c23
      ! counters
      integer :: nnx, nny, nnz, nnt

      integer :: t1, t2
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
         ! See hep-lat/0203008
         complex(kind=wc), dimension(3, 3, 4, nt, nx, ny, nz), intent(in) :: data
         integer, intent(in) :: nt, nx, ny, nz, mu, nu
         integer, dimension(4) :: coord
         complex(kind=wc), dimension(3, 3, nt, nx, ny, nz) :: u_xd
         complex(kind=wc), dimension(3, 3) :: clovleaf, thisu, worku
         integer, dimension(4, 4) :: path1x1
         integer, dimension(4, 8) :: path2x2
         integer, dimension(4, 12) :: path3x3
         integer, dimension(8, 6) :: path1x2
         integer, dimension(8, 8) :: path1x3
         integer, dimension(4) :: coordbase
         complex(kind=wc) :: c12, c13, c23

         real(kind=wp), parameter :: k5 = 1.0_WP / 180.0_WP  ! 5 loop
         !real(kind=WP), parameter :: k5 = 1.0_WP / 90.0_WP   ! 3 loop
         real(kind=wp), parameter :: k1 = (19.0_WP / 9.0_WP) - 55.0_WP * k5
         real(kind=wp), parameter :: k2 = (1.0_WP / 36.0_WP) - 16.0_WP * k5
         real(kind=wp), parameter :: k3 = 64.0_WP * k5 - (32.0_WP / 45.0_WP)
         real(kind=wp), parameter :: k4 = (1.0_WP / 15.0_WP) - 6.0_WP * k5
         !            real(kind=WP), parameter :: k5 = 0.0_WP
         !           real(kind=WP), parameter :: k1 = 0.0_WP
         !          real(kind=WP), parameter :: k2 = 0.0_WP
         !         real(kind=WP), parameter :: k3 = 0.0_WP
         !        real(kind=WP), parameter :: k4 = 1.0_WP
         ! counters
         integer :: nnx, nny, nnz, nnt
         integer :: t1, t2

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
                        complex(kind=wc), dimension(3, 3, 4, nt, nx, ny, nz), intent(in) :: data
                        integer, intent(in) :: nt, nx, ny, nz, mu, nu
                        integer, dimension(4) :: coord
                        complex(kind=wc), dimension(3, 3, nt, nx, ny, nz) :: u_xd
                        complex(kind=wc), dimension(3, 3) :: clovleaf
                        integer, dimension(4) :: plaqpath
                        integer, dimension(4) :: coordbase
                        complex(kind=wc) :: c12, c13, c23
                        ! counters
                        integer :: nnx, nny, nnz, nnt
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
                              complex(kind=wc), dimension(3, 3, 4, nt, nx, ny, nz), intent(in) :: data
                              integer, intent(in) :: nt, nx, ny, nz
                              real(kind=wp) :: b
                              !
                              complex(kind=wc), dimension(3, 3, nt, nx, ny, nz, 3) :: bfield
                              !complex(kind=WC), dimension(NT, NX, NY, NZ, 3, 8) :: BFieldAdjoint
                              !complex(kind=WP), dimension(NT,NX,NY,NZ,3,3,3) :: BDag
                              !complex(kind=WP), dimension(NT, NX, NY, NZ, 3) :: BBDagTraceless
                              !real(kind=WP), dimension(NT, NX, NY, NZ, 3, 3, 3) :: BBDagTraceless
                              complex(kind=wc), dimension(3, 3) :: tempeval
                              complex(kind=wc), dimension(8) :: com
                              real(kind=wp) :: cloverplaq, temppl
                              real(kind=wp) :: iiclov
                              ! counter
                              integer :: nnx, nny, nnz, nnt, ii

                              complex(kind=wc) :: ztemp12, ztemp23, ztemp31
                              real(kind=wp) :: trace
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
                                 local(tempeval, com, temppl, ztemp12, ztemp23, ztemp31, trace) &
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
                                                            call realtracemultmatmat(temppl, bfield(nnt, nnx, nny, nnz, :, :, ii), &
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
                                                integer, dimension(7), intent(in) :: datashape
                                                integer, dimension(4), intent(in) :: coord
                                                integer, dimension(4) :: periodcoord
                                                integer :: cc
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
                                                complex(kind=wc), dimension(:, :, :, :, :, :, :), intent(in) :: data
                                                integer, intent(in) :: tau, mu, nu
                                                integer, intent(in) :: nmu, nnu
                                                real(kind=wp) :: loopval
                                                !logical(kind=c_bool),                                intent(in) :: verbose
                                                !"""
                                                !Calculates the nMuxnNu loop at t = tau
                                                ! Nmu link in mu direction
                                                ! nNu links in nu direction
                                                !data is [nt, nx, ny, nz, mu, colour, colour] complex
                                                ! Returns the average value of the loop
                                                !"""
                                                integer, dimension(7) :: datashape
                                                integer, dimension(4) :: mucoord, nucoord, coordbase, coord
                                                integer :: nx, ny, nz, nn  ! Counters
                                                ! For intermediate calculating plaquette
                                                complex(kind=wc), dimension(3, 3) :: umu_x, unu_xpamuh, umu_xpanuh, unu_x
                                                complex(kind=wc), dimension(3, 3) :: utemp, utemp2
                                                complex(kind=wc), dimension(3, 3) :: umuunu, umudagunudag
                                                real(kind=wp) :: p
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
                                                complex(kind=wc), dimension(:, :, :, :, :, :, :), intent(in) :: data
                                                integer, intent(in) :: mu, nu
                                                integer, intent(in) :: nmu, nnu
                                                real(kind=wp) :: loopval
                                                !logical(kind=c_bool),                                intent(in) :: verbose
                                                !"""
                                                !Calculates the nMuxnNu loop
                                                ! Nmu link in mu direction
                                                ! nNu links in nu direction
                                                !data is [nt, nx, ny, nz, mu, colour, colour] complex
                                                ! Returns the average value of the loop
                                                !"""
                                                integer, dimension(7) :: datashape
                                                integer, dimension(4) :: mucoord, nucoord, coordbase, coord
                                                integer :: nx, ny, nz, nt, nn  ! Counters
                                                ! For intermediate calculating plaquette
                                                complex(kind=wc), dimension(3, 3) :: umu_x, unu_xpamuh, umu_xpanuh, unu_x
                                                complex(kind=wc), dimension(3, 3) :: utemp, utemp2
                                                complex(kind=wc), dimension(3, 3) :: umuunu, umudagunudag
                                                real(kind=wp) :: p
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
                                                complex(kind=wc), dimension(:, :, :, :, :, :, :), intent(in) :: data
                                                real(kind=wp), intent(out) :: sumtrp, time
                                                integer, intent(out) :: np
                                                !"""
                                                !Calculates the polyakov loop
                                                !data is [colour, colour, mu, nt, nx, ny, nz] complex
                                                !returns the sum of polyakov, number of loops measured,
                                                !the average polyakov and the time taken to calculate it
                                                !"""
                                                integer, dimension(7) :: datashape
                                                integer, dimension(4) :: coord
                                                integer :: nx, ny, nz, nt  ! Counters
                                                complex(kind=wc), dimension(3, 3) :: unu_x, unu_xt, unu_temp
                                                real(kind=wp) :: p
                                                ! Timers
                                                real(kind=wp) :: start, end
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
