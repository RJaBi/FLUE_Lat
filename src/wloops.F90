!! Calculate, spatial, temporal plaquettes
!! Compiles for use with python using
!! 'f2py -c fortPlaq.f90 -m fortPlaq'
!!
!! Which can then be imported into python as
!! import fortPlaq
!! fPlaq = fortPlaq.plaq.plaquette
!! which would give you the plaquette function

module FLUE_wloops
   use FLUE_constants, only: WP, WC, SP
   use FLUE_SU3MatrixOps, only: MultiplyMatMat, MultiplyMatdagMatdag, &
                                RealTraceMultMatMat, Ident, tracelessconjgsubtract, colourDecomp, RealTraceMat
   use M_stopwatch, only: watchtype, create_watch, start_watch, stop_watch, destroy_watch, read_watch
   use stdlib_io_npy, only: save_npy
   implicit none(external)
   public

contains

   pure function genericPath(data, coordBase, path) result(U_xd)
      complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(in) :: data
      integer, dimension(4), intent(in) :: coordBase
      integer, dimension(:), intent(in) :: path
      complex(kind=WC), dimension(3, 3) :: U_xd
      ! internal
      complex(kind=WC), dimension(3, 3) :: amat, bmat
      integer, dimension(4) :: coord, pCoord
      ! counters
      integer :: pp
      integer :: mu
      coord = coordBase
      U_xd = Ident
      do pp = 1, SIZE(path)
         bmat = U_xd
         mu = path(pp)
         pCoord = 0
         pCoord(ABS(mu)) = 1
         if (mu < 0) then
            ! the step is backwards so subtract off from coord first
            coord = coord - pCoord
            coord = periodCoord(coord, SHAPE(data))
            ! get the link here going forward in mu
            ! and dagger it
            amat = CONJG(TRANSPOSE(data(coord(1), coord(2), coord(3), coord(4), ABS(mu), :, :)))
         else
            amat = data(coord(1), coord(2), coord(3), coord(4), mu, :, :)
            coord = coord + pCoord
            coord = periodCoord(coord, SHAPE(data))
         end if
         ! now multiply it into U_xd from the right
         call MultiplyMatMat(U_xd, bmat, amat)
      end do
      U_xd = TRANSPOSE(U_xd)
   end function genericPath

   subroutine genPlaquette(data, NT, NX, NY, NZ, muStart, muEnd, nuEnd, sumTrP, nP, time)
      complex(kind=WC), dimension(NT, NX, NY, NZ, 4, 3, 3), intent(in) :: data
      integer, intent(in) :: muStart, muEnd, nuEnd
      integer, intent(in) :: NT, NX, NY, NZ
      real(kind=WP), intent(out) :: sumTrP, time
      integer, intent(out) :: nP
      ! Counters
      integer :: mu, nu, nnx, nny, nnz, nnt
      ! other variables
      complex(kind=WC), dimension(3, 3) :: plaq
      integer, dimension(7) :: dataShape
      integer, dimension(4) :: plaqPath, coordBase
      real(kind=WP) :: P
      ! Timers
      type(watchtype) :: watch
      real(Kind=SP) :: watchtime
      call create_watch(watch)
      call start_watch(watch)
      dataShape = (/NT, NX, NY, NZ, 4, 3, 3/)
      !# hold the sum
      sumTrP = 0.0_WP
      !# hold the number measured
      nP = 0
      do mu = muStart, muEnd
         do nu = mu + 1, nuEnd
            plaqPath = (/mu, nu, -mu, -nu/)
#ifdef LOCALITYSUPPORT
            do concurrent(nnx=1:nx, nny=1:ny, nnz=1:nz, nnt=1:nt) &
               local(coordBase, plaq, P) local_init(nnx, nny, nnz, nnt) reduce(+:sumTrP, nP)
#else
               do concurrent(nnx=1:nx, nny=1:ny, nnz=1:nz, nnt=1:nt)
#endif
                  coordBase = (/nnt, nnx, nny, nnz/)
                  !plaq = genericPath(data, coordBase, plaqPath)

                  call RealTraceMultMatMat(P, Ident, plaq)
                  sumTrP = sumTrP + P
                  nP = nP + 1
               end do
            end do
         end do
         call stop_watch(watch)
         call read_watch(watchTime, watch, 'wall')
         time = 0.0_WP
         time = watchTime
         !time = end - start
         end subroutine genPlaquette

         pure function cloverLoopMuNuCoord(data, NT, NX, NY, NZ, mu, nu) result(U_xd)
            complex(kind=WC), dimension(NT, NX, NY, NZ, 4, 3, 3), intent(in) :: data
            integer, intent(in) :: NT, NX, NY, NZ, mu, nu
            integer, dimension(4) :: coord
            complex(kind=WC), dimension(NT, NX, NY, NZ, 3, 3) :: U_xd
            complex(kind=WC), dimension(3, 3) :: clovLeaf, thisU
            integer, dimension(4, 4) :: plaqPath
            integer, dimension(4) :: coordBase
            complex(kind=WC) :: c12, c13, c23
            ! counters
            integer :: nnx, nny, nnz, nnt

            integer :: t1, t2
            ! top left
            plaqPath(1, :) = (/nu, -mu, -nu, mu/)
            ! top right
            plaqPath(2, :) = (/mu, nu, -mu, -nu/)
            ! bottom left
            plaqPath(3, :) = (/-mu, -nu, mu, nu/)
            ! bottom right
            plaqPath(4, :) = (/-nu, mu, nu, -mu/)
#ifdef LOCALITYSUPPORT
            do concurrent(nnx=1:nx, nny=1:ny, nnz=1:nz, nnt=1:nt) &
               default(none) local_init(nnx, nny, nnz, nnt, plaqPath) &
               local(coordBase, clovLeaf, trcnsub, thisU, c12, c13, c23) shared(U_xd)
#else
               do concurrent(nnx=1:nx, nny=1:ny, nnz=1:nz, nnt=1:nt)
#endif
                  ! Calculate the clover for this coordinate
                  coordBase = (/nnt, nnx, nny, nnz/)
                  ! 1
                  !clovLeaf = genericPath(data, coordBase, plaqPath(1, :))
                  thisU = +clovLeaf
                  ! 2
                  !clovLeaf = genericPath(data, coordBase, plaqPath(2, :))
                  thisU = thisU + clovLeaf
                  ! 3
                  !clovLeaf = genericPath(data, coordBase, plaqPath(3, :))
                  thisU = thisU + clovLeaf
                  ! 4
                  !clovLeaf = genericPath(data, coordBase, plaqPath(4, :))
                  thisU = thisU + clovLeaf
                  ! Projection & factors
                  ! projection to anti-hermition but not traceless matrix
                  ! as in Borsanyi wilson_flow.c
                  ! a_chm = 0.5*( a - conj(b))
                  c12 = 0.5_WP * (thisU(1, 2) - CONJG(thisU(2, 1)))
                  c13 = 0.5_WP * (thisU(1, 3) - CONJG(thisU(3, 1)))
                  c23 = 0.5_WP * (thisU(2, 3) - CONJG(thisU(3, 2)))
                  thisU(1, 2) = c12
                  thisU(1, 3) = c13
                  thisU(2, 3) = c23
                  thisU(2, 1) = -CONJG(c12)
                  thisU(3, 1) = -CONJG(c13)
                  thisU(3, 2) = -CONJG(c23)
                  thisU(1, 1) = CMPLX(0.0_WP, AIMAG(thisU(1, 1)))
                  thisU(2, 2) = CMPLX(0.0_WP, AIMAG(thisU(2, 2)))
                  thisU(3, 3) = CMPLX(0.0_WP, AIMAG(thisU(3, 3)))
                  U_xd(nnt, nnx, nny, nnz, :, :) = thisU * 0.25_WP
               end do
               end function cloverLoopMuNuCoord

               pure function Loop5MuNuCoord(data, NT, NX, NY, NZ, mu, nu) result(U_xd)
                  ! See hep-lat/0203008
                  complex(kind=WC), dimension(NT, NX, NY, NZ, 4, 3, 3), intent(in) :: data
                  integer, intent(in) :: NT, NX, NY, NZ, mu, nu
                  integer, dimension(4) :: coord
                  complex(kind=WC), dimension(NT, NX, NY, NZ, 3, 3) :: U_xd
                  complex(kind=WC), dimension(3, 3) :: clovLeaf, thisU, workU
                  integer, dimension(4, 4) :: path1x1
                  integer, dimension(4, 8) :: path2x2
                  integer, dimension(4, 12) :: path3x3
                  integer, dimension(8, 6) :: path1x2
                  integer, dimension(8, 8) :: path1x3
                  integer, dimension(4) :: coordBase
                  complex(kind=WC) :: c12, c13, c23

                  real(kind=WP), parameter :: k5 = 1.0_WP / 180.0_WP  ! 5 loop
                  !real(kind=WP), parameter :: k5 = 1.0_WP / 90.0_WP   ! 3 loop
                  real(kind=WP), parameter :: k1 = (19.0_WP / 9.0_WP) - 55.0_WP * k5
                  real(kind=WP), parameter :: k2 = (1.0_WP / 36.0_WP) - 16.0_WP * k5
                  real(kind=WP), parameter :: k3 = 64.0_WP * k5 - (32.0_WP / 45.0_WP)
                  real(kind=WP), parameter :: k4 = (1.0_WP / 15.0_WP) - 6.0_WP * k5
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

#ifdef LOCALITYSUPPORT
                  do concurrent(nnx=1:nx, nny=1:ny, nnz=1:nz, nnt=1:nt) &
                     default(none) &
                     local_init(nnx, nny, nnz, nnt, path1x1, path1x2, path2x2, path1x3, path3x3) &
                     local(coordBase, clovLeaf, trcnsub, thisU, c12, c13, c23, workU, t1, t2) shared(U_xd)
#else
                     do nnx = 1, nx; do nny = 1, ny; do nnz = 1, nz; do nnt = 1, nt;
#endif

                                 ! Calculate the clover 1x1 for this coordinate
                                 coordBase = (/nnt, nnx, nny, nnz/)
                                 ! 1
                                 clovLeaf = genericPath(data, coordBase, path1x1(1, :))
                                 thisU = +clovLeaf
                                 ! 2
                                 clovLeaf = genericPath(data, coordBase, path1x1(2, :))
                                 thisU = thisU + clovLeaf
                                 ! 3
                                 clovLeaf = genericPath(data, coordBase, path1x1(3, :))
                                 thisU = thisU + clovLeaf
                                 ! 4
                                 clovLeaf = genericPath(data, coordBase, path1x1(4, :))
                                 thisU = thisU + clovLeaf
                                 thisU = thisU * k1
                                 ! 2x2
                                 ! 1
                                 clovLeaf = genericPath(data, coordBase, path2x2(1, :))
                                 workU = +clovLeaf
                                 ! 2
                                 clovLeaf = genericPath(data, coordBase, path2x2(2, :))
                                 workU = workU + clovLeaf
                                 ! 3
                                 clovLeaf = genericPath(data, coordBase, path2x2(3, :))
                                 workU = workU + clovLeaf
                                 ! 4
                                 clovLeaf = genericPath(data, coordBase, path2x2(4, :))
                                 workU = workU + clovLeaf
                                 thisU = thisU + workU * k2

                        !!!!!
                        !!!!! 3x3
                        !!!!!
                                 ! 1
                                 clovLeaf = genericPath(data, coordBase, path3x3(1, :))
                                 workU = +clovLeaf
                                 ! 2
                                 clovLeaf = genericPath(data, coordBase, path3x3(2, :))
                                 workU = workU + clovLeaf
                                 ! 3
                                 clovLeaf = genericPath(data, coordBase, path3x3(3, :))
                                 workU = workU + clovLeaf
                                 ! 4
                                 clovLeaf = genericPath(data, coordBase, path3x3(4, :))
                                 workU = workU + clovLeaf
                                 thisU = thisU + workU * k5

                                 ! 1x2
                                 ! 1
                                 clovLeaf = genericPath(data, coordBase, path1x2(1, :))
                                 workU = +clovLeaf
                                 clovLeaf = genericPath(data, coordBase, path1x2(5, :))
                                 workU = workU + clovLeaf
                                 ! 2
                                 clovLeaf = genericPath(data, coordBase, path1x2(2, :))
                                 workU = workU + clovLeaf
                                 clovLeaf = genericPath(data, coordBase, path1x2(6, :))
                                 workU = workU + clovLeaf
                                 ! 3
                                 clovLeaf = genericPath(data, coordBase, path1x2(3, :))
                                 workU = workU + clovLeaf
                                 clovLeaf = genericPath(data, coordBase, path1x2(7, :))
                                 workU = workU + clovLeaf
                                 ! 4
                                 clovLeaf = genericPath(data, coordBase, path1x2(4, :))
                                 workU = workU + clovLeaf
                                 clovLeaf = genericPath(data, coordBase, path1x2(8, :))

                                 workU = workU + clovLeaf
                                 thisU = thisU + workU * k3 * 0.5_WP

                                 ! 1x3
                                 ! 1
                                 clovLeaf = genericPath(data, coordBase, path1x3(1, :))
                                 workU = +clovLeaf
                                 clovLeaf = genericPath(data, coordBase, path1x3(5, :))
                                 workU = workU + clovLeaf
                                 ! 2
                                 clovLeaf = genericPath(data, coordBase, path1x3(2, :))
                                 workU = workU + clovLeaf
                                 clovLeaf = genericPath(data, coordBase, path1x3(6, :))
                                 workU = workU + clovLeaf
                                 ! 3
                                 clovLeaf = genericPath(data, coordBase, path1x3(3, :))
                                 workU = workU + clovLeaf
                                 clovLeaf = genericPath(data, coordBase, path1x3(7, :))
                                 workU = workU + clovLeaf
                                 ! 4
                                 clovLeaf = genericPath(data, coordBase, path1x3(4, :))
                                 workU = workU + clovLeaf
                                 clovLeaf = genericPath(data, coordBase, path1x3(8, :))
                                 workU = workU + clovLeaf
                                 thisU = thisU + workU * k4 * 0.5_WP

                                 ! Projection & factors
                                 ! projection to anti-hermition but not traceless matrix
                                 ! as in Borsanyi wilson_flow.c
                                 ! a_chm = 0.5*( a - conj(b))
                                 c12 = 0.5_WP * (thisU(1, 2) - CONJG(thisU(2, 1)))
                                 c13 = 0.5_WP * (thisU(1, 3) - CONJG(thisU(3, 1)))
                                 c23 = 0.5_WP * (thisU(2, 3) - CONJG(thisU(3, 2)))
                                 thisU(1, 2) = c12
                                 thisU(1, 3) = c13
                                 thisU(2, 3) = c23
                                 thisU(2, 1) = -CONJG(c12)
                                 thisU(3, 1) = -CONJG(c13)
                                 thisU(3, 2) = -CONJG(c23)
                                 thisU(1, 1) = CMPLX(0.0_WP, AIMAG(thisU(1, 1)))
                                 thisU(2, 2) = CMPLX(0.0_WP, AIMAG(thisU(2, 2)))
                                 thisU(3, 3) = CMPLX(0.0_WP, AIMAG(thisU(3, 3)))
                                 U_xd(nnt, nnx, nny, nnz, :, :) = thisU * 0.25_WP
#ifdef LOCALITYSUPPORT
                              end do
#else
                           end do; end do; end do; end do
#endif
               end function Loop5MuNuCoord

               function plaquetteMuNuCoord(data, NT, NX, NY, NZ, mu, nu) result(U_xd)
                  complex(kind=WC), dimension(NT, NX, NY, NZ, 4, 3, 3), intent(in) :: data
                  integer, intent(in) :: NT, NX, NY, NZ, mu, nu
                  integer, dimension(4) :: coord
                  complex(kind=WC), dimension(NT, NX, NY, NZ, 3, 3) :: U_xd
                  complex(kind=WC), dimension(3, 3) :: clovLeaf
                  integer, dimension(4) :: plaqPath
                  integer, dimension(4) :: coordBase
                  complex(kind=WC) :: c12, c13, c23
                  ! counters
                  integer :: nnx, nny, nnz, nnt
                  ! top left
                  ! Untested
                  stop
                  plaqPath = (/nu, -mu, -nu, mu/)
#ifdef LOCALITYSUPPORT
                  do concurrent(nnx=1:nx, nny=1:ny, nnz=1:nz, nnt=1:nt) &
                     default(none) local_init(nnx, nny, nnz, nnt, plaqPath) &
                     local(coordBase, clovLeaf, c12, c13, c23) shared(U_xd)
#else
                     do concurrent(nnx=1:nx, nny=1:ny, nnz=1:nz, nnt=1:nt)
#endif
                        ! Calculate the clover for this coordinate
                        coordBase = (/nnt, nnx, nny, nnz/)
                        !clovLeaf = genericPath(data, coordBase, plaqPath)
                        ! Projection & factors
                        ! projection to anti-hermition but not traceless matrix
                        ! as in Borsanyi wilson_flow.c
                        ! a_chm = 0.5*( a - conj(b))
                        c12 = 0.5_WP * (clovLeaf(1, 2) - CONJG(clovLeaf(2, 1)))
                        c13 = 0.5_WP * (clovLeaf(1, 3) - CONJG(clovLeaf(3, 1)))
                        c23 = 0.5_WP * (clovLeaf(2, 3) - CONJG(clovLeaf(3, 2)))
                        clovLeaf(1, 2) = c12
                        clovLeaf(1, 3) = c13
                        clovLeaf(2, 3) = c23
                        clovLeaf(2, 1) = -CONJG(c12)
                        clovLeaf(3, 1) = -CONJG(c13)
                        clovLeaf(3, 2) = -CONJG(c23)
                        clovLeaf(1, 1) = CMPLX(0.0_WP, AIMAG(clovLeaf(1, 1)))
                        clovLeaf(2, 2) = CMPLX(0.0_WP, AIMAG(clovLeaf(2, 2)))
                        clovLeaf(3, 3) = CMPLX(0.0_WP, AIMAG(clovLeaf(3, 3)))
                        U_xd(nnt, nnx, nny, nnz, :, :) = clovLeaf * 0.25_WP
                     end do
                     end function plaquetteMuNuCoord

                     function magnetic(data, NT, NX, NY, NZ) result(B)
                        complex(kind=WC), dimension(NT, NX, NY, NZ, 4, 3, 3), intent(in) :: data
                        integer, intent(in) :: NT, NX, NY, NZ
                        real(kind=WP) :: B
                        !
                        complex(kind=WC), dimension(NT, NX, NY, NZ, 3, 3, 3) :: BField
                        !complex(kind=WC), dimension(NT, NX, NY, NZ, 3, 8) :: BFieldAdjoint
                        !complex(kind=WP), dimension(NT,NX,NY,NZ,3,3,3) :: BDag
                        !complex(kind=WP), dimension(NT, NX, NY, NZ, 3) :: BBDagTraceless
                        !real(kind=WP), dimension(NT, NX, NY, NZ, 3, 3, 3) :: BBDagTraceless
                        complex(kind=WC), dimension(3, 3) :: tempEval
                        complex(kind=WC), dimension(8) :: com
                        real(kind=WP) :: cloverPlaq, tempPl
                        real(kind=WP) :: iiclov
                        ! counter
                        integer :: nnx, nny, nnz, nnt, ii

                        complex(kind=WC) :: ztemp12, ztemp23, ztemp31, temp11, temp22
                        real(kind=WP) :: trace
                        ! The indices are so that they match the wilson_flow.c values exactly
                        ! there t,x,y,z = 0,3,2,1
                        ! here t,x,y,z = 1,2,3,4
                        ! match 23
                        BField(:, :, :, :, :, :, 1) = Loop5MuNuCoord(data, NT, NX, NY, NZ, 3, 2)
                        ! match 31
                        BField(:, :, :, :, :, :, 2) = Loop5MuNuCoord(data, NT, NX, NY, NZ, 2, 4)
                        ! match 12
                        BField(:, :, :, :, :, :, 3) = Loop5MuNuCoord(data, NT, NX, NY, NZ, 4, 3)

                        cloverPlaq = 0.0_WP
#ifdef LOCALITYSUPPORT
                        do concurrent(nny=1:ny, nnz=1:nz, nnt=1:nt, nnx=1:nx, ii=1:3) &
                           default(none) local_init(nnx, nny, nnz, nnt, ii) &
                           local(tempEval, com, tempPl, ztemp12, ztemp23, ztemp31, temp11, temp22, temp33, trace) &
                           reduce(+:cloverPlaq)
#else
                           do concurrent(nny=1:ny, nnz=1:nz, nnt=1:nt, nnx=1:nx, ii=1:3)
#endif
                              call RealTraceMultMatMat(tempPl, BField(nnt, nnx, nny, nnz, :, :, ii), &
                                                       BField(nnt, nnx, nny, nnz, :, :, ii))
                              cloverPlaq = cloverPlaq - tempPl
                           end do
                           B = cloverPlaq / real(NT * NX * NY * NZ, kind=WP)
                           end function magnetic

                           pure function periodCoord(coord, dataShape)
                              integer, dimension(7), intent(in) :: dataShape
                              integer, dimension(4), intent(in) :: coord
                              integer, dimension(4) :: periodCoord
                              ! A lazy function to handle the periodic boundary conditions
                              ! dataShape is (NT, nx, ny, nz, 4, 3, 3)
                              ! coord is (nt, nx, ny, nz)
                              ! checks if the value in coord is greater than corresponding N in datashape
                              ! if so sets it to 1
                              ! i.e. only handles steps of 1
                              integer :: cc
                              periodCoord = coord
                              do cc = 1, SIZE(coord)
                                 if (coord(cc) > dataShape(cc)) then
                                    periodCoord(cc) = 1
                                 else if (coord(cc) == 0) then
                                    periodCoord(cc) = dataShape(cc)
                                 end if
                              end do
                           end function periodCoord

                           pure function wmunutau(data, tau, mu, nu, nMu, nNu) result(loopVal)
                              complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(in) :: data
                              integer, intent(in) :: tau, mu, nu
                              integer, intent(in) :: nMu, nNu
                              real(kind=WP) :: loopVal
                              !logical(kind=c_bool),                                intent(in) :: verbose
                              !"""
                              !Calculates the nMuxnNu loop at t = tau
                              ! Nmu link in mu direction
                              ! nNu links in nu direction
                              !data is [nt, nx, ny, nz, mu, colour, colour] complex
                              ! Returns the average value of the loop
                              !"""
                              integer, dimension(7) :: dataShape
                              integer, dimension(4) :: muCoord, nuCoord, coordBase, coord
                              integer :: nx, ny, nz, nn  ! Counters
                              ! For intermediate calculating plaquette
                              complex(kind=WC), dimension(3, 3) :: Umu_x, Unu_xpamuh, Umu_xpanuh, Unu_x
                              complex(kind=WC), dimension(3, 3) :: Utemp, Utemp2
                              complex(kind=WC), dimension(3, 3) :: UmuUnu, UmudagUnudag
                              real(kind=WP) :: P
                              dataShape = SHAPE(data)
                              !# hold the sum
                              loopVal = 0.0_WP

                              muCoord(:) = 0
                              nuCoord(:) = 0
                              !# This is a single shift in mu
                              muCoord(mu) = 1
                              !# This is a single shift in nu
                              nuCoord(nu) = 1
                              !# loop over all sites
                              do nx = 1, dataShape(2)
                                 do ny = 1, dataShape(3)
                                    do nz = 1, dataShape(4)
                                       ! The starting coordinate
                                       coordBase = (/tau, nx, ny, nz/)
                                       ! Do U_nMu(x)
                                       coord = coordBase
                                       Umu_x = Ident
                                       Utemp2 = Ident
                                       do nn = 1, nMu
                                          ! Get data
                                          Utemp = data(coord(1), coord(2), coord(3), coord(4), mu, :, :)
                                          ! Multiply it in
                                          call MultiplyMatMat(Utemp2, Umu_x, Utemp)
                                          ! and re-assign
                                          Umu_x = Utemp2
                                          ! Update the coordinate
                                          coord = periodCoord(coord + muCoord, dataShape)
                                       end do
                                       ! Do U_nNu(x+nMu*amu)
                                       coord = periodCoord(coordBase + muCoord, dataShape)
                                       do nn = 2, nMu
                                          coord = periodCoord(coord + muCoord, dataShape)
                                       end do
                                       Unu_xpamuh = Ident
                                       Utemp2 = Ident
                                       do nn = 1, nNu
                                          ! Get data
                                          Utemp = data(coord(1), coord(2), coord(3), coord(4), nu, :, :)
                                          ! Multiply it in
                                          call MultiplyMatMat(Utemp2, Unu_xpamuh, Utemp)
                                          ! and re-assign
                                          Unu_xpamuh = Utemp2
                                          ! Update the coordinate
                                          coord = periodCoord(coord + nuCoord, dataShape)
                                       end do
                                       ! U_nMu(x+nNu)
                                       coord = periodCoord(coordBase + nuCoord, dataShape)
                                       do nn = 2, nNu
                                          coord = periodCoord(coord + nuCoord, dataShape)
                                       end do
                                       Umu_xpanuh = Ident
                                       Utemp2 = Ident
                                       do nn = 1, nMu
                                          ! Get data
                                          Utemp = data(coord(1), coord(2), coord(3), coord(4), mu, :, :)
                                          ! Multiply it in
                                          call MultiplyMatMat(Utemp2, Umu_xpanuh, Utemp)
                                          ! and re-assign
                                          Umu_xpanuh = Utemp2
                                          ! Update the coordinate
                                          coord = periodCoord(coord + muCoord, dataShape)
                                       end do
                                       ! U_nNu(x)
                                       coord = coordBase
                                       Unu_x = Ident
                                       Utemp2 = Ident
                                       do nn = 1, nNu
                                          ! Get data
                                          Utemp = data(coord(1), coord(2), coord(3), coord(4), nu, :, :)
                                          ! Multiply it in
                                          call MultiplyMatMat(Utemp2, Unu_x, Utemp)
                                          ! and re-assign
                                          Unu_x = Utemp2
                                          ! Update the coordinate
                                          coord = periodCoord(coord + nuCoord, dataShape)
                                       end do
                                       !# Multiply bottom, right together
                                       call MultiplyMatMat(UmuUnu, Umu_x, Unu_xpamuh)
                                       !# Multiply left, top together, take dagger
                                       call MultiplyMatdagMatdag(UmudagUnudag, Umu_xpanuh, Unu_x)
                                       !# multiply two halves together, take trace
                                       call RealTraceMultMatMat(P, UmuUnu, UmudagUnudag)
                                       loopVal = loopVal + P
                                    end do
                                 end do
                              end do
                              loopVal = loopVal / real(dataShape(2) * dataShape(3) * dataShape(4), kind=WP)
                           end function wmunutau

                           pure function wmunu(data, mu, nu, nMu, nNu) result(loopVal)
                              complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(in) :: data
                              integer, intent(in) :: mu, nu
                              integer, intent(in) :: nMu, nNu
                              real(kind=WP) :: loopVal
                              !logical(kind=c_bool),                                intent(in) :: verbose
                              !"""
                              !Calculates the nMuxnNu loop
                              ! Nmu link in mu direction
                              ! nNu links in nu direction
                              !data is [nt, nx, ny, nz, mu, colour, colour] complex
                              ! Returns the average value of the loop
                              !"""
                              integer, dimension(7) :: dataShape
                              integer, dimension(4) :: muCoord, nuCoord, coordBase, coord
                              integer :: nx, ny, nz, nt, nn  ! Counters
                              ! For intermediate calculating plaquette
                              complex(kind=WC), dimension(3, 3) :: Umu_x, Unu_xpamuh, Umu_xpanuh, Unu_x
                              complex(kind=WC), dimension(3, 3) :: Utemp, Utemp2
                              complex(kind=WC), dimension(3, 3) :: UmuUnu, UmudagUnudag
                              real(kind=WP) :: P
                              dataShape = SHAPE(data)
                              !# hold the sum
                              loopVal = 0.0_WP

                              muCoord(:) = 0
                              nuCoord(:) = 0
                              !# This is a single shift in mu
                              muCoord(mu) = 1
                              !# This is a single shift in nu
                              nuCoord(nu) = 1
                              !# loop over all sites
                              do nx = 1, dataShape(2)
                                 do ny = 1, dataShape(3)
                                    do nz = 1, dataShape(4)
                                       do nt = 1, dataShape(1)
                                          ! The starting coordinate
                                          coordBase = (/nt, nx, ny, nz/)
                                          ! Do U_nMu(x)
                                          coord = coordBase
                                          Umu_x = Ident
                                          Utemp2 = Ident
                                          do nn = 1, nMu
                                             ! Get data
                                             Utemp = data(coord(1), coord(2), coord(3), coord(4), mu, :, :)
                                             ! Multiply it in
                                             call MultiplyMatMat(Utemp2, Umu_x, Utemp)
                                             ! and re-assign
                                             Umu_x = Utemp2
                                             ! Update the coordinate
                                             coord = periodCoord(coord + muCoord, dataShape)
                                          end do
                                          ! Do U_nNu(x+nMu*amu)
                                          coord = periodCoord(coordBase + muCoord, dataShape)
                                          do nn = 2, nMu
                                             coord = periodCoord(coord + muCoord, dataShape)
                                          end do
                                          Unu_xpamuh = Ident
                                          Utemp2 = Ident
                                          do nn = 1, nNu
                                             ! Get data
                                             Utemp = data(coord(1), coord(2), coord(3), coord(4), nu, :, :)
                                             ! Multiply it in
                                             call MultiplyMatMat(Utemp2, Unu_xpamuh, Utemp)
                                             ! and re-assign
                                             Unu_xpamuh = Utemp2
                                             ! Update the coordinate
                                             coord = periodCoord(coord + nuCoord, dataShape)
                                          end do
                                          ! U_nMu(x+nNu)
                                          coord = periodCoord(coordBase + nuCoord, dataShape)
                                          do nn = 2, nNu
                                             coord = periodCoord(coord + nuCoord, dataShape)
                                          end do
                                          Umu_xpanuh = Ident
                                          Utemp2 = Ident
                                          do nn = 1, nMu
                                             ! Get data
                                             Utemp = data(coord(1), coord(2), coord(3), coord(4), mu, :, :)
                                             ! Multiply it in
                                             call MultiplyMatMat(Utemp2, Umu_xpanuh, Utemp)
                                             ! and re-assign
                                             Umu_xpanuh = Utemp2
                                             ! Update the coordinate
                                             coord = periodCoord(coord + muCoord, dataShape)
                                          end do
                                          ! U_nNu(x)
                                          coord = coordBase
                                          Unu_x = Ident
                                          Utemp2 = Ident
                                          do nn = 1, nNu
                                             ! Get data
                                             Utemp = data(coord(1), coord(2), coord(3), coord(4), nu, :, :)
                                             ! Multiply it in
                                             call MultiplyMatMat(Utemp2, Unu_x, Utemp)
                                             ! and re-assign
                                             Unu_x = Utemp2
                                             ! Update the coordinate
                                             coord = periodCoord(coord + nuCoord, dataShape)
                                          end do
                                          !# Multiply bottom, right together
                                          call MultiplyMatMat(UmuUnu, Umu_x, Unu_xpamuh)
                                          !# Multiply left, top together, take dagger
                                          call MultiplyMatdagMatdag(UmudagUnudag, Umu_xpanuh, Unu_x)
                                          !# multiply two halves together, take trace
                                          call RealTraceMultMatMat(P, UmuUnu, UmudagUnudag)
                                          loopVal = loopVal + P
                                       end do
                                    end do
                                 end do
                              end do
                              loopVal = loopVal / real(SIZE(data) / 36, kind=WP)
                           end function wmunu

!DOESNTWORK!   function WxNMuNu(data, NT, NX, NY, NZ, mu, nu, WSize, nSize) result(loopAverage)
!DOESNTWORK!     complex(kind=WC), dimension(NT, NX, NY, NZ, 4, 3, 3), intent(in) :: data
!DOESNTWORK!     integer, intent(in) :: NT, NX, NY, NZ
!DOESNTWORK!     integer, intent(in) :: mu, nu, wSize, nSize
!DOESNTWORK!     real(kind=WP), dimension(NT) :: loopAverage
!DOESNTWORK!     integer, dimension(4) :: coord
!DOESNTWORK!     ! counters
!DOESNTWORK!     integer :: nnx, nny, nnz, nnt
!DOESNTWORK!     loopAverage = 0.0_WP
!DOESNTWORK!     do nnt=1, NT
!DOESNTWORK!        do concurrent(nnx=1: NX, nny=1: NY, nnz=1: NZ)
!DOESNTWORK!           coord = (/nnt, nnx, nny, nnz/)
!DOESNTWORK!           loopAverage(nnt) = loopAverage(nnt) + WxNMuNuCoord(data, coord, mu, nu, WSize, NSize)
!DOESNTWORK!        end do
!DOESNTWORK!        loopAverage(nnt) = loopAverage(nnt) / real(NX*NY*NZ, kind=WP)
!DOESNTWORK!     end do
!DOESNTWORK!
!DOESNTWORK!   end function WxNMuNu
!DOESNTWORK!
!DOESNTWORK!
!DOESNTWORK!   pure function WxNMuNuCoord(data, coordBase, mu, nu, wSize, nSize) result(loopVal)
!DOESNTWORK!     complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(in) :: data
!DOESNTWORK!     integer, dimension(4), intent(in) :: coordBase
!DOESNTWORK!     integer, intent(in) :: mu, nu, wSize, nSize
!DOESNTWORK!     real(kind=WP) :: loopVal
!DOESNTWORK!     !
!DOESNTWORK!     integer, dimension(7) :: dataShape
!DOESNTWORK!     complex(kind=WC), dimension(3, 3) :: amat, bmat, cmat, dmat
!DOESNTWORK!     integer, dimension(4) :: muCoord, nuCoord, coord
!DOESNTWORK!     ! counters
!DOESNTWORK!     integer :: nn, ww, cc
!DOESNTWORK!     !"""
!DOESNTWORK!     ! Calculates the (real) loop in the mu-nu plane of size w in mu and n in nu
!DOESNTWORK!     ! at location coord = [t,x,y,z]
!DOESNTWORK!     ! data is [nt, nx, ny, nz, mu, colour, colour] complex
!DOESNTWORK!      dataShape = SHAPE(data)
!DOESNTWORK!
!DOESNTWORK!      ! Setting up the mu, nu incrementers
!DOESNTWORK!      muCoord(:) = 0
!DOESNTWORK!      nuCoord(:) = 0
!DOESNTWORK!      muCoord(mu) = 1
!DOESNTWORK!      nuCoord(nu) = 1
!DOESNTWORK!      ! Starting coordinate
!DOESNTWORK!      coord = coordBase
!DOESNTWORK!      ! Store the ongoing result in cmat
!DOESNTWORK!      cmat = data(coord(1), coord(2), coord(3), coord(4), mu, :, :)
!DOESNTWORK!      ! Do the links forward in mu
!DOESNTWORK!      do ww=1, wSize
!DOESNTWORK!         coord = coord + muCoord
!DOESNTWORK!         ! Respectiving periodic BC
!DOESNTWORK!         do cc = 1, SIZE(coord)
!DOESNTWORK!            if (coord(cc) > dataShape(cc)) then
!DOESNTWORK!               coord(cc) = 1
!DOESNTWORK!            end if
!DOESNTWORK!         end do
!DOESNTWORK!         ! Get the link for this coord
!DOESNTWORK!         amat = data(coord(1), coord(2), coord(3), coord(4), mu, :, :)
!DOESNTWORK!         ! Multiply onto existing links
!DOESNTWORK!         ! output, left, right
!DOESNTWORK!         call MultiplyMatMat(bmat, cmat, amat)
!DOESNTWORK!         cmat = bmat
!DOESNTWORK!      end do
!DOESNTWORK!      ! Do the links forward in nu
!DOESNTWORK!      cmat = data(coord(1), coord(2), coord(3), coord(4), nu, :, :)
!DOESNTWORK!      do nn=1, nSize
!DOESNTWORK!         coord = coord + nuCoord
!DOESNTWORK!         ! Respectiving periodic BC
!DOESNTWORK!         do cc = 1, SIZE(coord)
!DOESNTWORK!            if (coord(cc) > dataShape(cc)) then
!DOESNTWORK!               coord(cc) = 1
!DOESNTWORK!            end if
!DOESNTWORK!         end do
!DOESNTWORK!         amat = data(coord(1), coord(2), coord(3), coord(4), nu, :, :)
!DOESNTWORK!         call MultiplyMatMat(bmat, cmat, amat)
!DOESNTWORK!         cmat = bmat
!DOESNTWORK!      end do
!DOESNTWORK!      ! Do the links backwards in mu
!DOESNTWORK!      do nn=1, nSize
!DOESNTWORK!         coord = coordBase + nuCoord
!DOESNTWORK!         do cc = 1, SIZE(coord)
!DOESNTWORK!            if (coord(cc) > dataShape(cc)) then
!DOESNTWORK!               coord(cc) = 1
!DOESNTWORK!            end if
!DOESNTWORK!         end do
!DOESNTWORK!      end do
!DOESNTWORK!      dmat = data(coord(1), coord(2), coord(3), coord(4), mu, :, :)
!DOESNTWORK!      do ww=1, wSize
!DOESNTWORK!         coord = coord +  muCoord
!DOESNTWORK!         ! Respectiving periodic BC
!DOESNTWORK!         do cc = 1, SIZE(coord)
!DOESNTWORK!            if (coord(cc) > dataShape(cc)) then
!DOESNTWORK!               coord(cc) = 1
!DOESNTWORK!            end if
!DOESNTWORK!         end do
!DOESNTWORK!         amat = data(coord(1), coord(2), coord(3), coord(4), mu, :, :)
!DOESNTWORK!         call MultiplyMatMat(bmat, dmat, amat)
!DOESNTWORK!         dmat = bmat
!DOESNTWORK!      end do
!DOESNTWORK!      call MultiplyMatMat(bmat, cmat, conjg(transpose(dmat)))
!DOESNTWORK!      cmat = bmat
!DOESNTWORK!      ! Do the links backwards in nu
!DOESNTWORK!      coord = coordBase
!DOESNTWORK!      dmat = data(coord(1), coord(2), coord(3), coord(4), nu, :, :)
!DOESNTWORK!      do nn=1, nSize
!DOESNTWORK!         coord = coord + nuCoord
!DOESNTWORK!         ! Respectiving periodic BC
!DOESNTWORK!         do cc = 1, SIZE(coord)
!DOESNTWORK!            if (coord(cc) > dataShape(cc)) then
!DOESNTWORK!               coord(cc) = 1
!DOESNTWORK!            end if
!DOESNTWORK!         end do
!DOESNTWORK!         amat = data(coord(1), coord(2), coord(3), coord(4), nu, :, :)
!DOESNTWORK!         call MultiplyMatMat(bmat, dmat, amat)
!DOESNTWORK!         dmat = bmat
!DOESNTWORK!      end do
!DOESNTWORK!      call MultiplyMatMat(bmat, cmat, conjg(transpose(dmat)))
!DOESNTWORK!      cmat = bmat
!DOESNTWORK!      ! Do the links backwards in nu
!DOESNTWORK!      call RealTraceMultMatMat(loopVal, cmat, Ident)
!DOESNTWORK!
!DOESNTWORK!    end function WxNMuNuCoord

                           subroutine plaquette(data, muStart, muEnd, nuEnd, sumTrP, nP, time)
                              complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(in) :: data
                              integer, intent(in) :: muStart, muEnd, nuEnd
                              real(kind=WP), intent(out) :: sumTrP, time
                              integer, intent(out) :: nP
                              !"""
                              !Calculates the plaquette over muStart to muEnd
                              !data is [nt, nx, ny, nz, mu, colour, colour] complex
                              !the plaquette over all lattice is muStart=1, muEnd=4, nuEnd=4
                              !the spatial plaquette is muStart=2, muEnd=4, nuEnd=4
                              !the temporal plaquette is muStart=1, muEnd=1, nuEnd=4
                              !returns the sum of plaquettes, number of plaquettes measured,
                              !the average plaquette and the time taken to calculate it
                              !"""
                              integer, dimension(7) :: dataShape
                              integer, dimension(4) :: muCoord, nuCoord, coordBase, coord
                              integer :: mu, nu, nx, ny, nz, nt, cc  ! Counters
                              ! For intermediate calculating plaquette
                              complex(kind=WC), dimension(3, 3) :: Umu_x, Unu_xmu, UmuUnu
                              complex(kind=WC), dimension(3, 3) :: Umu_xnu, Unu_x, UmudagUnudag
                              real(kind=WP) :: P
                              ! Timers
                              type(watchtype) :: watch
                              real(Kind=SP) :: watchtime
                              call create_watch(watch)
                              call start_watch(watch)
                              dataShape = SHAPE(data)
                              !# hold the sum
                              sumTrP = 0.0_WP
                              !# hold the number measured
                              nP = 0
                              ! write(*,*) data(1,1,1,1,1,1,1)
                              do mu = muStart, muEnd
                                 muCoord(:) = 0
                                 !# This is the shift in mu
                                 muCoord(mu) = 1
                                 do nu = mu + 1, nuEnd
                                    nuCoord(:) = 0
                                    !# This is the shift in nu
                                    nuCoord(nu) = 1
                                    !# loop over all sites
                                    do nx = 1, dataShape(2)
                                       do ny = 1, dataShape(3)
                                          do nz = 1, dataShape(4)
                                             do nt = 1, dataShape(1)
                                                !# U_mu(x)
                                                coordBase = (/nt, nx, ny, nz/)
                                                coord = coordBase
                                                Umu_x = data(coord(1), coord(2), coord(3), coord(4), mu, :, :)
                                                if (ALL(coordBase == (/1, 1, 1, 1/))) then
                                                   !   write(*,*) 'Umu_x', mu, Umu_x(1,1)
                                                end if
                                                !# U_nu(x + amu)
                                                coord = coordBase + muCoord
                                                !# respect periodic boundary conditions
                                                do cc = 1, SIZE(coord)
                                                   if (coord(cc) > dataShape(cc)) then
                                                      coord(cc) = 1
                                                   end if
                                                end do
                                                Unu_xmu = data(coord(1), coord(2), coord(3), coord(4), nu, :, :)
                                                if (ALL(coordBase == (/1, 1, 1, 1/))) then
                                                   !  write(*,*) 'Unu_xmu', nu, Unu_xmu(1,1)
                                                end if
                                                !# U_mu(x + anu)
                                                coord = coordBase + nuCoord
                                                do cc = 1, SIZE(coord)
                                                   if (coord(cc) > dataShape(cc)) then
                                                      coord(cc) = 1
                                                   end if
                                                end do
                                                Umu_xnu = data(coord(1), coord(2), coord(3), coord(4), mu, :, :)
                                                if (ALL(coordBase == (/1, 1, 1, 1/))) then
                                                   !write(*,*) 'Umu_xnu', mu, Umu_xnu(1,1)
                                                end if
                                                !# U_nu(x)
                                                coord = coordBase
                                                Unu_x = data(coord(1), coord(2), coord(3), coord(4), nu, :, :)
                                                if (ALL(coordBase == (/1, 1, 1, 1/))) then
                                                   ! write(*,*) 'Unu_x', nu, Unu_x(1,1)
                                                end if
                                                !# Multiply bottom, right together
                                                call MultiplyMatMat(UmuUnu, Umu_x, Unu_xmu)
                                                !# Multiply left, top together, take dagger
                                                call MultiplyMatdagMatdag(UmudagUnudag, Umu_xnu, Unu_x)
                                                !# multiply two halves together, take trace
                                                call RealTraceMultMatMat(P, UmuUnu, UmudagUnudag)
                                                !if (nx == 1 .and. ny == 1 .and. nz == 1 .and. nt == 1 .and. mu == 2 .and. nu == 3) then
                                                !write(*,*) 'plaquette 9P', 9.0 * P
                                                !end if
                                                sumTrP = sumTrP + P
                                                nP = nP + 1
                                             end do
                                          end do
                                       end do
                                    end do
                                 end do
                              end do
                              call stop_watch(watch)
                              call read_watch(watchTime, watch, 'wall')
                              time = 0.0_WP
                              time = watchTime
                           end subroutine plaquette

                           subroutine polyakov(data, sumTrP, nP, time)
                              complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(in) :: data
                              real(kind=WP), intent(out) :: sumTrP, time
                              integer, intent(out) :: nP
                              !"""
                              !Calculates the polyakov loop
                              !data is [nt, nx, ny, nz, mu, colour, colour] complex
                              !returns the sum of polyakov, number of loops measured,
                              !the average polyakov and the time taken to calculate it
                              !"""
                              integer, dimension(7) :: dataShape
                              integer, dimension(4) :: coord
                              integer :: nx, ny, nz, nt  ! Counters
                              complex(kind=WC), dimension(3, 3) :: Unu_x, Unu_xt, Unu_Temp
                              real(kind=WP) :: P
                              ! Timers
                              real(kind=WP) :: start, end
                              call CPU_TIME(start)
                              dataShape = SHAPE(data)
                              !# hold the sum
                              sumTrP = 0.0_WP
                              !# hold the number measured
                              nP = 0
                              !# loop over all sites
                              do nx = 1, dataShape(2)
                                 do ny = 1, dataShape(3)
                                    do nz = 1, dataShape(4)
                                       ! Get first link U_t(x, 1)
                                       coord = (/1, nx, ny, nz/)
                                       Unu_x = data(coord(1), coord(2), coord(3), coord(4), 1, :, :)
                                       do nt = 2, dataShape(1) - 1
                                          ! Get middle links U_t(x, t)
                                          coord = (/nt, nx, ny, nz/)
                                          Unu_xt = data(coord(1), coord(2), coord(3), coord(4), 1, :, :)
                                          call MultiplyMatMat(Unu_Temp, Unu_x, Unu_xt)
                                          Unu_x = Unu_Temp
                                       end do
                                       ! get final link U_t(x, NT)
                                       coord = (/dataShape(1), nx, ny, nz/)
                                       Unu_x = data(coord(1), coord(2), coord(3), coord(4), 1, :, :)
                                       ! Multiply and trace
                                       call RealTraceMultMatMat(P, Unu_Temp, Unu_x)
                                       call MultiplyMatMat(Unu_xt, Unu_Temp, Unu_x)
                                       sumTrP = sumTrP + P
                                       nP = nP + 1
                                    end do
                                 end do
                              end do
                              call CPU_TIME(end)
                              time = end - start
                           end subroutine polyakov

                           end module FLUE_wloops
