!! Calculate, spatial, temporal plaquettes


module FLUE_SU2_wloops
  use FLUE_constants, only: WP, WC, SP
  use FLUE_wloops, only: periodCoord
  use M_stopwatch, only: watchtype, create_watch, start_watch, stop_watch, destroy_watch, read_watch
  implicit none(type, external)

  complex(kind=WC), dimension(2, 2), parameter :: Ident_SU2 = RESHAPE(source=[ &
       (1.0_WP, 0.0_WP), (0.0_WP, 0.0_WP), &
       (0.0_WP, 0.0_WP), (1.0_WP, 0.0_WP)], &
       shape=[2, 2])
  public
contains

  function SU2_genericPath(data, coordBase, path) result(U_xd)
      complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(in) :: data
      integer, dimension(4), intent(in) :: coordBase
      integer, dimension(:), intent(in) :: path
      complex(kind=WC), dimension(2, 2) :: U_xd
      ! internal
      complex(kind=WC), dimension(2, 2) :: amat, bmat
      integer, dimension(4) :: coord, pCoord
      integer, dimension(7) :: dataShape
      ! counters
      integer :: pp
      integer :: mu
      coord = coordBase
      U_xd = Ident_SU2
      ! Do this to prevent the array temporary in periodCoord
      dataShape = SHAPE(data)
      do pp = 1, SIZE(path)
         bmat = U_xd
         mu = path(pp)
         pCoord = 0
         pCoord(ABS(mu)) = 1
         if (mu < 0) then
            ! the step is backwards so subtract off from coord first
            coord = coord - pCoord
            coord = periodCoord(coord, datashape)
            ! get the link here going forward in mu
            ! and dagger it
            amat = CONJG(TRANSPOSE(data(:, :, ABS(mu), coord(1), coord(2), coord(3), coord(4))))
         else
            amat = data(:, :, mu, coord(1), coord(2), coord(3), coord(4))
            coord = coord + pCoord
            coord = periodCoord(coord, datashape)
         end if
         ! now multiply it into U_xd from the right
         U_xd = matmul(bmat, amat)
      end do
    end function SU2_genericPath

   subroutine SU2_genPlaquette(data, NT, NX, NY, NZ, muStart, muEnd, nuEnd, sumTrP, nP, time)
      complex(kind=WC), dimension(2, 2, 4, NT, NX, NY, NZ), intent(in) :: data
      integer, intent(in) :: muStart, muEnd, nuEnd
      integer, intent(in) :: NT, NX, NY, NZ
      real(kind=WP), intent(out) :: sumTrP, time
      integer, intent(out) :: nP
      ! Counters
      integer :: mu, nu, nnx, nny, nnz, nnt
      ! other variables
      complex(kind=WC), dimension(2, 2) :: plaq
      integer, dimension(7) :: dataShape
      integer, dimension(4) :: plaqPath, coordBase
      real(kind=WP) :: P
      ! Timers
      type(watchtype) :: watch
      real(Kind=SP) :: watchtime
      call create_watch(watch)
      call start_watch(watch)
      dataShape = (/2, 2, 4, NT, NX, NY, NZ/)
      !# hold the sum
      sumTrP = 0.0_WP
      !# hold the number measured
      nP = 0
      do mu = muStart, muEnd
         do nu = mu + 1, nuEnd
            plaqPath = (/mu, nu, -mu, -nu/)
            do nnx = 1, nx
               do nny = 1, ny
                  do nnz = 1, nz
                     do nnt = 1, nt
                        coordBase = (/nnt, nnx, nny, nnz/)
                        plaq = SU2_genericPath(data, coordBase, plaqPath)
                        P = real(plaq(1,1) + plaq(2,2), kind=WP)
                        ! Account for colours
                        P = 0.5_WP * P
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
    end subroutine SU2_genPlaquette

  end module FLUE_SU2_wloops
