module FLUE_SU2_wloops
  !< For wilson-loop (or staple) calculations in SU2
   use FLUE_constants, only: WP, WC, SP
   use FLUE_wloops, only: periodCoord
   use M_stopwatch, only: watchtype, create_watch, start_watch, stop_watch, destroy_watch, read_watch
   implicit none(type, external)

   complex(kind=WC), dimension(2, 2), parameter :: Ident_SU2 = RESHAPE(source=[ &
                                                                       (1.0_WP, 0.0_WP), (0.0_WP, 0.0_WP), &
                                                                       (0.0_WP, 0.0_WP), (1.0_WP, 0.0_WP)], &
                                                                       shape=[2, 2])
   !< 2x2 Identity matrix in complex variable
   public
contains

  pure function SU2_genericPath(data, coordBase, path) result(U_xd)
    !< Returns the links starting from coordbase multiplied together along path
    !< Used to construct wilson loops, staples, etc
    complex(kind=WC), dimension(:, :, :, :, :, :, :), intent(in) :: data
    !< The gaugefield of shape 2,2,4,nt,nx,ny,nz
    integer, dimension(4), intent(in) :: coordBase
    !< The starting coordinate [it, ix, iy, iz]
    integer, dimension(:), intent(in) :: path
    !< The path to follow in terms of +ve/-ve mu/nu/etc
    complex(kind=WC), dimension(2, 2) :: U_xd
    !< The resulting path ordered multiplication of links
      ! internal
    complex(kind=WC), dimension(2, 2) :: amat, bmat
    !< Working SU2 variables
    integer, dimension(4) :: coord, pCoord
    !< coordinate holders for moving around [it, ix, iy, iz]
    integer, dimension(7) :: dataShape
    !< The shape of the data [2,2,4,nt,nx,ny,nz]
    integer :: pp
    !< counter for step along path
    integer :: mu
    !< direction to move in
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
         U_xd = MATMUL(bmat, amat)
      end do
   end function SU2_genericPath

   subroutine SU2_genPlaquette(data, NT, NX, NY, NZ, muStart, muEnd, nuEnd, sumTrP, nP, time)
     !< Calculates the plaquette using genplaquette
     !< Outputs the sum over all real trace plaquettes, the number of plaquettes and the time taken
     !< Do all plaquettes using mustart=1, muend=4, nuend=4
     !< Do spatial plaquettes using mustart=2, muend=4, nuend=4
     !< Do temporal plaquettes using mustart=1, muend=1, nuend=4
     complex(kind=WC), dimension(2, 2, 4, NT, NX, NY, NZ), intent(in) :: data
     !< The gaugefield
     integer, intent(in) :: muStart, muEnd, nuEnd
     !< which dimensions to loop over
     integer, intent(in) :: NT, NX, NY, NZ
     !< The size of the lattice
     real(kind=WP), intent(out) :: sumTrP
     !< Sum of ReTr(P) of all plaquettes
     real(kind=WP), intent(out) :: time
     !< Time taken for calculation
     integer, intent(out) :: nP
     !< Number of plaquettes
     integer :: mu, nu, nnx, nny, nnz, nnt
     !< counters
      ! other variables
     complex(kind=WC), dimension(2, 2) :: plaq
     !< plaquette variable holder
     integer, dimension(7) :: dataShape
     !< Size of the gaugefield
     integer, dimension(4) :: plaqPath, coordBase
     !< Variables for constructing all plaquettes using genericPath_SU2
     real(kind=WP) :: P
     !< real trace of single plaquette
      ! Timers
     type(watchtype) :: watch
     !< holds multiple types of timer (cpu, wall, etc)
     real(Kind=SP) :: watchtime
     !< output wall time
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
                        P = real(plaq(1, 1) + plaq(2, 2), kind=WP)
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
