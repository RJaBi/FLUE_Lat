!! Calculate, spatial, temporal plaquettes
!! Compiles for use with python using
!! 'f2py -c fortPlaq.f90 -m fortPlaq'
!!
!! Which can then be imported into python as
!! import fortPlaq
!! fPlaq = fortPlaq.plaq.plaquette
!! which would give you the plaquette function

module FLUE_wloops
   use FLUE_constants, only: WP, WC
   use FLUE_SU3MatrixOps, only: MultiplyMatMat, MultiplyMatdagMatdag, RealTraceMultMatMat
   implicit none(external)
   public

contains

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
      real(kind=WP) :: start, end
      call CPU_TIME(start)
      dataShape = SHAPE(data)
      !# hold the sum
      sumTrP = 0.0_WP
      !# hold the number measured
      nP = 0
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
                        !# U_nu(x + amu)
                        coord = coordBase + muCoord
                        !# respect periodic boundary conditions
                        do cc = 1, SIZE(coord)
                           if (coord(cc) > dataShape(cc)) then
                              coord(cc) = 1
                           end if
                        end do
                        Unu_xmu = data(coord(1), coord(2), coord(3), coord(4), nu, :, :)
                        !# U_mu(x + anu)
                        coord = coordBase + nuCoord
                        do cc = 1, SIZE(coord)
                           if (coord(cc) > dataShape(cc)) then
                              coord(cc) = 1
                           end if
                        end do
                        Umu_xnu = data(coord(1), coord(2), coord(3), coord(4), mu, :, :)
                        !# U_nu(x)
                        coord = coordBase
                        Unu_x = data(coord(1), coord(2), coord(3), coord(4), nu, :, :)
                        !# Multiply bottom, right together
                        call MultiplyMatMat(UmuUnu, Umu_x, Unu_xmu)
                        !# Multiply left, top together, take dagger
                        call MultiplyMatdagMatdag(UmudagUnudag, Umu_xnu, Unu_x)
                        !# multiply two halves together, take trace
                        call RealTraceMultMatMat(P, UmuUnu, UmudagUnudag)
                        sumTrP = sumTrP + P
                        nP = nP + 1
                     end do
                  end do
               end do
            end do
         end do
      end do
      call CPU_TIME(end)
      time = end - start
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
