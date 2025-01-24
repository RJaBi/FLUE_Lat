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
   use FLUE_SU3MatrixOps, only: MultiplyMatMat, MultiplyMatdagMatdag, RealTraceMultMatMat, Ident
   implicit none(external)
   public

contains

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
