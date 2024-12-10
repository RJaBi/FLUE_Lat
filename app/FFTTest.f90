program multidim_fft
   use FLUE, only: WP, WC, ReadGaugeField_OpenQCD, RealTraceMat, plaquette, &
                   ReadGaugeField_ILDG, ReadGaugeTransformation_cola
   use fftpack, only: zfftf, zffti, zfftb
   use stdlib_math, only: linspace
   !use stdlib_linalg, only: det
   use stdlib_io_npy, only: save_npy
   implicit none(external)
   ! Gauge field variables
   integer, parameter :: NS = 32, NT = 64
   complex(kind=WC), dimension(NT, NS, NS, NS, 4, 3, 3) :: U
   !(56, 32, 32, 32, 3, 3)
   complex(kind=WC), dimension(NT, NS, NS, NS, 3, 3) :: UTrans
   integer, parameter :: NX = NS, NY = NS, NZ = NS
   integer :: i, j, k, l, m, n, o, cc
   integer, dimension(4) :: muCoord, coord
   integer, dimension(4), parameter :: dataShape = (/NT, NX, NY, NZ/)
   real(kind=WP) :: sLink, tLink
   ! FFT variables
   complex(kind=WC), dimension(:), allocatable :: temp
   real(kind=WP), dimension(4*NT + 15) :: wsaveNT
   real(kind=WP), dimension(4*NS + 15) :: wsaveNS
   ! plaquette params
   real(kind=WP) :: sumTrP, time, plaq
   integer :: nP

   ! Read the gaugefield in
!   U = ReadGaugeField_OpenQCD('/home/ryan/Documents/2024/conf/Gen2L/56x32/Gen2l_56x32n1', &
!        NS, NS, NS, NT)
   U = ReadGaugeField_ILDG('/home/ryan/Documents/2024/conf/PACS-CS/RCNF2+1/RC32x64_B1900Kud01370000Ks013640' &
                           //'00C1715/RC32x64_B1900Kud01370000Ks01364000C1715-b-002510', &
                           NS, NS, NS, NT)

   !call RealTraceMat(sLink, U(1,1,1,1,2,:,:))
   !call RealTraceMat(tLink, U(1,1,1,1,1,:,:))
   !write(*,*) 'space', (1.0_WP/3.0_WP) * sLink, det(U(1,1,1,1,1,:,:))
   !write(*,*) 'time', (1.0_WP/3.0_WP) * tLink, det(U(1,1,1,1,2,:,:))
   call plaquette(U, 1, 4, 4, sumTrP, nP, time)
   plaq = 9.0_WP * sumTrp / nP
   write (*, *) 'U Plaquette for', ' is ', plaq, 'and took', time, 'seconds'
   write (*, *) 'u0bar is (1/3*plaq)**0.25', ((1.0_WP / 3.0_WP) * plaq)**0.25_WP

   write (*, *)
   write (*, *) "AFTER THIS IS GAUGE-FIXED"
   write (*, *)

   ! Read gauge transform in
   UTrans = ReadGaugeTransformation_cola('/home/ryan/Documents/2024/conf/PACS-CS/RCNF2+1/' &
                                         //'RC32x64_B1900Kud01370000Ks01364000C1715/' &
                                         //'RC32x64_B1900Kud01370000Ks01364000C1715-b-002510' &
                                         //'FFTA.fix.gag', &
                                         NS, NS, NS, NT)

   ! Read gauge transform in
   !open (101, file='/home/ryan/Documents/Git/vortex-utils/transBin/Gen2l_56x32n//Gen2l_56x32n1.trans_2.lime.fort', &
   !      form='unformatted', action='read', access='stream')
   !read (101) Utrans
   !close (101)

   ! Apply gauge transform
   do m = 1, 4
      muCoord(:) = 0
      muCoord(m) = 1
      do i = 1, NT
         do j = 1, NX
            do k = 1, NY
               do l = 1, NZ
                  coord = (/i, j, k, l/) + muCoord
                  !# respect periodic boundary conditions
                  do cc = 1, SIZE(coord)
                     if (coord(cc) > dataShape(cc)) then
                        coord(cc) = 1
                     end if
                  end do

                  U(i, j, k, l, m, :, :) = MATMUL(UTrans(i, j, k, l, :, :), &
                                                  MATMUL(U(i, j, k, l, m, :, :), &
                                                         CONJG(TRANSPOSE(UTrans(coord(1), coord(2), coord(3), coord(4), :, :)))))
               end do
            end do
         end do
      end do
   end do

   !call RealTraceMat(sLink, U(1,1,1,1,2,:,:))
   !call RealTraceMat(tLink, U(1,1,1,1,1,:,:))
   !write(*,*) 'space', (1.0_WP/3.0_WP) * sLink, det(U(1,1,1,1,1,:,:))
   !write(*,*) 'time', (1.0_WP/3.0_WP) * tLink, det(U(1,1,1,1,2,:,:))
   call plaquette(U, 1, 4, 4, sumTrP, nP, time)
   plaq = 9.0_WP * sumTrp / nP
   write (*, *) 'U Plaquette for', ' is ', plaq, 'and took', time, 'seconds'

   !stop

   tLink = 0.0
   do m = 1, 4
      do i = 1, NT
         do j = 1, NX
            do k = 1, NY
               do l = 1, NZ
                  !call RealTraceMat(sLink, U(i,j,k,l,m,:,:))
                  tLink = tLink + real(U(i, j, k, l, m, 1, 1) + U(i, j, k, l, m, 2, 2) + U(i, j, k, l, m, 3, 3))
               end do
            end do
         end do
      end do
   end do
   tLink = tLink / real(4 * NT * NX * NY * NZ, kind=WP)
   write (*, *) 'u0bar is (1/N)TrU_mu(x) = ', tLink
   stop

   ! setup working storage
   call zffti(NT, wsaveNT)
   call zffti(NS, wsaveNS)

   ! Perform FFT over the first dimension
   allocate (temp(NT))
   do j = 1, NX
      do k = 1, NY
         do l = 1, NZ
            do m = 1, 4
               do n = 1, 3
                  do o = 1, 3
                     temp = U(:, j, k, l, m, n, o)
                     call zfftf(NT, temp, wsaveNT)
                     U(:, j, k, l, m, n, o) = temp
                  end do
               end do
            end do
         end do
      end do
   end do
   deallocate (temp)

   ! Perform FFT over the second dimension
   allocate (temp(NX))
   do i = 1, NT
      do k = 1, NY
         do l = 1, NZ
            do m = 1, 4
               do n = 1, 3
                  do o = 1, 3
                     temp = U(i, :, k, l, m, n, o)
                     call zfftf(NS, temp, wsaveNS)
                     U(i, :, k, l, m, n, o) = temp
                  end do
               end do
            end do
         end do
      end do
   end do
   deallocate (temp)

   ! Perform FFT over the third dimension
   allocate (temp(NY))
   do i = 1, NT
      do j = 1, NX
         do l = 1, NZ
            do m = 1, 4
               do n = 1, 3
                  do o = 1, 3
                     temp = U(i, j, :, l, m, n, o)
                     call zfftf(NS, temp, wsaveNS)
                     U(i, j, :, l, m, n, o) = temp
                  end do
               end do
            end do
         end do
      end do
   end do
   deallocate (temp)

   ! Perform FFT over the fourth dimension
   allocate (temp(NZ))
   do i = 1, NT
      do j = 1, NX
         do k = 1, NY
            do m = 1, 4
               do n = 1, 3
                  do o = 1, 3
                     temp = U(i, j, k, :, m, n, o)
                     call zfftf(NS, temp, wsaveNS)
                     U(i, j, k, :, m, n, o) = temp
                  end do
               end do
            end do
         end do
      end do
   end do
   deallocate (temp)
   write (*, *) U

end program multidim_fft
