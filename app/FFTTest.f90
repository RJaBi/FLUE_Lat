program multidim_fft
   use FLUE, only: WP, WC, ReadGaugeField_OpenQCD
   use fftpack, only: zfftf, zffti, zfftb
   use stdlib_math, only: linspace
   implicit none
   ! Gauge field variables
   integer, parameter :: NS = 32, NT = 56
   complex(kind=WC), dimension(NT, NS, NS, NS, 4, 3, 3) :: U
   !(56, 32, 32, 32, 3, 3)
   complex(kind=WC), dimension(NT, NS, NS, NS, 3, 3) :: UTrans
   integer, parameter :: NX = NS, NY = NS, NZ = NS
   integer :: i, j, k, l, m, n, o
   ! FFT variables
   complex(kind=WC), dimension(:), allocatable :: temp
   real(kind=WP), dimension(4*NT + 15) :: wsaveNT
   real(kind=WP), dimension(4*NS + 15) :: wsaveNS

   ! Read the gaugefield in
   U = ReadGaugeField_OpenQCD('/home/ryan/Documents/2024/conf/Gen2L/56x32/Gen2l_56x32n1', &
                              NS, NS, NS, NT)

   ! Read gauge transform in
   open (101, file='/home/ryan/Documents/Git/vortex-utils/transBin/Gen2l_56x32n//Gen2l_56x32n1.trans_2.lime.fort', &
         form='unformatted', action='read', access='stream')
   read (101) Utrans
   close (101)

   ! Apply gauge transform
   do i = 1, NT
      do j = 1, NX
         do k = 1, NY
            do l = 1, NZ
               do m = 1, 4
                  U(i, j, k, l, m, :, :) = matmul(UTrans(i, j, k, l, :, :), matmul(U(i, j, k, l, m, :, :), &
                                                                                   conjg(transpose(UTrans(i, j, k, l, :, :)))))
               end do
            end do
         end do
      end do
   end do

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
