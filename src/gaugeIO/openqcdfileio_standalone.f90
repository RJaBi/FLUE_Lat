!! Functions to read openQCD format gaugefields
module FLUE_openQCDFileIO_SA
   use FLUE_constants, only: WP, WC
   use FLUE_ILDG_bin, only: FixSU3Matrix
   implicit none(external)
   private
   public :: ReadGaugeField_OpenQCD

contains

  !! Maps an integer a to the set of integers [1,b] i.e. positive integers with cycle length b.
   elemental function modc(a, b) result(c)
    !! Taken directly from COLA
      integer, intent(in) :: a, b
      integer :: c
      !c = a - ((a-1)/b)*b
      c = MODULO(a - 1, b) + 1
   end function modc

!  function determinant(matrix) result(det)
!    complex(c_double_complex), dimension(3,3), intent(in) :: matrix
!    complex(c_double_complex) :: det
!
!    det = matrix(1,1)*(matrix(2,2)*matrix(3,3) - matrix(2,3)*matrix(3,2)) - &
!         matrix(1,2)*(matrix(2,1)*matrix(3,3) - matrix(2,3)*matrix(3,1)) + &
!         matrix(1,3)*(matrix(2,1)*matrix(3,2) - matrix(2,2)*matrix(3,1))
!  end function determinant

   function ReadGaugeField_OpenQCD(filename, NX, NY, NZ, NT) result(U_xd)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: NX, NY, NZ, NT
      complex(kind=WC), dimension(NT, NX, NY, NZ, 4, 3, 3) :: U_xd
      !complex(kind=WP), dimension(:,:,:,:,:,:,:), allocatable :: U_xd

      complex(kind=WC), dimension(3, 3) :: UTmp
      integer, parameter :: infl = 107
      ! Header info
      real(kind=WP) :: plaq
      integer :: ntdim, nxdim, nydim, nzdim
      ! counters
      integer :: it, ix, iy, iz, mu, id
      integer :: jx, jy, jz, jt
      integer, dimension(4) :: dmu
      write (*, *) 'here ', TRIM(filename)
      open (infl, file=TRIM(filename), form='unformatted', access='stream', status='old', action='read', convert='little_endian')
      read (infl) ntdim, nxdim, nydim, nzdim, plaq

      !allocate(U_xd(nxdim,nydim,nzdim,ntdim,4,3,3))

      ! z varies quickest, then y, then x, then t
      do it = 1, ntdim; do ix = 1, nxdim; do iy = 1, nydim; do iz = 1, nzdim
                  if (MODULO(ix + iy + iz + it - 4, 2) == 0) cycle  ! Format only considers odd points

                  do id = 1, 4
                     mu = modc(id - 1, 4)  ! Time dimension first: mu = 4, 1, 2, 3

                     dmu(:) = 0
                     dmu(mu) = 1

                     ! Get the backward site under periodic boundary conditions
                     jx = modc(ix - dmu(1), nxdim)
                     jy = modc(iy - dmu(2), nydim)
                     jz = modc(iz - dmu(3), nzdim)
                     jt = modc(it - dmu(4), ntdim)
                     ! Read the forward and backward links in mu direction
                     !read(infl) U_g(mu,ix,iy,iz,it)%cl(:,:)
                     !read(infl) U_g(mu,jx,jy,jz,jt)%cl(:,:)
                     read (infl) UTmp
                     U_xd(it, ix, iy, iz, mu, :, :) = UTmp
                     read (infl) UTmp
                     U_xd(jt, jx, jy, jz, mu, :, :) = UTmp

                     U_xd(it, ix, iy, iz, mu, :, :) = TRANSPOSE(U_xd(it, ix, iy, iz, mu, :, :))
                     U_xd(jt, jx, jy, jz, mu, :, :) = TRANSPOSE(U_xd(jt, jx, jy, jz, mu, :, :))
                     ! Welll FixSU3Matrix did nothing to the average plaquette value
                     call FixSU3Matrix(U_xd(it, ix, iy, iz, mu, :, :))
                     call FixSU3Matrix(U_xd(jt, jx, jy, jz, mu, :, :))
                  end do
               end do; end do; end do; end do

      close (infl)

      U_xd = CSHIFT(U_xd, -1, dim=5)

   end function ReadGaugeField_OpenQCD

end module FLUE_openQCDFileIO_SA
