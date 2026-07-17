!! Function read and write 'cssm' binary data format as from cola
!! Reads gaugefields and also cola output gauge transforms
!!
!! Essentially, these often have all the
!!        do ic = 1, nc-1
!!          read (101) ReU(:,:,:,:,:,ic,:)
!!          read (101) ImU(:,:,:,:,:,ic,:)
!!       end do
!! Whereas ILDG would have
!!
!!       do it=1,nlt
!!          read(101,rec=it) U_lxd(:,:,:,:,:,:,it)
!!       end do
!!

module FLUE_CSSM_bin
  use FLUE_constants, only: WP, WC, C_INT
  use FLUE_endianIO, only: need_endian_swap, endianSwap, endian_swap_i32
   use FLUE_SU3MatrixOps, only: FixSU3Matrix, orthogonalise_vectors, vector_product
   implicit none(type, external)
   private
   public :: ReadGaugeField_CSSM
   public :: ReadGaugeTransformation_cola

contains

   function ReadGaugeTransformation_cola(filename, NX, NY, NZ, NT) result(G_x)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: NX, NY, NZ, NT
      complex(kind=WC), dimension(3, 3, NT, NX, NY, NZ) :: G_x

      real(kind=WP), dimension(:, :, :, :, :, :), allocatable :: ReG, ImG
      integer :: ix, iy, iz, it, ic, irank
      integer :: jx, jy, jz, jt

      complex(kind=WC), dimension(:, :, :, :, :, :), allocatable :: G_tr
      !type(colour_vector) :: v1,v2,v3
      complex(kind=WC), dimension(3) :: v1, v2, v3

      integer, parameter :: nc = 3
      logical :: needEndianSwap

      allocate (G_tr(nx, ny, nz, nt, nc, nc - 1))
      allocate (ReG(nx, ny, nz, nt, nc, nc - 1))
      allocate (ImG(nx, ny, nz, nt, nc, nc - 1))
      open (101, file=filename, form='unformatted', status='old', action='read')
      ! File format means we read the first two rows of G_x.
      ! Instead we read the first two columns of G_tr (the transpose of G_x) for better memory alignment.
      do ic = 1, nc - 1
         read (101) ReG(:, :, :, :, :, ic)
         read (101) ImG(:, :, :, :, :, ic)
      end do
      close (101)

      ! check if we need to swap endianess
      ! .true. means data is big
      needEndianSwap = need_endian_swap(.true.)
      if (needEndianSwap) then
         ReG = endianSwap(ReG)
         ImG = endianSwap(ImG)
      end if
      
      G_tr(1:nx, 1:ny, 1:nz, 1:nt, :, :) = CMPLX(ReG(1:nx, 1:ny, 1:nz, 1:nt, :, :), ImG(1:nx, 1:ny, 1:nz, 1:nt, :, :), kind=WC)
      do it = 1, nt
         do iz = 1, nz
            do iy = 1, ny
               do ix = 1, nx
                  v1 = G_tr(ix, iy, iz, it, :, 1)
                  v2 = G_tr(ix, iy, iz, it, :, 2)
                  call orthogonalise_vectors(v2, v1)
                  call vector_product(v3, v1, v2)
                  ! G_x is the transpose of G_tr
                  G_x(1, :, it, ix, iy, iz) = v1
                  G_x(2, :, it, ix, iy, iz) = v2
                  G_x(3, :, it, ix, iy, iz) = v3
               end do
            end do
         end do
      end do
      deallocate (G_tr, ReG, ImG)

   end function ReadGaugeTransformation_cola

   function ReadGaugeField_CSSM(filename, NX, NY, NZ, NT, fixSU3) result(U_xd)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: NX, NY, NZ, NT
      logical, optional, intent(in) :: fixSU3
      complex(kind=WC), dimension(3, 3, 4, NT, NX, NY, NZ) :: U_xd
      !
      integer, parameter :: infl = 101
      ! header
      integer, parameter :: i32 = SELECTED_INT_KIND(9)  ! 32 bit integer.
      integer(i32) :: nconfig, nxdim, nydim, nzdim, ntdim
      integer, parameter :: dp = KIND(1.0D0) !! Double precision real scalars.
      real(kind=dP) :: beta
      ! for loading
      real(kind=wP), dimension(:, :, :, :, :, :, :), allocatable :: ReU, ImU
      ! counters
      integer :: ic, mu, it, iz, ix, iy
      ! For reconstructing from two rows
      complex(kind=WC), dimension(3) :: v1, v2, v3
      ! endian
      logical :: needEndianSwap

      open (infl, file=filename, form='unformatted', status='old', action='read')
      read (infl) nconfig, beta, nxdim, nydim, nzdim, ntdim
      !write(*,*) nconfig, beta, nxdim, nydim, nzdim, ntdim

      ! check if we need to swap endianess
      ! .true. means data is big
      needEndianSwap = need_endian_swap(.true.)
      
      allocate (ReU(NX, NY, NZ, NT, 4, 2, 3))
      allocate (ImU(NX, NY, NZ, NT, 4, 2, 3))

      do ic = 1, 3 - 1
         read (infl) ReU(:, :, :, :, :, ic, :)
         read (infl) ImU(:, :, :, :, :, ic, :)
      end do

      if (needEndianSwap) then
         nconfig = endian_swap_i32(nconfig)
         beta = endianSwap(beta)
         nxdim = endian_swap_i32(nxdim)
         nydim = endian_swap_i32(nydim)
         nzdim = endian_swap_i32(nzdim)
         ntdim = endian_swap_i32(ntdim)
         ReU = endianSwap(ReU)
         ImU = endianSwap(ImU)
      end if
      !read(infl) lastPlaq, plaqbarAvg, uzero
      close (infl)
      ! Get the two rows of the SU(3) matrix
      ! Reconstruct the third then put it in the gaugefield variable
      do mu = 1, 4
         do it = 1, NT
            do iz = 1, NZ
               do iy = 1, NY
                  do ix = 1, NX
                     v1 = CMPLX(ReU(ix, iy, iz, it, mu, 1, :), ImU(ix, iy, iz, it, mu, 1, :), kind=WC)
                     v2 = CMPLX(ReU(ix, iy, iz, it, mu, 2, :), ImU(ix, iy, iz, it, mu, 2, :), kind=WC)
                     call orthogonalise_vectors(v2, v1)
                     call vector_product(v3, v1, v2)
                     U_xd(1, :, mu, it, ix, iy, iz) = v1
                     U_xd(2, :, mu, it, ix, iy, iz) = v2
                     U_xd(3, :, mu, it, ix, iy, iz) = v3
                  end do
               end do
            end do
         end do
      end do
      U_xd = CSHIFT(U_xd, -1, dim=5)
   end function ReadGaugeField_CSSM

end module FLUE_CSSM_bin
