module FLUE_SU2_CSSM
  use FLUE_constants, only: WP, WC
  use FLUE_endianIO, only: need_endian_swap, endianSwap
   use FLUE_SU2_wloops, only: SU2_genPlaquette
   implicit none(type, external)
   private

   public :: writeGaugeField_SU2_CSSM

contains

   subroutine WriteGaugeField_SU2_CSSM(filename, NX, NY, NZ, NT, U_xd, nfig, beta)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: NX, NY, NZ, NT
      complex(kind=WC), dimension(2, 2, 4, NT, NX, NY, NZ), intent(in) :: U_xd
      integer, intent(in) :: nfig
      real(kind=WP), intent(in) :: beta
      ! for writing
      complex(kind=WC), dimension(NX, NY, NZ, NT, 4, 2, 2) :: UWrite
      real(kind=WP), dimension(NX, NY, NZ, NT, 4, 2, 2) :: UWReal, UWImag
      !complex(kind=WC), dimension(4) :: UTmp
      integer, parameter :: infl = 107
      ! counters
      integer :: it, ix, iy, iz, mu
      ! plaq like variables
      real(kind=WP) :: lastPlaq, plaqbarAvg, uzero
      real(kind=WP) :: time, sumtrp
      integer :: nP
      ! endian
      logical :: needEndianSwap
      integer :: nx_little, ny_little, nz_little, nt_little, nfig_little
      real(kind=WP) :: beta_little

      call SU2_genPlaquette(U_xd, NT, NX, NY, NZ, 1, 4, 4, sumTrP, nP, time)

      lastPlaq = sumtrp / real(nP, kind=WP)
      plaqbarAvg = lastplaq

      uzero = plaqbarAvg**0.25_WP

      ! Re-order links
      do ix = 1, NX
         do iy = 1, NY
            do iz = 1, NZ
               do it = 1, NT
                  do mu = 1, 4
                     UWrite(ix, iy, iz, it, mu, :, :) = U_xd(:, :, mu, it, ix, iy, iz)
                  end do
               end do
            end do
         end do
      end do
      ! Cshift mu
      UWrite = CSHIFT(UWrite, 1, dim=5)
      UWReal = real(UWrite, kind=WP)
      UWImag = AIMAG(UWrite)

      ! check whether need endian swap
      ! .false. means big data
      needEndianSwap = need_endian_swap(.false.)
      if (needEndianSwap) then
         nt_little = endianSwap(nt)
         nx_little = endianSwap(nx)
         ny_little = endianSwap(ny)
         nz_little = endianSwap(nz)
         nfig_little = endianSwap(nfig)
         beta_little = endianSwap(beta)
         UWReal = endianSwap(UWReal)
         UWImag = endianSwap(UWImag)
         lastPlaq = endianSwap(lastPlaq)
         plaqbarAvg = endianSwap(plaqBarAvg)
         uzero = endianSwap(uzero)
      else
         nt_little = nt
         nx_little = nx
         ny_little = ny
         nz_little = nz
         nfig_little = nfig
         beta_little = beta
      end if
      
      open (infl, file=TRIM(filename), form="unformatted", access="stream", &
           status="replace", action="write")
      
      write (infl) nfig_little, beta_little, NX_little, NY_little, NZ_little, NT_little
      ! write links
      write (infl) UWReal
      write (infl) UWImag
      write (infl) lastPlaq, plaqbarAvg, uzero
      close (infl)
   end subroutine WriteGaugeField_SU2_CSSM

end module FLUE_SU2_CSSM
