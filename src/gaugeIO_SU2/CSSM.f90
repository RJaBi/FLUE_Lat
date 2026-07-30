module FLUE_SU2_CSSM
  !< Functions to write CSSM format for SU2
  !< As used by COLA-based code
   use FLUE_constants, only: WP, WC
   use FLUE_SU2_wloops, only: SU2_genPlaquette
   implicit none(type, external)
   private

   public :: writeGaugeField_SU2_CSSM

contains

  subroutine WriteGaugeField_SU2_CSSM(filename, NX, NY, NZ, NT, U_xd, nfig, beta)
    !< Write a gaugefield to CSSM format (from cola)
    character(len=*), intent(in) :: filename
    !< filename to write to
    integer, intent(in) :: NX, NY, NZ, NT
    !< lattice dimensions
    complex(kind=WC), dimension(2, 2, 4, NT, NX, NY, NZ), intent(in) :: U_xd
    !< input, internal format gaugefield to write
    integer, intent(in) :: nfig
    !< config number to write
    real(kind=WP), intent(in) :: beta
    !< beta value to write
    ! for writing
    complex(kind=WC), dimension(NX, NY, NZ, NT, 4, 2, 2) :: UWrite
    !< re-ordered gaugefield
    real(kind=WP), dimension(NX, NY, NZ, NT, 4, 2, 2) :: UWReal, UWImag
    !< real and imaginary part of re-ordered gaugefield
      !complex(kind=WC), dimension(4) :: UTmp
    integer, parameter :: infl = 107
    !< file unit
    integer :: it, ix, iy, iz, mu
    !< counter
    real(kind=WP) :: lastPlaq, plaqbarAvg, uzero
    !< average plaquette, average plaquette, u_0 (average link)
    real(kind=WP) :: time, sumtrp
    !< for calculating plaquette
    integer :: nP
    !< number of plaquettes

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
      open (infl, file=TRIM(filename), form="unformatted", access="stream", &
            status="replace", action="write", convert="big_endian")
      write (infl) nfig, beta, NX, NY, NZ, NT
      ! write links
      write (infl) UWReal
      write (infl) UWImag
      write (infl) lastPlaq, plaqbarAvg, uzero
      close (infl)
   end subroutine WriteGaugeField_SU2_CSSM

end module FLUE_SU2_CSSM
