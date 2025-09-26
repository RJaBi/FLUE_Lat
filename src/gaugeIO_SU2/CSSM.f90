module FLUE_SU2_CSSM
  use FLUE_constants, only: WP, WC
  implicit none(external)
  private

contains

  subroutine WriteGaugeField_SU2_CSSM(filename, NX, NY, NZ, NT, U_xd, nfig, beta)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: NX, NY, NZ, NT
      complex(kind=WC), dimension(NT, NX, NY, NZ, 4, 2, 2), intent(in) :: U_xd
      integer, intent(in) :: nfig
      real(kind=WP), intent(in) :: beta
      ! for writing
      complex(kind=WC), dimension(NX, NY, NZ, NT, 4, 2, 2) :: UWrite
      real(kind=WP), , dimension(NX, NY, NZ, NT, 4, 2, 2) :: UWReal, UWImag
      !complex(kind=WC), dimension(4) :: UTmp
      integer, parameter :: infl = 107
      ! counters
      integer :: it, ix, iy, iz
      ! plaq like variables
      real(kind=WP) :: lastPlaq, plaqbarAvg, uzero


      ! Re-order links
      do ix=1, NX
         do iy=1, NY
            do iz=1, NZ
               do it=1, NT
                  UWrite(ix,iy,iz,it,:,:,:) = U_xd(it,ix,iy,iz,:,:,:)
               end do
            end do
         end do
      end do
      ! Cshift mu
      UWrite = CSHIFT(UWrite, 1, dim=5)
      UWReal = real(UWrite, kind=WP)
      UWImag = aimag(UWrite)
      open (infl, file=TRIM(filename), form="unformatted", access="stream", &
           status="replace", action="write", convert="big_endian")
      write(infl) nfig, beta, NX, NY, NZ, NT
      ! write links
      write(infl) UWReal
      write(infl) UWImag
      write(infl) lastPlaq, plabarAvg, uzero
      close(infl)
    end subroutine WriteGaugeField_SU2_CSSM


end module FLUE_SU2_CSSM
