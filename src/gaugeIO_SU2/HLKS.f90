!! Functions to read & write SU2 HKLS format
module FLUE_SU2_HKLS
   use, intrinsic :: ISO_C_BINDING, only: C_INT, C_DOUBLE_COMPLEX, C_LONG
   use FLUE_constants, only: WP, WC
   !use FLUE_SU3MatrixOps, only: FixSU3Matrix
   !use FLUE_wloops, only: genPlaquette
   implicit none(type, external)
   private
   public :: ReadGaugeField_HKLS
   public :: writeGaugeField_HKLS

contains

   subroutine ReadGaugeField_HKLS(filename, NX, NY, NZ, NT, U_xd, seed, old_nproc)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: NX, NY, NZ, NT
      complex(kind=C_DOUBLE_COMPLEX), dimension(2, 2, 4, NT, NX, NY, NZ), intent(out) :: U_xd
      integer(kind=C_INT), intent(out) :: old_nproc
      integer(kind=C_LONG), intent(out) :: seed
      ! for reading
      complex(kind=C_DOUBLE_COMPLEX), dimension(NX, NY, NZ, NT, 4, 2) :: UHold
      !complex(kind=C_DOUBLE_COMPLEX), dimension(4) :: UTmp
      complex(kind=C_DOUBLE_COMPLEX) :: UTmp
      integer, parameter :: infl = 107
      ! counters
      integer :: it, ix, iy, iz, mu, ab, aa, bb

      open (infl, file=TRIM(filename), form="unformatted", access="stream", status="old", action="read", convert="little_endian")

      read (infl) old_nproc
      do ab = 1, 2
         do mu = 1, 4
            do it = 1, NT
               do iz = 1, NZ
                  do iy = 1, NY
                     do ix = 1, NX
                        read (infl) UTmp
                        UHold(ix, iy, iz, it, mu, ab) = UTmp
                     end do
                  end do
               end do
            end do
         end do
      end do
      read (infl) seed
      close (infl)
      ! put it in 2x2
      do it = 1, NT
         do iz = 1, NZ
            do iy = 1, NY
               do ix = 1, NX
                  ! Unpack
                  U_xd(1, 1, :, it, ix, iy, iz) = UHold(ix, iy, iz, it, :, 1)
                  U_xd(2, 2, :, it, ix, iy, iz) = CONJG(UHold(ix, iy, iz, it, :, 1))
                  U_xd(1, 2, :, it, ix, iy, iz) = UHold(ix, iy, iz, it, :, 2)
                  U_xd(2, 1, :, it, ix, iy, iz) = -CONJG(UHold(ix, iy, iz, it, :, 2))

               end do
            end do
         end do
      end do
      ! Shift the mu ordering
      U_xd = CSHIFT(U_xd, -1, dim=3)
   end subroutine ReadGaugeField_HKLS

   subroutine WriteGaugeField_HKLS(filename, NX, NY, NZ, NT, U_xd, seed, old_nproc)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: NX, NY, NZ, NT
      complex(kind=C_DOUBLE_COMPLEX), dimension(2, 2, 4, NT, NX, NY, NZ), intent(in) :: U_xd
      integer(kind=C_INT), intent(in) :: old_nproc
      integer(kind=C_LONG), intent(in) :: seed
      ! for writing
      complex(kind=C_DOUBLE_COMPLEX), dimension(NX, NY, NZ, NT, 4, 2) :: UHold
      complex(kind=C_DOUBLE_COMPLEX) :: UTmp
      integer, parameter :: infl = 107
      ! counters
      integer :: it, ix, iy, iz, mu, ab, aa, bb

      ! put it in 1x2
      do it = 1, NT
         do iz = 1, NZ
            do iy = 1, NY
               do ix = 1, NX
                  UHold(ix, iy, iz, it, :, 1) = U_xd(1, 1, :, it, ix, iy, iz)
                  UHold(ix, iy, iz, it, :, 2) = U_xd(1, 2, :, it, ix, iy, iz)
               end do
            end do
         end do
      end do
      ! Shift the mu ordering
      UHold = CSHIFT(UHold, 1, dim=3)
      open (infl, file=TRIM(filename), form="unformatted", access="stream", &
            status="replace", action="write", convert="little_endian")

      write (infl) old_nproc
      do ab = 1, 2
         do mu = 1, 4
            do it = 1, NT
               do iz = 1, NZ
                  do iy = 1, NY
                     do ix = 1, NX
                        UTmp = UHold(ix, iy, iz, it, mu, ab)
                        write (infl) UTmp
                     end do
                  end do
               end do
            end do
         end do
      end do
      write (infl) seed
      close (infl)
   end subroutine WriteGaugeField_HKLS

end module FLUE_SU2_HKLS
