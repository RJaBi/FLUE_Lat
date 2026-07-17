!! Functions to read & write SU2 HKLS format
module FLUE_SU2_HKLS
   use, intrinsic :: ISO_C_BINDING, only: C_INT, C_DOUBLE_COMPLEX, C_LONG
   use FLUE_constants, only: WP, WC
   use FLUE_endianIO, only: need_endian_swap, endianSwap, endian_swap_c_double_complex
   !use FLUE_SU3MatrixOps, only: FixSU3Matrix
   !use FLUE_wloops, only: genPlaquette
   implicit none(type, external)
   private
   public :: ReadGaugeField_HKLS
   public :: writeGaugeField_HKLS

contains

   subroutine ReadGaugeField_HKLS(filename, NX, NY, NZ, NT, U_xd, seed, old_nproc, bigEndian)
      character(len=*), intent(IN) :: filename
      integer, intent(IN) :: NX, NY, NZ, NT
      complex(kind=C_DOUBLE_COMPLEX), dimension(2, 2, 4, NT, NX, NY, NZ), intent(OUT) :: U_xd
      integer(kind=C_INT), intent(OUT) :: old_nproc
      integer(kind=C_LONG), dimension(:), allocatable, intent(OUT) :: seed
      logical, intent(IN), optional :: bigEndian
      ! for reading
      complex(kind=C_DOUBLE_COMPLEX), dimension(NX, NY, NZ, NT, 4, 2) :: UHold
      logical :: littleEndian
      !complex(kind=C_DOUBLE_COMPLEX), dimension(4) :: UTmp
      complex(kind=C_DOUBLE_COMPLEX) :: UTmp
      integer, parameter :: infl = 107
      ! counters
      integer :: it, ix, iy, iz, mu, ab, aa, bb
      integer :: ioStatus
      ! endian
      logical :: needEndianSwap
      
      if (PRESENT(bigEndian)) then
         littleEndian = .NOT. bigEndian
      else
         littleEndian = .TRUE.
      end if
      if (littleEndian) then
         ! .true. means little data
         needEndianSwap = need_endian_swap(.true.)
      else
         ! .false. means big data
         needEndianSwap = need_endian_swap(.false.)
      end if
      ! Open in host
      open (infl, file=TRIM(filename), form="unformatted", access="stream", status="old", action="read")
      
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

      allocate (seed(old_nproc))
      read (infl, iostat=IOStatus) seed
      if (IOStatus /= 0) then
         ! Fortran produced
         ! rewind
         rewind (infl)
         ! make it length 1 instead
         deallocate (seed)
         allocate (seed(1))
         read (infl) seed
      end if
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
      if (needEndianSwap) then
         U_xd = endian_swap_c_double_complex(U_xd)
         old_nproc = endianSwap(old_nproc)
         seed = endianSwap(seed)
      end if
   end subroutine ReadGaugeField_HKLS

   subroutine WriteGaugeField_HKLS(filename, NX, NY, NZ, NT, U_xd, seed, old_nproc)
      character(len=*), intent(IN) :: filename
      integer, intent(IN) :: NX, NY, NZ, NT
      complex(kind=C_DOUBLE_COMPLEX), dimension(2, 2, 4, NT, NX, NY, NZ), intent(IN) :: U_xd
      integer(kind=C_INT), intent(IN) :: old_nproc
      integer(kind=C_LONG), dimension(:), intent(IN) :: seed
      ! for writing
      complex(kind=C_DOUBLE_COMPLEX), dimension(NX, NY, NZ, NT, 4, 2) :: UHold
      complex(kind=C_DOUBLE_COMPLEX) :: UTmp
      integer, parameter :: infl = 107
      ! counters
      integer :: it, ix, iy, iz, mu, ab, aa, bb
      ! endian
      logical :: needEndianSwap
      integer(kind=C_INT) :: old_nproc_little
      integer(kind=C_LONG), dimension(size(seed)) :: seed_little

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
      UHold = CSHIFT(UHold, 1, dim=5)
      open (infl, file=TRIM(filename), form="unformatted", access="stream", &
            status="replace", action="write")

      ! check if need to swap endian
      ! .true. means we want little
      needEndianSwap = need_endian_swap(.true.)
      if (needEndianSwap) then
         old_nproc_little = endianSwap(old_nproc)
         seed_little = endianSwap(seed)
         UHold = endian_swap_c_double_complex(UHold)
      else
         old_nproc_little = old_nproc
         seed_little = seed
      end if
      
      write (infl) old_nproc_little
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
