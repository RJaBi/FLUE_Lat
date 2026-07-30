module FLUE_SU2_HKLS
  !< Functions to read and write HKLS format for SU2
  !< As used by [su2hmc](https://doi.org/10.5281/zenodo.7164406)
   use, intrinsic :: ISO_C_BINDING, only: C_INT, C_DOUBLE_COMPLEX, C_LONG
   use FLUE_constants, only: WP, WC
   !use FLUE_SU3MatrixOps, only: FixSU3Matrix
   !use FLUE_wloops, only: genPlaquette
   implicit none(type, external)
   private
   public :: ReadGaugeField_HKLS
   public :: writeGaugeField_HKLS

contains

  subroutine ReadGaugeField_HKLS(filename, NX, NY, NZ, NT, U_xd, seed, old_nproc, bigEndian)
    !< Read a gaugefield from HKLS format
    !< Supports both old Fortran (single seed) and new C (seed per rank) formats
    !< Supports old cray Fortran format (bigEndian=true)
    character(len=*), intent(IN) :: filename
    !< filename to read from
    integer, intent(IN) :: NX, NY, NZ, NT
    !< lattice dimensions
    complex(kind=C_DOUBLE_COMPLEX), dimension(2, 2, 4, NT, NX, NY, NZ), intent(OUT) :: U_xd
    !< output gaugefield
    integer(kind=C_INT), intent(OUT) :: old_nproc
    !< number of procs in header
    integer(kind=C_LONG), dimension(:), allocatable, intent(OUT) :: seed
    !< seeds - either single number or one for each proc
    logical, intent(IN), optional :: bigEndian
    !< Whether we are reading a big-Endian (cray-era) gaugefield
    complex(kind=C_DOUBLE_COMPLEX), dimension(NX, NY, NZ, NT, 4, 2) :: UHold
    !< For holding as we read
    logical :: littleEndian
    !< whether it is littleEndian or not
    complex(kind=C_DOUBLE_COMPLEX) :: UTmp
    !< Read single complex number a time
    integer, parameter :: infl = 107
    !< file unit
    integer :: it, ix, iy, iz, mu, ab, aa, bb
    !< counters
    integer :: ioStatus
    !< Used to test for multiple-seeds

      if (PRESENT(bigEndian)) then
         littleEndian = .NOT. bigEndian
      else
         littleEndian = .TRUE.
      end if
      if (littleEndian) then
         open (infl, file=TRIM(filename), form="unformatted", access="stream", status="old", action="read", convert="little_endian")
      else
         open (infl, file=TRIM(filename), form="unformatted", access="stream", status="old", action="read", convert="big_endian")
      end if
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
   end subroutine ReadGaugeField_HKLS

   subroutine WriteGaugeField_HKLS(filename, NX, NY, NZ, NT, U_xd, seed, old_nproc)
     !< Writes a gaugefield to HKLS format
     !< supports little-endian only
     !< Does support multiple or single seed
     character(len=*), intent(IN) :: filename
     !< filename to write to
     integer, intent(IN) :: NX, NY, NZ, NT
     !< lattice dimensions
     complex(kind=C_DOUBLE_COMPLEX), dimension(2, 2, 4, NT, NX, NY, NZ), intent(IN) :: U_xd
     !< input, internal representation of gaugefield
     integer(kind=C_INT), intent(IN) :: old_nproc
     !< number of processes for header
     integer(kind=C_LONG), dimension(:), intent(IN) :: seed
     !< seed(s) for header
      ! for writing
     complex(kind=C_DOUBLE_COMPLEX), dimension(NX, NY, NZ, NT, 4, 2) :: UHold
     !< re-ordered gaugefield
     complex(kind=C_DOUBLE_COMPLEX) :: UTmp
     !< write single complex number at a time
     integer, parameter :: infl = 107
     !< file unit
     integer :: it, ix, iy, iz, mu, ab, aa, bb
     !< counters

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
