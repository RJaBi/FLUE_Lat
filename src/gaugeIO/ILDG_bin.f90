
module FLUE_ILDG_bin
  !< Functions to read and write ILDG (binary) format gaugefields
   use, intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT
   use FLUE_constants, only: WP, WC
   use FLUE_SU3MatrixOps, only: FixSU3Matrix
   implicit none(type, external)
   private
   public :: ReadGaugeField_ILDG
   public :: writeGaugeField_ILDG

contains

  function ReadGaugeField_ILDG(filename, NX, NY, NZ, NT, fixSU3) result(U_xd)
    !< Read a gaugefield in ILDG (binary) format
    character(len=*), intent(IN) :: filename
    !< filename to read from
    integer, intent(IN) :: NX, NY, NZ, NT
    !< lattice dimensions
    logical, optional, intent(IN) :: fixSU3
    !< optionally re-project to SU3
    complex(kind=WC), dimension(3, 3, 4, NT, NX, NY, NZ) :: U_xd
    !< Output gaugefield
    complex(kind=WC), dimension(3, 3, 4, NX, NY, NZ, NT) :: URead
    !< gaugefield read in
    integer, parameter :: infl = 101
    !< use this file unit
    integer :: matrix_len, irecl
    !< Step through data records
    integer :: it, ix, iy, iz, mu, nu
    !< counters
    logical :: fixSU3Set
    !< whether we are re-projecting each link

      if (PRESENT(fixSU3)) then
         fixSU3Set = fixSU3
      else
         fixSU3Set = .TRUE.
      end if
      ! First read the gaugefield
      write (OUTPUT_UNIT, *) TRIM(filename)
      matrix_len = 16 * 3 * 3
      irecl = matrix_len * 4 * NX * NY * NZ
      ! write(*,*) matrix_len, irecl, 3, 3, 4, nx, ny, nz, nt
      open (infl, file=TRIM(filename), form='unformatted', access='direct', &
            status='old', action='read', recl=irecl, convert='BIG_ENDIAN')
      do it = 1, NT
         read (infl, rec=it) URead(:, :, :, :, :, :, it)
      end do
      close (infl)
      ! THen really dumbly re-order the indices
      do it = 1, NT
         do ix = 1, NX
            do iy = 1, NY
               do iz = 1, NZ
                  do mu = 1, 4
                     ! Mu is the openqcd index
                     ! i.e. starts t,x,y,z
                     ! nu is the ILDG index
                     ! i.e. starts x,y,z,t
                     select case (mu)
                     case (1)
                        nu = 2
                     case (2)
                        nu = 3
                     case (3)
                        nu = 4
                     case (4)
                        nu = 1
                     end select
                     U_xd(:, :, nu, it, ix, iy, iz) = TRANSPOSE(URead(:, :, mu, ix, iy, iz, it))
                     if (fixSU3Set) then
                        call FixSU3Matrix(U_xd(:, :, nu, it, ix, iy, iz))
                     end if
                  end do
               end do
            end do
         end do
      end do

   end function ReadGaugeField_ILDG

   subroutine writeGaugeField_ILDG(filename, U_xd, NX, NY, NZ, NT, fixSU3)
     !< Write a gaugefield in ILDG (binary) format
     character(len=*), intent(IN) :: filename
     !< filename to read to
     integer, intent(IN) :: NX, NY, NZ, NT
     !< lattice dimensions
     logical, optional, intent(IN) :: fixSU3
     !< optionally re-project to SU3
     complex(kind=WC), dimension(3, 3, 4, NT, NX, NY, NZ), intent(IN) :: U_xd
     !< input gaugefield
     complex(kind=WC), dimension(3, 3, 4, NX, NY, NZ, NT) :: URead
     !< re-ordered gaugefield to write
     integer, parameter :: infl = 101
     !< use this file unit
     integer :: matrix_len, irecl
     !< step through records in data
     integer :: it, ix, iy, iz, mu, nu
     !< counters
     logical :: fixSU3Set
     !< whether we are re-projecting each link to SU3

      if (PRESENT(fixSU3)) then
         fixSU3Set = fixSU3
      else
         fixSU3Set = .TRUE.
      end if

      ! First really dumbly re-order the indices
      do it = 1, NT
         do ix = 1, NX
            do iy = 1, NY
               do iz = 1, NZ
                  do mu = 1, 4
                     ! Mu is the openqcd index
                     ! i.e. starts t,x,y,z
                     ! nu is the ILDG index
                     ! i.e. starts x,y,z,t
                     select case (mu)
                     case (1)
                        nu = 2
                     case (2)
                        nu = 3
                     case (3)
                        nu = 4
                     case (4)
                        nu = 1
                     end select
                     URead(:, :, mu, ix, iy, iz, it) = TRANSPOSE(U_xd(:, :, nu, it, ix, iy, iz))
                     if (fixSU3Set) then
                        call FixSU3Matrix(URead(:, :, mu, ix, iy, iz, it))
                     end if

                  end do
               end do
            end do
         end do
      end do

      ! Then write the gaugefield
      matrix_len = 16 * 3 * 3
      irecl = matrix_len * 4 * NX * NY * NZ
      ! write(*,*) matrix_len, irecl, 3, 3, 4, nx, ny, nz, nt
      open (infl, file=TRIM(filename), form='unformatted', access='direct', &
            status='replace', action='write', recl=irecl, convert='BIG_ENDIAN')
      do it = 1, NT
         write (infl, rec=it) URead(:, :, :, :, :, :, it)
      end do
      close (infl)

   end subroutine WriteGaugeField_ILDG

end module FLUE_ILDG_bin
