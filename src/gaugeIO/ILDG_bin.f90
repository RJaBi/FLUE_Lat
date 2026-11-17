
module FLUE_ILDG_bin
  !< Functions to read and write ILDG (binary) format gaugefields
   use, intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT, int64
   use FLUE_constants, only: WP, WC
   use FLUE_SU3MatrixOps, only: FixSU3Matrix
   use FLUE_endianIO, only: need_endian_swap, endianSwap
   implicit none(type, external)
   private
   public :: ReadGaugeField_ILDG
   public :: writeGaugeField_ILDG

   public :: writeGaugeField_ILDG_unit

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

   subroutine writeGaugeField_ILDG_Unit( &
        myUnit, U_xd, NX, NY, NZ, NT, fixSU3)
     !< Writes the binary portion of the ILDG gaugefield to an open unit
     integer, intent(in) :: myUnit
     !< Already-open big-endian unformatted stream unit.
     integer, intent(in) :: NX, NY, NZ, NT
     !< Lattice dimensions.
     complex(kind=WC), dimension(3, 3, 4, NT, NX, NY, NZ), intent(in) :: U_xd
     !< Input gauge field in FLUE indexing
     logical, optional, intent(in) :: fixSU3
     !< Reproject links to SU(3), default false
     complex(kind=WC), dimension(3, 3, 4, NX, NY, NZ, NT) :: UWrite
     !< re-ordered gaugefield to write
     integer, dimension(4), parameter :: Directionmap = (/2, 3, 4, 1/)
     !< For re-ordering
     integer :: it, ix, iy, iz, mu
     !< integers
     logical :: fixSU3Set
     !< Whether we are re-projecting each link to SU3
     fixSU3Set = .false.
     if (present(fixSU3)) fixSU3Set = fixSU3
     ! Reorder directions t,x,y,z -> x,y,z,t and transpose each matrix.
     !
     ! directionMap(mu) returns the corresponding input direction:
     ! output x <- input x (2)
     ! output y <- input y (3)
     ! output z <- input z (4)
     ! output t <- input t (1)
     do concurrent ( &
          it = 1:NT, &
          ix = 1:NX, &
          iy = 1:NY, &
          iz = 1:NZ, &
          mu = 1:4)
        UWrite(:, :, mu, ix, iy, iz, it) = transpose( &
             U_xd(:, :, directionMap(mu), it, ix, iy, iz))
     end do
     if (fixSU3Set) then
        do concurrent(it=1:nt, ix=1:nx, iy=1:ny, iz=1:nz, mu=1:4) default(none) shared(UWrite)
           call FixSU3Matrix(UWrite(:, :, mu, ix, iy, iz, it))
        end do
     end if
     ! Write timeslice by timeslice to avoid 2GB limit
     do it = 1, NT
        write (myUnit) UWrite(:, :, :, :, :, :, it)
     end do
 end subroutine writeGaugeField_ILDG_Unit


   subroutine writeGaugeField_ILDG(filename, U_xd, NX, NY, NZ, NT, fixSU3)
     !< Write a gaugefield in ILDG (binary) format
     character(len=*), intent(IN) :: filename
     !< filename to write to
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
