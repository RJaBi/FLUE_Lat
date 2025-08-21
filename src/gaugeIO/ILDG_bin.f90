!! Functions to read and write ILDG binary data formats as from cola
module FLUE_ILDG_bin
  use FLUE_constants, only: WP, WC
  use FLUE_SU3MatrixOps, only: FixSU3Matrix
   !use stdlib_linalg, only: det
   use, intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT
   implicit none(external)
   private
   public :: ReadGaugeField_ILDG

contains

   function ReadGaugeField_ILDG(filename, NX, NY, NZ, NT, fixSU3) result(U_xd)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: NX, NY, NZ, NT
      logical, optional, intent(in) :: fixSU3
      complex(kind=WC), dimension(NT, NX, NY, NZ, 4, 3, 3) :: U_xd
      complex(kind=WC), dimension(3, 3, 4, NX, NY, NZ, NT) :: URead
      integer, parameter :: infl = 101
      integer :: matrix_len, irecl
      ! counters
      integer :: it, ix, iy, iz, mu, nu
      logical :: fixSU3Set

      if (present(fixSU3)) then
         fixSU3Set = fixSU3
      else
         fixSU3Set = .true.
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
                     U_xd(it, ix, iy, iz, nu, :, :) = TRANSPOSE(URead(:, :, mu, ix, iy, iz, it))
                     if (fixSU3Set) then
                        call FixSU3Matrix(U_xd(it, ix, iy, iz, nu, :, :))
                     end if
                  end do
               end do
            end do
         end do
      end do

   end function ReadGaugeField_ILDG

end module FLUE_ILDG_bin
