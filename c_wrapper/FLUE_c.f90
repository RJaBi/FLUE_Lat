module FLUE_c
  !! Note that the c_wrapper and python uses order [NT, NX, NY, NZ, mu, colour, colour]
  !! This contrasts the internal ordering of [colour, colour, mu, NT, NX, NY, NZ]
  !! This is for compatibility with i.e. lyncs_io
   use, intrinsic :: ISO_C_BINDING, only: C_DOUBLE, C_INT, C_DOUBLE_COMPLEX, C_CHAR
   !gp!use FLUE, only: Q_Average, get_qhat, cone_cut, MultiplyMatMat, calc_mom_space_scalarD, scalarGluonProp
   use FLUE, only: ReadGaugeField_ILDG, writeGaugeField_ILDG, &
                   ReadGaugeField_OpenQCD, writeGaugeField_OpenQCD, &
                   ReadGaugefield_CSSM, &
                   genplaquette, StoutSmearLinks
   implicit none(external)
   public

contains

   !gaugeIO
   ! Writers
   subroutine writeGaugeField_ILDG_c(filename, U_xd, NX, NY, NZ, NT)
      character(kind=C_CHAR, len=*), intent(in) :: filename
      integer(kind=C_INT), intent(in) :: NX, NY, NZ, NT
      complex(kind=C_DOUBLE_COMPLEX), dimension(NT, NX, NY, NZ, 4, 3, 3), intent(in) :: U_xd
      complex(kind=C_DOUBLE_COMPLEX), dimension(3, 3, 4, NT, NX, NY, NZ) :: uinternal
      integer :: mu, ac, bc
      do ac = 1, 3
         do bc = 1, 3
            do mu = 1, 4
               uinternal(ac, bc, mu, :, :, :, :) = U_xd(:, :, :, :, mu, ac, bc)
            end do
         end do
      end do
      call writeGaugeField_ILDG(filename, uinternal, NX, NY, NZ, NT, .true.)
   end subroutine writeGaugeField_ILDG_c

   subroutine writeGaugeField_OpenQCD_c(filename, U_xd, NX, NY, NZ, NT)
      character(kind=C_CHAR, len=*), intent(in) :: filename
      integer(kind=C_INT), intent(in) :: NX, NY, NZ, NT
      complex(kind=C_DOUBLE_COMPLEX), dimension(NT, NX, NY, NZ, 4, 3, 3), intent(in) :: U_xd
      complex(kind=C_DOUBLE_COMPLEX), dimension(3, 3, 4, NT, NX, NY, NZ) :: uinternal
      integer :: mu, ac, bc
      do ac = 1, 3
         do bc = 1, 3
            do mu = 1, 4
               uinternal(ac, bc, mu, :, :, :, :) = U_xd(:, :, :, :, mu, ac, bc)
            end do
         end do
      end do
      call writeGaugeField_OpenQCD(filename, uinternal, NX, NY, NZ, NT)
   end subroutine writeGaugeField_OpenQCD_c

   ! Readers
   function ReadGaugeField_ILDG_c(filename, NX, NY, NZ, NT) result(U_xd)
      character(kind=C_CHAR, len=*), intent(in) :: filename
      integer(kind=C_INT), intent(in) :: NX, NY, NZ, NT
      complex(kind=C_DOUBLE_COMPLEX), dimension(NT, NX, NY, NZ, 4, 3, 3) :: U_xd
      complex(kind=C_DOUBLE_COMPLEX), dimension(3, 3, 4, NT, NX, NY, NZ) :: uinternal
      integer :: mu, ac, bc
      uinternal = ReadGaugeField_ILDG(filename, NX, NY, NZ, NT, .true.)
      do ac = 1, 3
         do bc = 1, 3
            do mu = 1, 4
               U_xd(:, :, :, :, mu, ac, bc) = uinternal(ac, bc, mu, :, :, :, :)
            end do
         end do
      end do
   end function ReadGaugeField_ILDG_c

   function ReadGaugeField_OpenQCD_c(filename, NX, NY, NZ, NT) result(U_xd)
      character(kind=C_CHAR, len=*), intent(in) :: filename
      integer(kind=C_INT), intent(in) :: NX, NY, NZ, NT
      complex(kind=C_DOUBLE_COMPLEX), dimension(NT, NX, NY, NZ, 4, 3, 3) :: U_xd
      complex(kind=C_DOUBLE_COMPLEX), dimension(3, 3, 4, NT, NX, NY, NZ) :: uinternal
      integer :: mu, ac, bc
      uinternal = ReadGaugeField_OpenQCD(filename, NX, NY, NZ, NT, .true.)
      do ac = 1, 3
         do bc = 1, 3
            do mu = 1, 4
               U_xd(:, :, :, :, mu, ac, bc) = uinternal(ac, bc, mu, :, :, :, :)
            end do
         end do
      end do
   end function ReadGaugeField_OpenQCD_c

   function ReadGaugeField_CSSM_c(filename, NX, NY, NZ, NT) result(U_xd)
      character(kind=C_CHAR, len=*), intent(in) :: filename
      integer(kind=C_INT), intent(in) :: NX, NY, NZ, NT
      complex(kind=C_DOUBLE_COMPLEX), dimension(NT, NX, NY, NZ, 4, 3, 3) :: U_xd
      complex(kind=C_DOUBLE_COMPLEX), dimension(3, 3, 4, NT, NX, NY, NZ) :: uinternal
      integer :: mu, ac, bc
      uinternal = ReadGaugeField_CSSM(filename, NX, NY, NZ, NT, .true.)
      do ac = 1, 3
         do bc = 1, 3
            do mu = 1, 4
               U_xd(:, :, :, :, mu, ac, bc) = uinternal(ac, bc, mu, :, :, :, :)
            end do
         end do
      end do
   end function ReadGaugeField_CSSM_c

   ! wloops
   subroutine genplaquette_c(data, nt, nx, ny, nz, mustart, muend, nuend, sumtrp, np, time)
      complex(kind=C_DOUBLE_COMPLEX), dimension(nt, nx, ny, nz, 4, 3, 3), intent(in) :: data
      integer(kind=C_INT), intent(in) :: mustart, muend, nuend
      integer(kind=C_INT), intent(in) :: nt, nx, ny, nz
      real(kind=C_DOUBLE), intent(out) :: sumtrp, time
      integer(kind=C_INT), intent(out) :: np
      complex(kind=C_DOUBLE_COMPLEX), dimension(3, 3, 4, NT, NX, NY, NZ) :: uinternal
      integer :: mu, ac, bc
      do ac = 1, 3
         do bc = 1, 3
            do mu = 1, 4
               uinternal(ac, bc, mu, :, :, :, :) = data(:, :, :, :, mu, ac, bc)
            end do
         end do
      end do
      call genplaquette(uinternal, nt, nx, ny, nz, mustart, muend, nuend, sumtrp, np, time)
   end subroutine genplaquette_c

   subroutine stoutsmearlinks_c(data, rho, nSweeps, nt, nx, ny, nz, usmeared)
      complex(kind=C_DOUBLE_COMPLEX), dimension(nt, nx, ny, nz, 4, 3, 3), intent(in) :: data
      real(kind=C_DOUBLE), intent(in) :: rho
      integer(kind=C_INT), intent(in) :: nSweeps
      integer(kind=C_INT), intent(in) :: nx, ny, nz, nt
      complex(kind=C_DOUBLE_COMPLEX), dimension(nt, nx, ny, nz, 4, 3, 3), intent(out) :: usmeared
      complex(kind=C_DOUBLE_COMPLEX), dimension(3, 3, 4, NT, NX, NY, NZ) :: uinternal, uinternalsmeared
      integer :: mu, ac, bc
      do ac = 1, 3
         do bc = 1, 3
            do mu = 1, 4
               uinternal(ac, bc, mu, :, :, :, :) = data(:, :, :, :, mu, ac, bc)
            end do
         end do
      end do
      call StoutSmearLinks(uinternal, rho, nsweeps, uinternalsmeared)
      do ac = 1, 3
         do bc = 1, 3
            do mu = 1, 4
               usmeared(:, :, :, :, mu, ac, bc) = uinternalsmeared(ac, bc, mu, :, :, :, :)
            end do
         end do
      end do
   end subroutine stoutsmearlinks_c

end module FLUE_c
