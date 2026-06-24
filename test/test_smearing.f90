module test_smearing
   use FLUE_constants, only: WP, WC
   use FLUE_matrixConstants, only: Ident3x3
   use FLUE_stoutSmearing, only: StoutSmearLinks
   use test_helpers, only: assert_close_mat3, &
                           assert_is_unitary3, assert_det_one3, fill_identity_su3
   use testdrive, only: new_unittest, unittest_type, error_type, check

   implicit none(type, external)
   private
   public :: collect_smearing

contains

   !=========================================================
   ! Collect tests
   !=========================================================
   subroutine collect_smearing(testsuite)
      type(unittest_type), allocatable, intent(OUT) :: testsuite(:)
      testsuite = [ &
                  new_unittest("stout_identity_preserved", test_stout_identity_preserved), &
                  new_unittest("stout_preserves_su3", test_stout_preserves_su3), &
                  new_unittest("stout_zero_sweeps_is_identity", test_stout_zero_sweeps_is_identity) &
                  ]
   end subroutine collect_smearing
   !=========================================================
   ! Fillers
   !=========================================================
   subroutine fill_commuting_su3(U, theta)
      ! A simple exact SU(3) background:
      ! A = diag(exp(i theta), exp(-i theta), 1)
      complex(kind=WC), intent(OUT) :: U(:, :, :, :, :, :, :)
      real(kind=WP), intent(IN) :: theta
      complex(kind=WC) :: A(3, 3)
      integer :: mu, nt, nx, ny, nz

      A = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      A(1, 1) = EXP(CMPLX(0.0_WP, theta, kind=WC))
      A(2, 2) = EXP(CMPLX(0.0_WP, -theta, kind=WC))
      A(3, 3) = CMPLX(1.0_WP, 0.0_WP, kind=WC)

      do concurrent(mu=1:SIZE(U, 3), nt=1:SIZE(U, 4), nx=1:SIZE(U, 5), &
                    ny=1:SIZE(U, 6), nz=1:SIZE(U, 7))
         U(:, :, mu, nt, nx, ny, nz) = A
      end do
   end subroutine fill_commuting_su3

   subroutine assert_all_links_su3(error, U, tol, message)
      type(error_type), allocatable, intent(OUT) :: error
      complex(kind=WC), intent(IN) :: U(:, :, :, :, :, :, :)
      real(kind=WP), intent(IN) :: tol
      character(len=*), intent(IN) :: message
      integer :: mu, nt, nx, ny, nz

      do mu = 1, SIZE(U, 3)
         do nt = 1, SIZE(U, 4)
            do nx = 1, SIZE(U, 5)
               do ny = 1, SIZE(U, 6)
                  do nz = 1, SIZE(U, 7)
                     call assert_is_unitary3(error, U(:, :, mu, nt, nx, ny, nz), tol, TRIM(message)//": unitarity")
                     if (ALLOCATED(error)) return
                     call assert_det_one3(error, U(:, :, mu, nt, nx, ny, nz), tol, TRIM(message)//": determinant")
                     if (ALLOCATED(error)) return
                  end do
               end do
            end do
         end do
      end do
   end subroutine assert_all_links_su3

   !=========================================================
   ! 1) Identity field should stay identity under stout smearing
   !=========================================================
   subroutine test_stout_identity_preserved(error)
      type(error_type), allocatable, intent(OUT) :: error
      complex(kind=WC) :: U(3, 3, 4, 2, 2, 2, 2)
      complex(kind=WC) :: USmear(3, 3, 4, 2, 2, 2, 2)
      integer :: mu, nt, nx, ny, nz
      real(kind=WP), parameter :: rho = 0.10_WP
      integer, parameter :: nSweeps = 1
      real(kind=WP), parameter :: tol = 1.0E-13_WP

      call fill_identity_su3(U)
      call StoutSmearLinks(U, rho, nSweeps, USmear)

      do mu = 1, SIZE(USmear, 3)
         do nt = 1, SIZE(USmear, 4)
            do nx = 1, SIZE(USmear, 5)
               do ny = 1, SIZE(USmear, 6)
                  do nz = 1, SIZE(USmear, 7)
                     call assert_close_mat3(error, USmear(:, :, mu, nt, nx, ny, nz), Ident3x3, tol, &
                                            "StoutSmearLinks should preserve the identity field")
                     if (ALLOCATED(error)) return
                  end do
               end do
            end do
         end do
      end do
   end subroutine test_stout_identity_preserved

   !=========================================================
   ! 2) Smearing output should remain in SU(3)
   !=========================================================
   subroutine test_stout_preserves_su3(error)
      type(error_type), allocatable, intent(OUT) :: error
      complex(kind=WC) :: U(3, 3, 4, 2, 2, 2, 2)
      complex(kind=WC) :: USmear(3, 3, 4, 2, 2, 2, 2)
      real(kind=WP), parameter :: rho = 0.10_WP
      integer, parameter :: nSweeps = 2
      real(kind=WP), parameter :: tol = 1.0E-10_WP

      call fill_commuting_su3(U, 0.23_WP)
      call StoutSmearLinks(U, rho, nSweeps, USmear)

      call assert_all_links_su3(error, USmear, tol, "StoutSmearLinks should preserve SU(3)")
   end subroutine test_stout_preserves_su3

   !=========================================================
   ! 3) Zero sweeps should return the input unchanged
   !=========================================================
   subroutine test_stout_zero_sweeps_is_identity(error)
      type(error_type), allocatable, intent(OUT) :: error
      complex(kind=WC) :: U(3, 3, 4, 2, 2, 2, 2)
      complex(kind=WC) :: USmear(3, 3, 4, 2, 2, 2, 2)
      real(kind=WP), parameter :: rho = 0.10_WP
      integer, parameter :: nSweeps = 0
      real(kind=WP), parameter :: tol = 1.0E-14_WP

      call fill_commuting_su3(U, 0.31_WP)
      call StoutSmearLinks(U, rho, nSweeps, USmear)

      call assert_close_field_su3(error, USmear, U, tol, &
                                  "StoutSmearLinks with nSweeps=0 should return the input unchanged")
   end subroutine test_stout_zero_sweeps_is_identity

   subroutine assert_close_field_su3(error, A, B, tol, message)
      type(error_type), allocatable, intent(OUT) :: error
      complex(kind=WC), intent(IN) :: A(:, :, :, :, :, :, :), B(:, :, :, :, :, :, :)
      real(kind=WP), intent(IN) :: tol
      character(len=*), intent(IN) :: message

      call testdrive_check_close_field(error, MAXVAL(ABS(A - B)) < tol, message)
   end subroutine assert_close_field_su3

   subroutine testdrive_check_close_field(error, ok, message)
      type(error_type), allocatable, intent(OUT) :: error
      logical, intent(IN) :: ok
      character(len=*), intent(IN) :: message
      call check(error, ok, message)
   end subroutine testdrive_check_close_field

end module test_smearing
