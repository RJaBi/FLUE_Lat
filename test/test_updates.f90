module test_updates
   use FLUE_constants, only: WP, WC
   use FLUE_heatbath, only: updateLinks, build_colour_sites
   use FLUE_SU2_heatbath, only: SU2_updateLinks, constructXMatrix
   use Philox, only: C64
   use stdlib_random, only: dist_rand
   use test_helpers, only: assert_close_mat2, &
                           assert_is_unitary2, assert_is_unitary3, &
                           assert_det_one2, assert_det_one3, seed_rng_fixed, &
                           fill_identity_su2, fill_identity_su3
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none(type, external)
   private

   public :: collect_updates

contains

   !=========================================================
   ! Collect tests
   !=========================================================
   subroutine collect_updates(testsuite)
      type(unittest_type), allocatable, intent(OUT) :: testsuite(:)

      testsuite = (/ &
                  new_unittest("updateLinks_preserves_su3", test_updateLinks_preserves_su3), &
                  new_unittest("su2_updateLinks_preserves_su2", test_su2_updateLinks_preserves_su2), &
                  new_unittest("constructXMatrix_returns_su2", test_constructXMatrix_returns_su2)/)
   end subroutine collect_updates

   !=========================================================
   ! Batch invariant checks
   !=========================================================
   subroutine assert_all_links_su3(error, U, tol, message)
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC), intent(IN) :: U(:, :, :, :, :, :, :)
      real(WP), intent(IN) :: tol
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

   subroutine assert_all_links_su2(error, U, tol, message)
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC), intent(IN) :: U(:, :, :, :, :, :, :)
      real(WP), intent(IN) :: tol
      character(len=*), intent(IN) :: message
      integer :: mu, nt, nx, ny, nz

      do mu = 1, SIZE(U, 3)
         do nt = 1, SIZE(U, 4)
            do nx = 1, SIZE(U, 5)
               do ny = 1, SIZE(U, 6)
                  do nz = 1, SIZE(U, 7)
                     call assert_is_unitary2(error, U(:, :, mu, nt, nx, ny, nz), tol, TRIM(message)//": unitarity")
                     if (ALLOCATED(error)) return
                     call assert_det_one2(error, U(:, :, mu, nt, nx, ny, nz), tol, TRIM(message)//": determinant")
                     if (ALLOCATED(error)) return
                  end do
               end do
            end do
         end do
      end do
   end subroutine assert_all_links_su2

   !=========================================================
   ! 1) SU(3) update keeps all links in SU(3)
   !=========================================================
   subroutine test_updateLinks_preserves_su3(error)
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC) :: U(3, 3, 4, 2, 2, 2, 2)
      complex(WC) :: UUpdated(3, 3, 4, 2, 2, 2, 2)
      real(WP), parameter :: beta = 6.0_WP
      real(WP), parameter :: xi = 1.0_WP
      real(WP), parameter :: tol = 1.0E-10_WP
      integer(kind=C64), dimension(2), parameter :: key = (/1_C64, 2_C64/)
      ! Sites linearisation
      integer, allocatable, dimension(:, :, :) :: sites_t, sites_x, sites_y, sites_z
      integer, allocatable, dimension(:, :) :: counts
      call fill_identity_su3(U)
      UUpdated = U

      call build_colour_sites(2, 2, 2, 2, .TRUE., sites_t, sites_x, sites_y, sites_z, counts)

      call updateLinks(U, beta, UUpdated, key, 1, sites_t, sites_x, sites_y, sites_z, counts, 'Wilon')

      call assert_all_links_su3(error, UUpdated, tol, "updateLinks should keep all links in SU(3)")
   end subroutine test_updateLinks_preserves_su3

   !=========================================================
   ! 3) SU(2) update keeps all links in SU(2)
   !=========================================================
   subroutine test_su2_updateLinks_preserves_su2(error)
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC) :: U(2, 2, 4, 2, 2, 2, 2)
      complex(WC) :: UUpdated(2, 2, 4, 2, 2, 2, 2)
      real(WP), parameter :: beta = 2.3_WP
      real(WP), parameter :: tol = 1.0E-10_WP
      integer(kind=C64), dimension(2), parameter :: key = (/1_C64, 2_C64/)

      call fill_identity_su2(U)
      UUpdated = U

      call SU2_updateLinks(U, beta, UUpdated, key, 1)

      call assert_all_links_su2(error, UUpdated, tol, "SU2_updateLinks should keep all links in SU(2)")
   end subroutine test_su2_updateLinks_preserves_su2

   !=========================================================
   ! 4) constructXMatrix(alpha,beta) should return SU(2)
   !=========================================================
   subroutine test_constructXMatrix_returns_su2(error)
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC) :: X(2, 2)
      real(WP), parameter :: tol = 1.0E-12_WP
      real(WP), dimension(4) :: alpha_vals, beta_vals
      integer :: i

      integer(kind=C64), dimension(4) :: counter
      integer(kind=C64), dimension(2) :: key
      integer(kind=C64) :: randInt

      alpha_vals = [0.10_WP, 0.50_WP, 1.00_WP, 2.00_WP]
      beta_vals = [1.50_WP, 2.30_WP, 4.00_WP, 6.00_WP]

      call seed_rng_fixed(124)

      do i = 1, SIZE(alpha_vals)

         randInt = dist_rand(C64)
         counter(1) = randInt
         randInt = dist_rand(C64)
         counter(2) = randInt
         randInt = dist_rand(C64)
         counter(3) = randInt
         randInt = dist_rand(C64)
         counter(4) = randInt
         randInt = dist_rand(C64)
         key(1) = randInt
         randInt = dist_rand(C64)
         key(2) = randInt

         X = constructXMatrix(alpha_vals(i), beta_vals(i), key, counter, i)
         call assert_is_unitary2(error, X, tol, "constructXMatrix should return a unitary SU(2) matrix")
         if (ALLOCATED(error)) return
         call assert_det_one2(error, X, tol, "constructXMatrix should return determinant one")
         if (ALLOCATED(error)) return
      end do
   end subroutine test_constructXMatrix_returns_su2

end module test_updates
