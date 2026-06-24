module test_gauge_path
   use FLUE_constants, only: WP, WC
   use FLUE_matrixConstants, only: Ident3x3
   use FLUE_wloops, only: genericPath, periodCoord
   use test_helpers, only: assert_close_mat3, fill_identity_su3
   use testdrive, only: new_unittest, unittest_type, error_type, check
   implicit none(type, external)
   private
   public :: collect_gauge_path

contains

   subroutine collect_gauge_path(testsuite)
      type(unittest_type), allocatable, intent(OUT) :: testsuite(:)
      testsuite = [ &
                  new_unittest("periodcoord_no_wrap_interior", test_periodcoord_no_wrap_interior), &
                  new_unittest("periodcoord_wrap_forward_t", test_periodcoord_wrap_forward_t), &
                  new_unittest("periodcoord_wrap_backward_t", test_periodcoord_wrap_backward_t), &
                  new_unittest("periodcoord_wrap_all_axes", test_periodcoord_wrap_all_axes), &
                  new_unittest("genericpath_identity", test_genericpath_identity), &
                  new_unittest("genericpath_backtracking_su3", test_genericpath_backtracking_su3), &
                  new_unittest("genericpath_two_step_product", test_genericpath_two_step_product), &
                  new_unittest("genericpath_periodic_single_wrap", test_genericpath_periodic_single_wrap), &
                  new_unittest("genericpath_periodic_wrap_backtrack", test_genericpath_periodic_wrap_backtrack) &
                  ]
   end subroutine collect_gauge_path

   pure function datashape_from_field(U) result(ds)
      complex(WC), intent(IN) :: U(:, :, :, :, :, :, :)
      integer :: ds(7)
      ds = SHAPE(U)
   end function datashape_from_field
   !=========================================================
   ! Direct periodCoord tests
   !=========================================================
   subroutine test_periodcoord_no_wrap_interior(error)
    !! Check that periodcoord does not effect the interior points
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC) :: data(3, 3, 4, 2, 2, 2, 2)
      integer :: coord(4), got(4), expected(4), ds(7)
      call fill_identity_su3(data)
      ds = datashape_from_field(data)
      coord = [1, 2, 1, 2]
      expected = [1, 2, 1, 2]
      got = periodCoord(coord, ds)
      call check(error, ALL(got == expected), &
                 "periodCoord should leave interior coordinates unchanged")
   end subroutine test_periodcoord_no_wrap_interior

   subroutine test_periodcoord_wrap_forward_t(error)
    !! Check forward wrapping
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC) :: data(3, 3, 4, 2, 2, 2, 2)
      integer :: coord(4), got(4), expected(4), ds(7)
      call fill_identity_su3(data)
      ds = datashape_from_field(data)
      ! NT = ds(4) = 2, so stepping to 3 should wrap to 1
      coord = [3, 1, 1, 1]
      expected = [1, 1, 1, 1]
      got = periodCoord(coord, ds)
      call check(error, ALL(got == expected), &
                 "periodCoord should wrap forward across the temporal boundary")
   end subroutine test_periodcoord_wrap_forward_t

   subroutine test_periodcoord_wrap_backward_t(error)
    !! Check backwards wrapping
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC) :: data(3, 3, 4, 2, 2, 2, 2)
      integer :: coord(4), got(4), expected(4), ds(7)
      call fill_identity_su3(data)
      ds = datashape_from_field(data)
      ! 0 should wrap to NT = ds(4) = 2
      coord = [0, 1, 1, 1]
      expected = [2, 1, 1, 1]
      got = periodCoord(coord, ds)
      call check(error, ALL(got == expected), &
                 "periodCoord should wrap backward across the temporal boundary")
   end subroutine test_periodcoord_wrap_backward_t

   subroutine test_periodcoord_wrap_all_axes(error)
      ! Do a wrap on all axes
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC) :: data(3, 3, 4, 2, 2, 2, 2)
      integer :: coord(4), got(4), expected(4), ds(7)
      call fill_identity_su3(data)
      ds = datashape_from_field(data)
      ! On a 2x2x2x2 lattice:
      ! nt=3 -> 1, nx=0 -> 2, ny=3 -> 1, nz=0 -> 2
      coord = [3, 0, 3, 0]
      expected = [1, 2, 1, 2]
      got = periodCoord(coord, ds)
      call check(error, ALL(got == expected), &
                 "periodCoord should wrap correctly on all four lattice axes")
   end subroutine test_periodcoord_wrap_all_axes
   !=========================================================
   ! Existing / strengthened genericPath tests
   !=========================================================
   subroutine test_genericpath_identity(error)
    !! Generic path on an identity field multiplies identity together
    !! I.e. should return identity
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC) :: data(3, 3, 4, 2, 2, 2, 2)
      integer :: coord(4), path(4)
      complex(WC) :: P(3, 3)
      real(WP), parameter :: tol = 1.0E-14_WP
      call fill_identity_su3(data)
      coord = [1, 1, 1, 1]
      path = [1, -1, 1, -1]
      P = genericPath(data, coord, path)
      call assert_close_mat3(error, P, Ident3x3, tol, &
                             "genericPath on identity field should return identity")
   end subroutine test_genericpath_identity

   subroutine test_genericpath_backtracking_su3(error)
    !! Make non-identity field
    !! Make sure that going forward then back cancels to give identity
    !! As gauge links are unitary
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC) :: data(3, 3, 4, 2, 2, 2, 2)
      complex(WC) :: A(3, 3), P(3, 3)
      integer :: coord(4), path(2)
      integer :: i, j, k, l
      real(WP), parameter :: tol = 1.0E-12_WP
      call fill_identity_su3(data)
      A = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      A(1, 1) = EXP(CMPLX(0.0_WP, 0.3_WP, kind=WC))
      A(2, 2) = EXP(CMPLX(0.0_WP, -0.3_WP, kind=WC))
      A(3, 3) = CMPLX(1.0_WP, 0.0_WP, kind=WC)
      do i = 1, SIZE(data, 4)
         do j = 1, SIZE(data, 5)
            do k = 1, SIZE(data, 6)
               do l = 1, SIZE(data, 7)
                  data(:, :, 1, i, j, k, l) = A
               end do
            end do
         end do
      end do
      coord = [1, 1, 1, 1]
      path = [1, -1]
      P = genericPath(data, coord, path)
      call assert_close_mat3(error, P, Ident3x3, tol, &
                             "genericPath(U_mu then -U_mu) should backtrack to identity")
   end subroutine test_genericpath_backtracking_su3

   subroutine test_genericpath_two_step_product(error)
    !! Check that going forward is working as intended
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC) :: data(3, 3, 4, 2, 2, 2, 2)
      complex(WC) :: A(3, 3), P(3, 3), ref(3, 3)
      integer :: coord(4), path(2)
      integer :: i, j, k, l
      real(WP), parameter :: tol = 1.0E-12_WP
      call fill_identity_su3(data)
      A = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      A(1, 1) = EXP(CMPLX(0.0_WP, 0.2_WP, kind=WC))
      A(2, 2) = EXP(CMPLX(0.0_WP, -0.2_WP, kind=WC))
      A(3, 3) = CMPLX(1.0_WP, 0.0_WP, kind=WC)
      do i = 1, SIZE(data, 4)
         do j = 1, SIZE(data, 5)
            do k = 1, SIZE(data, 6)
               do l = 1, SIZE(data, 7)
                  data(:, :, 1, i, j, k, l) = A
               end do
            end do
         end do
      end do
      coord = [1, 1, 1, 1]
      path = [1, 1]
      P = genericPath(data, coord, path)
      ref = MATMUL(A, A)
      call assert_close_mat3(error, P, ref, tol, &
                             "genericPath on two identical forward links should equal A*A")
   end subroutine test_genericpath_two_step_product

   subroutine test_genericpath_periodic_single_wrap(error)
    !! Test generic path across boundary forward
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC) :: data(3, 3, 4, 2, 2, 2, 2)
      complex(WC) :: A(3, 3), P(3, 3)
      integer :: coord(4), path(1)
      integer :: i, j, k, l
      real(WP), parameter :: tol = 1.0E-12_WP
      call fill_identity_su3(data)
      A = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      A(1, 1) = EXP(CMPLX(0.0_WP, 0.17_WP, kind=WC))
      A(2, 2) = EXP(CMPLX(0.0_WP, -0.17_WP, kind=WC))
      A(3, 3) = CMPLX(1.0_WP, 0.0_WP, kind=WC)
      do i = 1, SIZE(data, 4)
         do j = 1, SIZE(data, 5)
            do k = 1, SIZE(data, 6)
               do l = 1, SIZE(data, 7)
                  data(:, :, 1, i, j, k, l) = A
               end do
            end do
         end do
      end do
      ! Start at boundary in t-direction on a 2-site lattice and step forward once.
      coord = [2, 1, 1, 1]
      path = [1]
      P = genericPath(data, coord, path)
      call assert_close_mat3(error, P, A, tol, &
                             "genericPath should wrap correctly across the periodic boundary for one forward step")
   end subroutine test_genericpath_periodic_single_wrap

   subroutine test_genericpath_periodic_wrap_backtrack(error)
    !! Test generic path across boundary backwards
      type(error_type), allocatable, intent(OUT) :: error
      complex(WC) :: data(3, 3, 4, 2, 2, 2, 2)
      complex(WC) :: A(3, 3), P(3, 3)
      integer :: coord(4), path(2)
      integer :: i, j, k, l
      real(WP), parameter :: tol = 1.0E-12_WP
      call fill_identity_su3(data)
      A = CMPLX(0.0_WP, 0.0_WP, kind=WC)
      A(1, 1) = EXP(CMPLX(0.0_WP, 0.29_WP, kind=WC))
      A(2, 2) = EXP(CMPLX(0.0_WP, -0.29_WP, kind=WC))
      A(3, 3) = CMPLX(1.0_WP, 0.0_WP, kind=WC)
      do i = 1, SIZE(data, 4)
         do j = 1, SIZE(data, 5)
            do k = 1, SIZE(data, 6)
               do l = 1, SIZE(data, 7)
                  data(:, :, 1, i, j, k, l) = A
               end do
            end do
         end do
      end do
      coord = [2, 1, 1, 1]
      path = [1, -1]
      P = genericPath(data, coord, path)
      call assert_close_mat3(error, P, Ident3x3, tol, &
                             "genericPath should backtrack correctly across a periodic boundary")
   end subroutine test_genericpath_periodic_wrap_backtrack
end module test_gauge_path
