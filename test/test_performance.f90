module test_performance
   use, intrinsic :: IEEE_ARITHMETIC, only: ieee_is_finite
   use FLUE_constants, only: WP, WC, SP
   use FLUE_heatbath, only: updateLinks, build_colour_sites
   use FLUE_matrixConstants, only: Ident3x3
   use FLUE_openQCDFileIO_SA, only: ReadGaugeField_OpenQCD
   use FLUE_wloops, only: genPlaquette
   use M_stopwatch, only: watchtype, create_watch, start_watch, stop_watch, &
                          destroy_watch, read_watch
   use Philox, only: C64
   use test_helpers, only: seed_rng_fixed, fill_identity_su3
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use tomlf, only: toml_table, toml_error, toml_load, get_value
   implicit none(type, external)
   private
   public :: collect_performance

contains

   subroutine collect_performance(testsuite)
      type(unittest_type), allocatable, intent(OUT) :: testsuite(:)

      testsuite = [ &
                  new_unittest("openqcd_performance", test_openqcd_performance), &
                  new_unittest("heatbath_wilson_xi1", test_heatbath_wilson_xi1), &
                  new_unittest("heatbath_wilson_xi10", test_heatbath_wilson_xi10), &
                  new_unittest("heatbath_symanzik_xi1", test_heatbath_symanzik_xi1), &
                  new_unittest("heatbath_symanzik_xi10", test_heatbath_symanzik_xi10) &
                  ]
   end subroutine collect_performance

   !=========================================================
   ! Existing OpenQCD read/plaquette benchmark
   !=========================================================
   subroutine test_openqcd_performance(error)
      type(error_type), allocatable, intent(OUT) :: error
      complex(kind=WC), allocatable :: U(:, :, :, :, :, :, :)
      type(watchtype) :: watch
      real(kind=WP) :: read_time, plaquette_time, baseline_read, baseline_plaq, tolerance
      real(kind=WP) :: sumtrp
      real(kind=SP) :: watchtime
      character(len=128) :: gauge_file, baseline_path
      integer :: NT, NS, np
      logical :: found_baseline
      NT = 8
      NS = 24
      gauge_file = 'testdata/Gen2_8x24n9'
      baseline_path = 'testdata/reference_values.toml'
      tolerance = 0.30_WP
      call read_openqcd_baseline(baseline_path, baseline_read, baseline_plaq, found_baseline)
      call check(error, found_baseline, &
                 'performance baseline section [performance.Gen2_8x24n9] must be present in '//TRIM(baseline_path))
      if (ALLOCATED(error)) return
      call create_watch(watch)
      call start_watch(watch)
      U = ReadGaugeField_OpenQCD(TRIM(gauge_file), NS, NS, NS, NT)
      call stop_watch(watch)
      call read_watch(watchtime, watch, 'wall')
      call destroy_watch(watch)
      read_time = real(watchtime, kind=WP)
      call genPlaquette(U, NT, NS, NS, NS, 1, 4, 4, sumtrp, np, plaquette_time)
      write (*, '(A,F10.6)') 'Measured ReadGaugeField_OpenQCD time: ', read_time
      write (*, '(A,F10.6)') 'Measured genPlaquette time:          ', plaquette_time
      call check(error, read_time <= baseline_read * (1.0_WP + tolerance), &
                 'ReadGaugeField_OpenQCD time within baseline tolerance')
      if (ALLOCATED(error)) return
      call check(error, plaquette_time <= baseline_plaq * (1.0_WP + tolerance), &
                 'genPlaquette time within baseline tolerance')
   end subroutine test_openqcd_performance
   !=========================================================
   ! Heatbath benchmark cases
   !=========================================================
   subroutine test_heatbath_wilson_xi1(error)
      type(error_type), allocatable, intent(OUT) :: error
      call benchmark_heatbath_case(error, "wilson_xi1", "wilson", 1.0_WP, 6.0_WP)
   end subroutine test_heatbath_wilson_xi1

   subroutine test_heatbath_wilson_xi10(error)
      type(error_type), allocatable, intent(OUT) :: error
      call benchmark_heatbath_case(error, "wilson_xi10", "wilson", 10.0_WP, 6.8_WP)
   end subroutine test_heatbath_wilson_xi10

   subroutine test_heatbath_symanzik_xi1(error)
      type(error_type), allocatable, intent(OUT) :: error
      call benchmark_heatbath_case(error, "symanzik_xi1", "symanzik", 1.0_WP, 6.0_WP)
   end subroutine test_heatbath_symanzik_xi1

   subroutine test_heatbath_symanzik_xi10(error)
      type(error_type), allocatable, intent(OUT) :: error
      call benchmark_heatbath_case(error, "symanzik_xi10", "symanzik", 10.0_WP, 6.8_WP)
   end subroutine test_heatbath_symanzik_xi10

   subroutine benchmark_heatbath_case(error, case_name, action_tag, xi, beta)
      type(error_type), allocatable, intent(OUT) :: error
      character(len=*), intent(IN) :: case_name, action_tag
      real(WP), intent(IN) :: xi, beta
      integer, parameter :: NS = 4
      integer, parameter :: NT = 8
      integer, parameter :: nTraj = 8
      integer, parameter :: nRepeat = 3
      real(WP), parameter :: tolerance = 0.30_WP
      complex(WC) :: U(3, 3, 4, NT, NS, NS, NS)
      complex(WC) :: UNew(3, 3, 4, NT, NS, NS, NS)
      real(WP) :: total_time(nRepeat), time_per_traj(nRepeat), median_t
      real(WP) :: baseline_update
      logical :: found_baseline
      integer :: irep, iTraj
      real(WP) :: aplaq, splaq, tplaq
      type(watchtype) :: watch
      real(kind=SP) :: watchtime
      real(WP) :: sumTrP, plaq_time
      integer :: nPlaq
      ! Sites linearisation
      integer, allocatable, dimension(:, :, :) :: sites_t, sites_x, sites_y, sites_z
      integer, allocatable, dimension(:, :) :: counts
      do irep = 1, nRepeat
         call fill_identity_su3(U)
         UNew = U
         ! Fixed stdlib_random seed for repeatable benchmark trajectories
         !CALL seed_rng_fixed(20240615 + irep)
         ! Warm-up call to reduce first-call effects in the timed region

         call build_colour_sites(NT, NS, NS, NS, .TRUE., sites_t, sites_x, sites_y, sites_z, counts)

         !CALL updateLinks(U, beta, xi, UNew, trim(action_tag))
         call updateLinks(U, beta, UNew, (/1_C64, 2_C64/), 0, sites_t, sites_x, sites_y, sites_z, counts, TRIM(action_tag))
         U = UNew

         call create_watch(watch)
         call start_watch(watch)
         do iTraj = 1, nTraj
            !CALL updateLinks(U, beta, xi, UNew, trim(action_tag))
            call updateLinks(U, beta, UNew, (/1_C64, 2_C64/), iTraj, sites_t, sites_x, sites_y, sites_z, counts, TRIM(action_tag))
            U = UNew
         end do
         call stop_watch(watch)
         call read_watch(watchtime, watch, 'wall')
         call destroy_watch(watch)
         total_time(irep) = real(watchtime, kind=WP)
         time_per_traj(irep) = total_time(irep) / real(nTraj, kind=WP)
      end do
      median_t = median3(time_per_traj)
      ! Final plaquette sanity values from the last run
      call genPlaquette(U, NT, NS, NS, NS, 1, 4, 4, sumTrP, nPlaq, plaq_time)
      aplaq = sumTrP / (3.0_WP * real(nPlaq, kind=WP))

      if (xi /= 1.0_WP) then
         call genPlaquette(U, NT, NS, NS, NS, 2, 4, 4, sumTrP, nPlaq, plaq_time)
         splaq = sumTrP / (3.0_WP * real(nPlaq, kind=WP))

         call genPlaquette(U, NT, NS, NS, NS, 1, 1, 4, sumTrP, nPlaq, plaq_time)
         tplaq = sumTrP / (3.0_WP * real(nPlaq, kind=WP))
      else
         splaq = 0.0_WP
         tplaq = 0.0_WP
      end if

      call read_heatbath_baseline("testdata/reference_values.toml", case_name, baseline_update, found_baseline)
      call check(error, found_baseline, &
                 "heatbath performance baseline [performance.heatbath."//TRIM(case_name) &
                 //"] must exist in testdata/reference_values.toml")
      if (ALLOCATED(error)) return

      !write(*,'(A,1X,A)')        'Benchmark case:', trim(case_name)
      !write(*,'(A,1X,A)')        'Action tag:   ', trim(action_tag)
      !write(*,'(A,F10.6)')       'xi:           ', xi
      !write(*,'(A,F10.6)')       'beta:         ', beta
      !write(*,'(A,3(1X,F10.6))') 'time/traj:    ', time_per_traj
      !write(*,'(A,F10.6)')       'median/traj:  ', median_t
      !write(*,'(A,F10.6)')       'baseline:     ', baseline_update
      !write(*,'(A,F10.6)')       'aplaq:        ', aplaq
      !if (xi /= 1.0_WP) then
      !   write(*,'(A,F10.6)')    'splaq:        ', splaq
      !   write(*,'(A,F10.6)')    'tplaq:        ', tplaq
      !end if

      call check(error, median_t <= baseline_update * (1.0_WP + tolerance), &
                 "heatbath update time per trajectory within baseline tolerance")
      if (ALLOCATED(error)) return

      call check(error, ieee_is_finite(aplaq) .AND. ABS(aplaq) <= 1.0_WP + 1.0E-12_WP, &
                 "average plaquette should be finite and bounded")
      if (ALLOCATED(error)) return

      if (xi /= 1.0_WP) then
         call check(error, ieee_is_finite(splaq) .AND. ABS(splaq) <= 1.0_WP + 1.0E-12_WP, &
                    "spatial plaquette should be finite and bounded")
         if (ALLOCATED(error)) return

         call check(error, ieee_is_finite(tplaq) .AND. ABS(tplaq) <= 1.0_WP + 1.0E-12_WP, &
                    "temporal plaquette should be finite and bounded")
      end if
   end subroutine benchmark_heatbath_case

   !=========================================================
   ! Helpers
   !=========================================================
   pure real(WP) function median3(x) result(m)
      real(WP), intent(IN) :: x(3)
      real(WP) :: a, b, c
      a = x(1)
      b = x(2)
      c = x(3)
      if ((a <= b .AND. b <= c) .OR. (c <= b .AND. b <= a)) then
         m = b
      else if ((b <= a .AND. a <= c) .OR. (c <= a .AND. a <= b)) then
         m = a
      else
         m = c
      end if
   end function median3

   subroutine read_openqcd_baseline(filename, read_baseline, plaquette_baseline, found)
      character(len=*), intent(IN) :: filename
      real(kind=WP), intent(OUT) :: read_baseline, plaquette_baseline
      logical, intent(OUT) :: found

      type(toml_table), allocatable :: root
      type(toml_table), pointer :: perf
      type(toml_table), pointer :: case_tbl
      type(toml_error), allocatable :: err
      integer :: stat
      read_baseline = -1.0_WP
      plaquette_baseline = -1.0_WP
      found = .FALSE.
      call toml_load(root, TRIM(filename), error=err)
      if (ALLOCATED(err)) return
      call get_value(root, "performance", perf, stat=stat)
      if (stat /= 0 .OR. .NOT. ASSOCIATED(perf)) return
      call get_value(perf, "Gen2_8x24n9", case_tbl, stat=stat)
      if (stat /= 0 .OR. .NOT. ASSOCIATED(case_tbl)) return
      call get_value(case_tbl, "read_gauge_time", read_baseline, stat=stat)
      if (stat /= 0) return
      call get_value(case_tbl, "gen_plaquette_time", plaquette_baseline, stat=stat)
      if (stat /= 0) return
      found = (read_baseline > 0.0_WP) .AND. (plaquette_baseline > 0.0_WP)
   end subroutine read_openqcd_baseline

   subroutine read_heatbath_baseline(filename, case_name, update_time_baseline, found)
      character(len=*), intent(IN) :: filename
      character(len=*), intent(IN) :: case_name
      real(kind=WP), intent(OUT) :: update_time_baseline
      logical, intent(OUT) :: found

      type(toml_table), allocatable :: root
      type(toml_table), pointer :: perf
      type(toml_table), pointer :: hb
      type(toml_table), pointer :: case_tbl
      type(toml_error), allocatable :: err
      integer :: stat
      update_time_baseline = -1.0_WP
      found = .FALSE.
      call toml_load(root, TRIM(filename), error=err)
      if (ALLOCATED(err)) return
      call get_value(root, "performance", perf, stat=stat)
      if (stat /= 0 .OR. .NOT. ASSOCIATED(perf)) return
      call get_value(perf, "heatbath", hb, stat=stat)
      if (stat /= 0 .OR. .NOT. ASSOCIATED(hb)) return
      call get_value(hb, TRIM(case_name), case_tbl, stat=stat)
      if (stat /= 0 .OR. .NOT. ASSOCIATED(case_tbl)) return
      call get_value(case_tbl, "update_time_per_traj", update_time_baseline, stat=stat)
      if (stat /= 0) return
      found = (update_time_baseline > 0.0_WP)
   end subroutine read_heatbath_baseline

end module test_performance
