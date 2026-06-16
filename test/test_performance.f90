MODULE test_performance
   USE FLUE_constants, ONLY : WP, WC, SP
   USE FLUE_heatbath, ONLY : updateLinks
   USE FLUE_matrixConstants, ONLY : Ident3x3
   USE FLUE_openQCDFileIO_SA, ONLY : ReadGaugeField_OpenQCD
   USE FLUE_wloops, ONLY : genPlaquette
   USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
   USE M_stopwatch, ONLY : watchtype, create_watch, start_watch, stop_watch, &
                           destroy_watch, read_watch
   USE test_helpers, ONLY: seed_rng_fixed, fill_identity_su3
   USE testdrive, ONLY : new_unittest, unittest_type, error_type, check
   USE tomlf, ONLY : toml_table, toml_error, toml_load, get_value
   IMPLICIT NONE(TYPE, EXTERNAL)
   PRIVATE
   PUBLIC :: collect_performance

CONTAINS

   SUBROUTINE collect_performance(testsuite)
      TYPE(unittest_type), ALLOCATABLE, INTENT(OUT) :: testsuite(:)

      testsuite = [ &
         new_unittest("openqcd_performance",    test_openqcd_performance), &
         new_unittest("heatbath_wilson_xi1",    test_heatbath_wilson_xi1), &
         new_unittest("heatbath_wilson_xi10",   test_heatbath_wilson_xi10), &
         new_unittest("heatbath_symanzik_xi1",  test_heatbath_symanzik_xi1), &
         new_unittest("heatbath_symanzik_xi10", test_heatbath_symanzik_xi10) &
      ]
   END SUBROUTINE collect_performance

   !=========================================================
   ! Existing OpenQCD read/plaquette benchmark
   !=========================================================
   SUBROUTINE test_openqcd_performance(error)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      COMPLEX(kind=WC), ALLOCATABLE :: U(:, :, :, :, :, :, :)
      TYPE(watchtype) :: watch
      REAL(kind=WP) :: read_time, plaquette_time, baseline_read, baseline_plaq, tolerance
      REAL(kind=WP) :: sumtrp
      REAL(kind=SP) :: watchtime
      CHARACTER(len=128) :: gauge_file, baseline_path
      INTEGER :: NT, NS, np
      LOGICAL :: found_baseline
      NT = 8
      NS = 24
      gauge_file = 'testdata/Gen2_8x24n9'
      baseline_path = 'testdata/reference_values.toml'
      tolerance = 0.30_WP
      CALL read_openqcd_baseline(baseline_path, baseline_read, baseline_plaq, found_baseline)
      CALL check(error, found_baseline, &
           'performance baseline section [performance.Gen2_8x24n9] must be present in '//trim(baseline_path))
      IF (allocated(error)) RETURN
      CALL create_watch(watch)
      CALL start_watch(watch)
      U = ReadGaugeField_OpenQCD(trim(gauge_file), NS, NS, NS, NT)
      CALL stop_watch(watch)
      CALL read_watch(watchtime, watch, 'wall')
      CALL destroy_watch(watch)
      read_time = real(watchtime, kind=WP)
      CALL genPlaquette(U, NT, NS, NS, NS, 1, 4, 4, sumtrp, np, plaquette_time)
      WRITE(*,'(A,F10.6)') 'Measured ReadGaugeField_OpenQCD time: ', read_time
      WRITE(*,'(A,F10.6)') 'Measured genPlaquette time:          ', plaquette_time
      CALL check(error, read_time <= baseline_read * (1.0_WP + tolerance), &
           'ReadGaugeField_OpenQCD time within baseline tolerance')
      IF (allocated(error)) RETURN
      CALL check(error, plaquette_time <= baseline_plaq * (1.0_WP + tolerance), &
           'genPlaquette time within baseline tolerance')
   END SUBROUTINE test_openqcd_performance
   !=========================================================
   ! Heatbath benchmark cases
   !=========================================================
   SUBROUTINE test_heatbath_wilson_xi1(error)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      CALL benchmark_heatbath_case(error, "wilson_xi1", "wilson", 1.0_WP, 6.0_WP)
   END SUBROUTINE test_heatbath_wilson_xi1

   SUBROUTINE test_heatbath_wilson_xi10(error)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      CALL benchmark_heatbath_case(error, "wilson_xi10", "wilson", 10.0_WP, 6.8_WP)
   END SUBROUTINE test_heatbath_wilson_xi10

   SUBROUTINE test_heatbath_symanzik_xi1(error)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      CALL benchmark_heatbath_case(error, "symanzik_xi1", "symanzik", 1.0_WP, 6.0_WP)
   END SUBROUTINE test_heatbath_symanzik_xi1

   SUBROUTINE test_heatbath_symanzik_xi10(error)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      CALL benchmark_heatbath_case(error, "symanzik_xi10", "symanzik", 10.0_WP, 6.8_WP)
   END SUBROUTINE test_heatbath_symanzik_xi10

   SUBROUTINE benchmark_heatbath_case(error, case_name, action_tag, xi, beta)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      CHARACTER(len=*), INTENT(IN) :: case_name, action_tag
      REAL(WP), INTENT(IN) :: xi, beta
      INTEGER, PARAMETER :: NS = 4
      INTEGER, PARAMETER :: NT = 8
      INTEGER, PARAMETER :: nTraj = 8
      INTEGER, PARAMETER :: nRepeat = 3
      REAL(WP), PARAMETER :: tolerance = 0.30_WP
      COMPLEX(WC) :: U(3,3,4,NT,NS,NS,NS)
      COMPLEX(WC) :: UNew(3,3,4,NT,NS,NS,NS)
      REAL(WP) :: total_time(nRepeat), time_per_traj(nRepeat), median_t
      REAL(WP) :: baseline_update
      LOGICAL  :: found_baseline
      INTEGER :: irep, iTraj
      REAL(WP) :: aplaq, splaq, tplaq
      TYPE(watchtype) :: watch
      REAL(kind=SP) :: watchtime
      REAL(WP) :: sumTrP, plaq_time
      INTEGER :: nPlaq
      DO irep = 1, nRepeat
         CALL fill_identity_su3(U)
         UNew = U
         ! Fixed stdlib_random seed for repeatable benchmark trajectories
         CALL seed_rng_fixed(20240615 + irep)
         ! Warm-up call to reduce first-call effects in the timed region
         CALL updateLinks(U, beta, xi, UNew, trim(action_tag))
         U = UNew
         CALL create_watch(watch)
         CALL start_watch(watch)
         DO iTraj = 1, nTraj
            CALL updateLinks(U, beta, xi, UNew, trim(action_tag))
            U = UNew
         END DO
         CALL stop_watch(watch)
         CALL read_watch(watchtime, watch, 'wall')
         CALL destroy_watch(watch)
         total_time(irep) = real(watchtime, kind=WP)
         time_per_traj(irep) = total_time(irep) / real(nTraj, kind=WP)
      END DO
      median_t = median3(time_per_traj)
      ! Final plaquette sanity values from the last run
      CALL genPlaquette(U, NT, NS, NS, NS, 1, 4, 4, sumTrP, nPlaq, plaq_time)
      aplaq = sumTrP / (3.0_WP * real(nPlaq, kind=WP))

      IF (xi /= 1.0_WP) THEN
         CALL genPlaquette(U, NT, NS, NS, NS, 2, 4, 4, sumTrP, nPlaq, plaq_time)
         splaq = sumTrP / (3.0_WP * real(nPlaq, kind=WP))

         CALL genPlaquette(U, NT, NS, NS, NS, 1, 1, 4, sumTrP, nPlaq, plaq_time)
         tplaq = sumTrP / (3.0_WP * real(nPlaq, kind=WP))
      ELSE
         splaq = 0.0_WP
         tplaq = 0.0_WP
      END IF

      CALL read_heatbath_baseline("testdata/reference_values.toml", case_name, baseline_update, found_baseline)
      CALL check(error, found_baseline, &
           "heatbath performance baseline [performance.heatbath."//trim(case_name) &
           //"] must exist in testdata/reference_values.toml")
      IF (allocated(error)) RETURN

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

      CALL check(error, median_t <= baseline_update * (1.0_WP + tolerance), &
           "heatbath update time per trajectory within baseline tolerance")
      IF (allocated(error)) RETURN

      CALL check(error, ieee_is_finite(aplaq) .and. abs(aplaq) <= 1.0_WP + 1.0e-12_WP, &
           "average plaquette should be finite and bounded")
      IF (allocated(error)) RETURN

      IF (xi /= 1.0_WP) THEN
         CALL check(error, ieee_is_finite(splaq) .and. abs(splaq) <= 1.0_WP + 1.0e-12_WP, &
              "spatial plaquette should be finite and bounded")
         IF (allocated(error)) RETURN

         CALL check(error, ieee_is_finite(tplaq) .and. abs(tplaq) <= 1.0_WP + 1.0e-12_WP, &
              "temporal plaquette should be finite and bounded")
      END IF
   END SUBROUTINE benchmark_heatbath_case


   !=========================================================
   ! Helpers
   !=========================================================
   PURE REAL(WP) FUNCTION median3(x) RESULT(m)
      REAL(WP), INTENT(IN) :: x(3)
      REAL(WP) :: a, b, c
      a = x(1)
      b = x(2)
      c = x(3)
      IF ((a <= b .and. b <= c) .or. (c <= b .and. b <= a)) THEN
         m = b
      ELSE IF ((b <= a .and. a <= c) .or. (c <= a .and. a <= b)) THEN
         m = a
      ELSE
         m = c
      END IF
   END FUNCTION median3

   SUBROUTINE read_openqcd_baseline(filename, read_baseline, plaquette_baseline, found)
      CHARACTER(len=*), INTENT(IN)  :: filename
      REAL(kind=WP),    INTENT(OUT) :: read_baseline, plaquette_baseline
      LOGICAL,          INTENT(OUT) :: found

      TYPE(toml_table), ALLOCATABLE :: root
      TYPE(toml_table), POINTER     :: perf
      TYPE(toml_table), POINTER     :: case_tbl
      TYPE(toml_error), ALLOCATABLE :: err
      INTEGER :: stat
      read_baseline      = -1.0_WP
      plaquette_baseline = -1.0_WP
      found              = .FALSE.
      CALL toml_load(root, trim(filename), error=err)
      IF (allocated(err)) RETURN
      CALL get_value(root, "performance", perf, stat=stat)
      IF (stat /= 0 .or. .not. associated(perf)) RETURN
      CALL get_value(perf, "Gen2_8x24n9", case_tbl, stat=stat)
      IF (stat /= 0 .or. .not. associated(case_tbl)) RETURN
      CALL get_value(case_tbl, "read_gauge_time", read_baseline, stat=stat)
      IF (stat /= 0) RETURN
      CALL get_value(case_tbl, "gen_plaquette_time", plaquette_baseline, stat=stat)
      IF (stat /= 0) RETURN
      found = (read_baseline > 0.0_WP) .and. (plaquette_baseline > 0.0_WP)
   END SUBROUTINE read_openqcd_baseline


   SUBROUTINE read_heatbath_baseline(filename, case_name, update_time_baseline, found)
      CHARACTER(len=*), INTENT(IN)  :: filename
      CHARACTER(len=*), INTENT(IN)  :: case_name
      REAL(kind=WP),    INTENT(OUT) :: update_time_baseline
      LOGICAL,          INTENT(OUT) :: found

      TYPE(toml_table), ALLOCATABLE :: root
      TYPE(toml_table), POINTER     :: perf
      TYPE(toml_table), POINTER     :: hb
      TYPE(toml_table), POINTER     :: case_tbl
      TYPE(toml_error), ALLOCATABLE :: err
      INTEGER :: stat
      update_time_baseline = -1.0_WP
      found = .FALSE.
      CALL toml_load(root, trim(filename), error=err)
      IF (allocated(err)) RETURN
      CALL get_value(root, "performance", perf, stat=stat)
      IF (stat /= 0 .or. .not. associated(perf)) RETURN
      CALL get_value(perf, "heatbath", hb, stat=stat)
      IF (stat /= 0 .or. .not. associated(hb)) RETURN
      CALL get_value(hb, trim(case_name), case_tbl, stat=stat)
      IF (stat /= 0 .or. .not. associated(case_tbl)) RETURN
      CALL get_value(case_tbl, "update_time_per_traj", update_time_baseline, stat=stat)
      IF (stat /= 0) RETURN
      found = (update_time_baseline > 0.0_WP)
   END SUBROUTINE read_heatbath_baseline

END MODULE test_performance
