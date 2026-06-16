PROGRAM tester
   USE, INTRINSIC :: iso_fortran_env, ONLY : error_unit
   USE test_gauge_path,       ONLY: collect_gauge_path
   USE test_group_generators, ONLY: collect_group_generators
   USE test_heatbath_observables,         ONLY: collect_heatbath_observables
   USE test_matrix_ops,       ONLY: collect_matrix_ops
   USE test_metadata_wrapper, ONLY: collect_metadata_wrapper
   USE test_performance,      ONLY: collect_performance
   USE test_random_su,        ONLY: collect_random_su
   USE test_smearing,         ONLY: collect_smearing
   USE test_updates,          ONLY: collect_updates
   USE testdrive,             ONLY: testsuite_type, run_testsuite, new_testsuite
   IMPLICIT NONE(TYPE, EXTERNAL)

   INTEGER :: stat, is
   TYPE(testsuite_type), ALLOCATABLE :: testsuites(:)

   stat = 0
   testsuites = [ &
      new_testsuite("group_generators", collect_group_generators), &
      new_testsuite("matrix_ops",       collect_matrix_ops),       &
      new_testsuite("random_su",        collect_random_su),        &
      new_testsuite("gauge_path",       collect_gauge_path),       &
      new_testsuite("smearing",         collect_smearing),         &
      new_testsuite("updates",          collect_updates),          &
      new_testsuite("heatbath_observables",  collect_heatbath_observables), &
      new_testsuite("metadata_wrapper", collect_metadata_wrapper), &
      new_testsuite("performance",      collect_performance)       &
   ]

   DO is = 1, size(testsuites)
      WRITE(error_unit, '(1x, a)') "Testing: " // testsuites(is)%name
      CALL run_testsuite(testsuites(is)%collect, error_unit, stat)
   END DO

   IF (stat > 0) THEN
      WRITE(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      ERROR STOP
   ELSE
      WRITE(error_unit, *) 'All tests succeeded'
   END IF
END PROGRAM tester
