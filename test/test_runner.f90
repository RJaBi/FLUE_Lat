program tester
   use, intrinsic :: ISO_FORTRAN_ENV, only: ERROR_UNIT
   use test_gauge_path, only: collect_gauge_path
   use test_group_generators, only: collect_group_generators
   use test_heatbath_observables, only: collect_heatbath_observables
   use test_matrix_ops, only: collect_matrix_ops
   use test_metadata_wrapper, only: collect_metadata_wrapper
   use test_performance, only: collect_performance
   use test_random_su, only: collect_random_su
   use test_smearing, only: collect_smearing
   use test_updates, only: collect_updates
   use testdrive, only: testsuite_type, run_testsuite, new_testsuite
   implicit none(type, external)

   integer :: stat, is
   type(testsuite_type), allocatable :: testsuites(:)

   stat = 0
   testsuites = [ &
                new_testsuite("group_generators", collect_group_generators), &
                new_testsuite("matrix_ops", collect_matrix_ops), &
                new_testsuite("random_su", collect_random_su), &
                new_testsuite("gauge_path", collect_gauge_path), &
                new_testsuite("smearing", collect_smearing), &
                new_testsuite("updates", collect_updates), &
                new_testsuite("heatbath_observables", collect_heatbath_observables), &
                new_testsuite("metadata_wrapper", collect_metadata_wrapper) &
                !      new_testsuite("performance",      collect_performance)       &
                ]

   do is = 1, SIZE(testsuites)
      write (ERROR_UNIT, '(1x, a)') "Testing: "//testsuites(is)%name
      call run_testsuite(testsuites(is)%collect, ERROR_UNIT, stat, parallel=.false.)
      write (ERROR_UNIT, *) 'These tests finished'
   end do

   if (stat > 0) then
      write (ERROR_UNIT, '(i0, 1x, a)') stat, "test(s) failed!"
      ERROR stop
   else
      write (ERROR_UNIT, *) 'All tests succeeded'
   end if
end program tester
