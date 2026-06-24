module test_metadata_wrapper
   use FLUE_version, only: writeCompiler, writeGit
   use testdrive, only: new_unittest, unittest_type, error_type, check, new_testsuite, testsuite_type
   implicit none(type, external)
   private

   public :: collect_metadata_wrapper

contains

   subroutine collect_metadata_wrapper(testsuite)
      type(unittest_type), allocatable, intent(OUT) :: testsuite(:)
      testsuite = [new_unittest("version_write", test_version_write)]
   end subroutine collect_metadata_wrapper

   subroutine test_version_write(error)
    !! Run the version-writing helpers to ensure they execute without error
      type(error_type), allocatable, intent(OUT) :: error
      call writeCompiler()
      call writeGit()
      call check(error, .TRUE., "version write ran without error")

   end subroutine test_version_write

end module test_metadata_wrapper
