MODULE test_metadata_wrapper
  USE FLUE_version, ONLY: writeCompiler, writeGit
  USE testdrive, ONLY: new_unittest, unittest_type, error_type, check, new_testsuite, testsuite_type
  IMPLICIT NONE(TYPE, EXTERNAL)
  PRIVATE

  PUBLIC :: collect_metadata_wrapper

CONTAINS

  SUBROUTINE collect_metadata_wrapper(testsuite)
    TYPE(unittest_type), ALLOCATABLE, INTENT(OUT) :: testsuite(:)
    testsuite = [ new_unittest("version_write", test_version_write) ]
  END SUBROUTINE collect_metadata_wrapper

  SUBROUTINE test_version_write(error)
    !! Run the version-writing helpers to ensure they execute without error
    TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
    CALL writeCompiler()
    CALL writeGit()
    CALL check(error, .TRUE., "version write ran without error")

  END SUBROUTINE test_version_write

END MODULE test_metadata_wrapper
