module FLUE_version
  !< Module FLUE_version
  !< Provides subroutines to write the compiler & compiler options used
  !< as well as the git Hash (via preprocessor)
   use, intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT, compiler_version, compiler_options
   implicit none(type, external)
   private

   public :: writeCompiler

   public :: writeGit

contains
  subroutine writeCompiler()
    !< Writes the compiler_version and compiler_options to output_unit
    write (OUTPUT_UNIT, *) 'This file was compiled by ', &
         COMPILER_VERSION(), ' using the options ', &
         COMPILER_OPTIONS()
    FLUSH (OUTPUT_UNIT)
   end subroutine writeCompiler

   subroutine writeGit()
     !< Via the SETGITHASH macro
     !< and the 'GITHASH.txt' file
     !< writes the contents of that file to output_unit
#ifdef SETGITHASH
#include "GITHASH.txt"
     character(len=*), parameter :: git_hash = GITHASH
     !< The included git hash.
#else
     character(len=*), parameter :: git_hash = "Unknown"
     !< The included git hash
#endif
      write (OUTPUT_UNIT, *) 'Git Commit: ', git_hash
      FLUSH (OUTPUT_UNIT)
   end subroutine writeGit
end module FLUE_version
