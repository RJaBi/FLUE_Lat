!! A module to write the git and compiler versions
!! Is pre-processed by the C pre-processor
module FLUE_version
#ifdef lFORTRAN
   use, intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT
#else
   use, intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT, compiler_version, compiler_options
#endif
   implicit none(external)
   private

   public :: writeCompiler

   public :: writeGit

contains
   subroutine writeCompiler()
      !if (this_image() == 1) then
      write (OUTPUT_UNIT, *) 'This file was compiled by ', &
         COMPILER_VERSION(), ' using the options ', &
         COMPILER_OPTIONS()
      FLUSH (OUTPUT_UNIT)
      !end if
   end subroutine writeCompiler

   subroutine writeGit()
#ifdef SETGITHASH
#include "GITHASH.txt"
      character(len=*), parameter :: git_hash = GITHASH
#else
      character(len=*), parameter :: git_hash = "Unknown"
#endif
      !if (this_image() == 1) then
      write (OUTPUT_UNIT, *) 'Git Commit: ', git_hash
      FLUSH (OUTPUT_UNIT)
      !end if
   end subroutine writeGit
end module FLUE_version
