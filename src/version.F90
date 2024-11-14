!! A module to write the git and compiler versions
!! Is pre-processed by the C pre-processor
module FLUE_version
   use, intrinsic :: iso_fortran_env, only: output_unit, compiler_version, compiler_options
   implicit none
   private

   public :: writeCompiler

   public :: writeGit

contains
   subroutine writeCompiler()
      if (this_image() == 1) then
         write (output_unit, *) 'This file was compiled by ', &
            compiler_version(), ' using the options ', &
            compiler_options()
         flush (output_unit)
      end if
   end subroutine writeCompiler

   subroutine writeGit()
#ifdef SETGITHASH
#include "GITHASH.txt"
      character(len=*), parameter :: git_hash = GITHASH
#else
      character(len=*), parameter :: git_hash = "Unknown"
#endif
      if (this_image() == 1) then
         write (output_unit, *) 'Git Commit: ', git_hash
         flush (output_unit)
      end if
   end subroutine writeGit
end module FLUE_version
