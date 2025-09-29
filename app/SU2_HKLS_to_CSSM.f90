!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Convert SU2 HKLS (su2hmc) to CSSM
!! Ryan Bignell 2025
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program SU2_HKLS_to_CSSM
  use, intrinsic :: ISO_C_BINDING, only: C_INT, C_DOUBLE_COMPLEX, C_LONG
  use FLUE, only: WP, WC, writeCompiler, writeGit, &
       ReadGaugeField_HKLS, WriteGaugeField_SU2_CSSM
  implicit none(external)
  ! IO vars
  character(len=256) :: inputFile, outputFile, NTS, NSS
  complex(C_DOUBLE_COMPLEX), dimension(:, :, :, :, :, :, :), allocatable :: U1, U2
  integer :: NT, NS

  integer(kind=C_LONG) :: seed
  integer(kind=C_INT) :: old_nproc

  integer(kind=C_LONG) :: seed2
  integer(kind=C_INT) :: old_nproc2

  real(kind=WP) :: sumTrP,time
  integer :: nP


  integer :: ix, iy, iz, it, mu


  call writeCompiler()
  call writeGit()
  write (*, *) ""

  if (COMMAND_ARGUMENT_COUNT() > 0) then
     call GET_COMMAND_ARGUMENT(1, inputFile)
     write (*, *) "Reading gauge file from ", TRIM(inputFile)
     call GET_COMMAND_ARGUMENT(2, outputFile)
     write (*, *) "Saving gauge file to ", TRIM(outputFile)
     call GET_COMMAND_ARGUMENT(3, NTS)
     call str2int(NTS, NT)
     write (*, *) "NT IS", NT
     call GET_COMMAND_ARGUMENT(4, NSS)
     call str2int(NSS, NS)
     write (*, *) "NS IS", NS
  else
     write (*, *) "Pass the full path to the input toml on the command line"
     write (*, *) "i.e. fpm run ILDG_to_OQCD -- inputFile outputFile NT NS"
     stop
  end if

  allocate(U1(NT, NS, NS, NS, 4, 2, 2))
  call ReadGaugeField_HKLS(trim(inputFile), NS, NS, NS, NT, U1, seed, old_nproc)

  ! call writeGaugeField_HKLS(trim(outputFile), NS, NS, NS, NT, U1, seed, old_nproc)

  call writeGaugeField_SU2_CSSM(trim(outputFile), NS, NS, NS, NT, U1, 5, 1.9_WP)


contains

  elemental subroutine str2int(str,int,stat)
    ! Modified from https://stackoverflow.com/a/24077338
    implicit none(external)
    ! Arguments
    character(len=*),intent(in) :: str
    integer,intent(out)         :: int
    integer,optional, intent(out)         :: stat
    character(len=25) :: mystr
    if (present(stat) )then
       read(str,*,iostat=stat)  int
       if (stat /= 0) then
          mystr = 'iostat was not equal to 0'
!          write(*,*) "iostat is ", stat
       end if
    else
       read(str,*)  int
    end if
  end subroutine str2int

end program SU2_HKLS_To_CSSM
