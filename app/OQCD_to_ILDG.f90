!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Convert openqcd gaugefield to ILDG format
!! Ryan Bignell 2025
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program OQCD_to_ILDG
  use FLUE, only: WP, WC, writeCompiler, writeGit, &
       writeGaugeField_ILDG, readGaugeField_OpenQCD, readGaugeField_ILDG
  implicit none(external)
  ! IO vars
  character(len=256) :: inputFile, outputFile, NTS, NSS
  complex(kind=WC), dimension(:, :, :, :, :, :, :), allocatable :: U1
!CHECKS!  complex(kind=WC), dimension(:, :, :, :, :, :, :), allocatable :: U2
  integer :: NT, NS

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

  allocate(U1(NT, NS, NS, NS, 4, 3, 3))
  U1 = ReadGaugeField_OpenQCD(trim(inputFile), NS, NS, NS, NT, fixSU3=.true.)

  call writeGaugeField_ILDG(trim(outputFile), U1, NS, NS, NS, NT, fixSU3=.false.)

!CHECKS!
!CHECKS!  allocate(U2(NT, NS, NS, NS, 4, 3, 3))
!CHECKS!  U2 = ReadGaugeField_ILDG(trim(outputFile), NS, NS, NS, NT, fixSU3=.false.)
!CHECKS!
!CHECKS!    if (any(U1 /= U2)) then
!CHECKS!     write(*,*) 'bad'
!CHECKS!  else
!CHECKS!     write(*,*) 'good'
!CHECKS!  end if

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

end program OQCD_To_ILDG
