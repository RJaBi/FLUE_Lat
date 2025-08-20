!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Convert ildg-bin big endian gaugefield to openqcd format
!! Ryan Bignell 2025
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program ILDG_To_OQCD
  use FLUE, only: WP, WC, writeCompiler, writeGit, &
       ReadGaugeField_ILDG, writeGaugeField_OpenQCD  !CHECKS! ,&
!Checks!       ReadGaugeField_OpenQCD, genPlaquette
  implicit none(external)
  ! IO vars
  character(len=256) :: inputFile, outputFile, NTS, NSS
  complex(kind=WC), dimension(:, :, :, :, :, :, :), allocatable :: U1
!Checks!  complex(kind=WC), dimension(:, :, :, :, :, :, :), allocatable :: U2
  integer :: NT, NS
!Checks!  ! Plaqeutte holders
!Checks!  real(kind=WP) :: sumTrP, time
  !Checks!  integer :: nP

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
  U1 = ReadGaugeField_ILDG(trim(inputFile), NS, NS, NS, NT, fixSU3=.false.)

!Checks!  call genPlaquette(U1, NT, NS, NS, NS, 1, 4, 4, sumTrP, nP, time)
!Checks!  write(*,*) sumTrp / real(nP, kind=WP)

  call writeGaugeField_OpenQCD(trim(outputFile), U1, NS, NS, NS, NT)


!Checks!  allocate(U2(NT, NS, NS, NS, 4, 3, 3))
!Checks!  U2 = ReadGaugeField_OpenQCD(trim(outputFile), NS, NS, NS, NT, fixSU3=.false.)
!Checks!
!Checks!  call genPlaquette(U2, NT, NS, NS, NS, 1, 4, 4, sumTrP, nP, time)
!Checks!  write(*,*) sumTrp / real(nP, kind=WP)
!Checks!
!Checks!  if (any(U1 /= U2)) then
!Checks!     write(*,*) 'bad'
!Checks!  else
!Checks!     write(*,*) 'good'
!Checks!  end if

contains

  elemental subroutine str2int(str,int,stat)
    ! Modified from https://stackoverflow.com/a/24077338
    implicit none(external)
    ! Arguments
    character(len=*),intent(in) :: str
    integer,intent(out)         :: int
    integer,optional, intent(out)         :: stat
    if (present(stat) )then
       read(str,*,iostat=stat)  int
       if (stat /= 0) then
          str = 'iostat was not equal to 0'
!          write(*,*) "iostat is ", stat
       end if
    else
       read(str,*)  int
    end if
  end subroutine str2int

end program ILDG_To_OQCD
