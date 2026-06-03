!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Convert ildg-bin big endian gaugefield to openqcd format
!! Ryan Bignell 2025
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM ILDG_To_OQCD
   USE FLUE, ONLY: WP, WC, writeCompiler, writeGit, &
                   ReadGaugeField_ILDG, writeGaugeField_OpenQCD  !CHECKS! ,&
!Checks!       ReadGaugeField_OpenQCD, genPlaquette
   IMPLICIT NONE(TYPE, EXTERNAL)
   ! IO vars
   CHARACTER(len=256) :: inputFile, outputFile, NTS, NSS
   COMPLEX(kind=WC), DIMENSION(:, :, :, :, :, :, :), ALLOCATABLE :: U1
!Checks!  complex(kind=WC), dimension(:, :, :, :, :, :, :), allocatable :: U2
   INTEGER :: NT, NS
!Checks!  ! Plaqeutte holders
!Checks!  real(kind=WP) :: sumTrP, time
   !Checks!  integer :: nP

   CALL writeCompiler()
   CALL writeGit()
   WRITE (*, *) ""

   IF (COMMAND_ARGUMENT_COUNT() > 0) THEN
      CALL GET_COMMAND_ARGUMENT(1, inputFile)
      WRITE (*, *) "Reading gauge file from ", TRIM(inputFile)
      CALL GET_COMMAND_ARGUMENT(2, outputFile)
      WRITE (*, *) "Saving gauge file to ", TRIM(outputFile)
      CALL GET_COMMAND_ARGUMENT(3, NTS)
      CALL str2int(NTS, NT)
      WRITE (*, *) "NT IS", NT
      CALL GET_COMMAND_ARGUMENT(4, NSS)
      CALL str2int(NSS, NS)
      WRITE (*, *) "NS IS", NS
   ELSE
      WRITE (*, *) "Pass the full path to the input toml on the command line"
      WRITE (*, *) "i.e. fpm run ILDG_to_OQCD -- inputFile outputFile NT NS"
      STOP
   END IF

   ALLOCATE (U1(3, 3, 4, NT, NS, NS, NS))
   U1 = ReadGaugeField_ILDG(TRIM(inputFile), NS, NS, NS, NT, fixSU3=.FALSE.)

!Checks!  call genPlaquette(U1, NT, NS, NS, NS, 1, 4, 4, sumTrP, nP, time)
!Checks!  write(*,*) sumTrp / real(nP, kind=WP)

   CALL writeGaugeField_OpenQCD(TRIM(outputFile), U1, NS, NS, NS, NT)

!Checks!  allocate(U2(3, 3, 4, NT, NS, NS, NS))
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

CONTAINS

   ELEMENTAL SUBROUTINE str2int(str, int, stat)
      ! Modified from https://stackoverflow.com/a/24077338
      IMPLICIT NONE(EXTERNAL)
      ! Arguments
      CHARACTER(len=*), INTENT(IN) :: str
      INTEGER, INTENT(OUT) :: int
      INTEGER, OPTIONAL, INTENT(OUT) :: stat
      CHARACTER(len=25) :: mystr
      IF (PRESENT(stat)) THEN
         READ (str, *, iostat=stat) int
         IF (stat /= 0) THEN
            mystr = 'iostat was not equal to 0'
!          write(*,*) "iostat is ", stat
         END IF
      ELSE
         READ (str, *) int
      END IF
   END SUBROUTINE str2int

END PROGRAM ILDG_To_OQCD
