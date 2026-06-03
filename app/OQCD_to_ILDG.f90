!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Convert openqcd gaugefield to ILDG format
!! Ryan Bignell 2025
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM OQCD_to_ILDG
   USE FLUE, ONLY: WP, WC, writeCompiler, writeGit, &
                   writeGaugeField_ILDG, readGaugeField_OpenQCD, readGaugeField_ILDG
   IMPLICIT NONE(TYPE, EXTERNAL)
   ! IO vars
   CHARACTER(len=256) :: inputFile, outputFile, NTS, NSS
   COMPLEX(kind=WC), DIMENSION(:, :, :, :, :, :, :), ALLOCATABLE :: U1
!CHECKS!  complex(kind=WC), dimension(:, :, :, :, :, :, :), allocatable :: U2
   INTEGER :: NT, NS

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
   U1 = ReadGaugeField_OpenQCD(TRIM(inputFile), NS, NS, NS, NT, fixSU3=.TRUE.)

   CALL writeGaugeField_ILDG(TRIM(outputFile), U1, NS, NS, NS, NT, fixSU3=.FALSE.)

!CHECKS!
!CHECKS!  allocate(U2(3, 3, 4, NT, NS, NS, NS))
!CHECKS!  U2 = ReadGaugeField_ILDG(trim(outputFile), NS, NS, NS, NT, fixSU3=.false.)
!CHECKS!
!CHECKS!    if (any(U1 /= U2)) then
!CHECKS!     write(*,*) 'bad'
!CHECKS!  else
!CHECKS!     write(*,*) 'good'
!CHECKS!  end if

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

END PROGRAM OQCD_To_ILDG
