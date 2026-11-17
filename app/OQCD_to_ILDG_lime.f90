Program OQCD_to_ILDG_lime
  use FLUE, only: WP, WC, writeCompiler, writeGit, &
       writeGaugeField_ILDG_lime, readGaugeField_openQCD
  implicit none(type, external)
  character(len=256) :: inputFile, outputFile, NTS, NSS
  complex(kind=WC), dimension(:,:,:,:,:,:,:), allocatable :: U1
  integer :: NT, NS

  call writeCompiler()
  call writeGit()
  write(*,*) ""

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

   CALL writeGaugeField_ILDG_lime(TRIM(outputFile), U1, NS, NS, NS, NT,lfn='lfn', fixSU3=.FALSE.)

 contains

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

 end Program OQCD_to_ILDG_lime
