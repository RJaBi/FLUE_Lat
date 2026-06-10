!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Convert SU2 HKLS (su2hmc) to CSSM
!! Ryan Bignell 2025
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM SU2_HKLS_to_CSSM
   USE FLUE, ONLY: WP, WC, writeCompiler, writeGit, &
                   ReadGaugeField_HKLS, WriteGaugeField_SU2_CSSM, writeGaugeField_HKLS, SU2_genPlaquette
   USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_DOUBLE_COMPLEX, C_LONG
   IMPLICIT NONE(TYPE, EXTERNAL)
   ! IO vars
   CHARACTER(len=256) :: inputFile, outputFile, NTS, NSS
   COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:, :, :, :, :, :, :), ALLOCATABLE :: U1
   INTEGER :: NT, NS

   INTEGER(kind=C_LONG), DIMENSION(:), ALLOCATABLE :: seed
   INTEGER(kind=C_INT) :: old_nproc

   INTEGER :: ix, iy, iz, it, mu

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
      WRITE (*, *) "i.e. fpm run SU2_HKLS_to_CSSM -- inputFile outputFile NT NS"
      STOP
   END IF

   ALLOCATE (U1(2, 2, 4, NT, NS, NS, NS))
   CALL ReadGaugeField_HKLS(TRIM(inputFile), NS, NS, NS, NT, U1, seed, old_nproc)

   CALL writeGaugeField_SU2_CSSM(TRIM(outputFile), NS, NS, NS, NT, U1, 5, 1.9_WP)

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

END PROGRAM SU2_HKLS_To_CSSM
