!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Convert SU2 HKLS (su2hmc) to CSSM
!! Ryan Bignell 2025
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM SU2_HKLS_to_CSSM
   USE FLUE, ONLY: WP, WC, writeCompiler, writeGit, &
                   ReadGaugeField_HKLS, writeGaugeField_NRQ2CD, SU2_genPlaquette
   USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_DOUBLE_COMPLEX, C_LONG
   IMPLICIT NONE(TYPE, EXTERNAL)
   ! IO vars
   CHARACTER(len=256) :: inputFile, outputFile, NTS, NSS, endianess
   COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(:, :, :, :, :, :, :), ALLOCATABLE :: U1
   INTEGER :: NT, NS
   LOGICAL :: bigEndian

   INTEGER(kind=C_LONG), DIMENSION(:), ALLOCATABLE :: seed
   INTEGER(kind=C_INT) :: old_nproc

   REAL(kind=WP) :: sumTrp, time
   INTEGER :: nPlaq
   REAL(kind=WP) :: aPlaq
   REAL(kind=WP) :: u0

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
      IF (COMMAND_ARGUMENT_COUNT() == 5) THEN
         CALL GET_COMMAND_ARGUMENT(5, endianess)
         IF (trim(endianess) == 'big') THEN
            bigEndian = .TRUE.
         ELSE
            bigEndian = .FALSE.
         END IF
      ELSE
         bigEndian = .FALSE.
      END IF
   ELSE
      WRITE (*, *) "Pass the full path to the input toml on the command line"
      WRITE (*, *) "i.e. fpm run SU2_HKLS_to_NRQ2CD -- inputFile outputFile NT NS [endianess=big|little (default)]"
      STOP
   END IF

   ALLOCATE (U1(2, 2, 4, NT, NS, NS, NS))
   CALL ReadGaugeField_HKLS(TRIM(inputFile), NS, NS, NS, NT, U1, seed, old_nproc, bigEndian=bigEndian)

   CALL writeGaugeField_NRQ2CD(TRIM(outputFile), NS, NS, NS, NT, U1)
   ! Calculate plaquette and mean link
   ! total
   CALL SU2_genPlaquette(U1, NT, NS, NS, NS, 1, 4, 4, sumtrp, nplaq, time)
   aplaq = sumtrp / real(nplaq, kind=WP)
   u0 = aplaq**0.25_WP
   WRITE(*,*) 'Plaquette', aplaq
   WRITE(*,*) 'u0', u0
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
