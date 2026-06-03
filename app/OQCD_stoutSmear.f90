!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Read openqcd gaugefield, do some stout smearing
!! Ryan Bignell 2025
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM OQCD_to_ILDG
   USE FLUE, ONLY: WP, WC, writeCompiler, writeGit, &
                   readGaugeField_OpenQCD, &
                   genPlaquette, &
                   StoutSmearLinks
   IMPLICIT NONE(TYPE, EXTERNAL)
   ! IO vars
   CHARACTER(len=256) :: inputFile, NTS, NSS, rhoS, nSmearS
   COMPLEX(kind=WC), DIMENSION(:, :, :, :, :, :, :), ALLOCATABLE :: U1, USmeared
   INTEGER :: nSmear
   REAL(kind=WP) :: rho
   REAL(kind=WP) :: plaq, time, plaqSmear
   INTEGER :: nplaq
   INTEGER :: NT, NS

   CALL writeCompiler()
   CALL writeGit()
   WRITE (*, *) ""

   IF (COMMAND_ARGUMENT_COUNT() > 0) THEN
      CALL GET_COMMAND_ARGUMENT(1, inputFile)
      WRITE (*, *) "Reading gauge file from ", TRIM(inputFile)
      CALL GET_COMMAND_ARGUMENT(2, rhoS)
      WRITE (*, *) "Smearing with rho =  ", TRIM(rhoS)
      READ (rhoS, *) rho
      CALL GET_COMMAND_ARGUMENT(3, nSmearS)
      WRITE (*, *) "Smearing sweeps = ", TRIM(nSmearS)
      CALL str2int(nSmearS, nSmear)
      CALL GET_COMMAND_ARGUMENT(4, NTS)
      CALL str2int(NTS, NT)
      WRITE (*, *) "NT IS", NT
      CALL GET_COMMAND_ARGUMENT(5, NSS)
      CALL str2int(NSS, NS)
      WRITE (*, *) "NS IS", NS
   ELSE
      WRITE (*, *) "Pass the full path to the input toml on the command line"
      WRITE (*, *) "i.e. fpm run OQCD_stoutSmear -- inputFile rho nSmear NT NS"
      STOP
   END IF

   ALLOCATE (U1(3, 3, 4, NT, NS, NS, NS))
   ALLOCATE (USmeared, mold=U1)
   U1 = ReadGaugeField_OpenQCD(TRIM(inputFile), NS, NS, NS, NT, fixSU3=.TRUE.)

   ! First calculate the unsmeared plaquette
   CALL genPlaquette(U1, NT, NS, NS, NS, 1, 4, 4, plaq, nplaq, time)
   ! Now do some smearing
   CALL StoutSmearLinks(U1, rho, nSmear, USmeared)
   ! Now calculate smeared plaquette
   CALL genPlaquette(USmeared, NT, NS, NS, NS, 1, 4, 4, plaqSmear, nplaq, time)
   plaq = plaq / real(nplaq, kind=WP)
   plaqSmear = plaqSmear / real(nplaq, kind=WP)
   WRITE (*, *) 'unsmeared'
   WRITE (*, *) plaq
   WRITE (*, *) 'smeared'
   WRITE (*, *) plaqSmear
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
