!! Functions to read and write ILDG binary data formats as from cola
MODULE FLUE_ILDG_bin
   USE FLUE_constants, ONLY: WP, WC
   USE FLUE_SU3MatrixOps, ONLY: FixSU3Matrix
   !use stdlib_linalg, only: det
   USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: OUTPUT_UNIT
   IMPLICIT NONE(TYPE, EXTERNAL)
   PRIVATE
   PUBLIC :: ReadGaugeField_ILDG
   PUBLIC :: writeGaugeField_ILDG

CONTAINS

   FUNCTION ReadGaugeField_ILDG(filename, NX, NY, NZ, NT, fixSU3) RESULT(U_xd)
      CHARACTER(len=*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: NX, NY, NZ, NT
      LOGICAL, OPTIONAL, INTENT(IN) :: fixSU3
      COMPLEX(kind=WC), DIMENSION(3, 3, 4, NT, NX, NY, NZ) :: U_xd
      COMPLEX(kind=WC), DIMENSION(3, 3, 4, NX, NY, NZ, NT) :: URead
      INTEGER, PARAMETER :: infl = 101
      INTEGER :: matrix_len, irecl
      ! counters
      INTEGER :: it, ix, iy, iz, mu, nu
      LOGICAL :: fixSU3Set

      IF (PRESENT(fixSU3)) THEN
         fixSU3Set = fixSU3
      ELSE
         fixSU3Set = .TRUE.
      END IF
      ! First read the gaugefield
      WRITE (OUTPUT_UNIT, *) TRIM(filename)
      matrix_len = 16 * 3 * 3
      irecl = matrix_len * 4 * NX * NY * NZ
      ! write(*,*) matrix_len, irecl, 3, 3, 4, nx, ny, nz, nt
      OPEN (infl, file=TRIM(filename), form='unformatted', access='direct', &
            status='old', action='read', recl=irecl, convert='BIG_ENDIAN')
      DO it = 1, NT
         READ (infl, rec=it) URead(:, :, :, :, :, :, it)
      END DO
      CLOSE (infl)
      ! THen really dumbly re-order the indices
      DO it = 1, NT
         DO ix = 1, NX
            DO iy = 1, NY
               DO iz = 1, NZ
                  DO mu = 1, 4
                     ! Mu is the openqcd index
                     ! i.e. starts t,x,y,z
                     ! nu is the ILDG index
                     ! i.e. starts x,y,z,t
                     SELECT CASE (mu)
                     CASE (1)
                        nu = 2
                     CASE (2)
                        nu = 3
                     CASE (3)
                        nu = 4
                     CASE (4)
                        nu = 1
                     END SELECT
                     U_xd(:, :, nu, it, ix, iy, iz) = TRANSPOSE(URead(:, :, mu, ix, iy, iz, it))
                     IF (fixSU3Set) THEN
                        CALL FixSU3Matrix(U_xd(:, :, nu, it, ix, iy, iz))
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO

   END FUNCTION ReadGaugeField_ILDG

   SUBROUTINE writeGaugeField_ILDG(filename, U_xd, NX, NY, NZ, NT, fixSU3)
      CHARACTER(len=*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: NX, NY, NZ, NT
      LOGICAL, OPTIONAL, INTENT(IN) :: fixSU3
      COMPLEX(kind=WC), DIMENSION(3, 3, 4, NT, NX, NY, NZ), INTENT(IN) :: U_xd
      COMPLEX(kind=WC), DIMENSION(3, 3, 4, NX, NY, NZ, NT) :: URead
      INTEGER, PARAMETER :: infl = 101
      INTEGER :: matrix_len, irecl
      ! counters
      INTEGER :: it, ix, iy, iz, mu, nu
      LOGICAL :: fixSU3Set

      IF (PRESENT(fixSU3)) THEN
         fixSU3Set = fixSU3
      ELSE
         fixSU3Set = .TRUE.
      END IF

      ! First really dumbly re-order the indices
      DO it = 1, NT
         DO ix = 1, NX
            DO iy = 1, NY
               DO iz = 1, NZ
                  DO mu = 1, 4
                     ! Mu is the openqcd index
                     ! i.e. starts t,x,y,z
                     ! nu is the ILDG index
                     ! i.e. starts x,y,z,t
                     SELECT CASE (mu)
                     CASE (1)
                        nu = 2
                     CASE (2)
                        nu = 3
                     CASE (3)
                        nu = 4
                     CASE (4)
                        nu = 1
                     END SELECT
                     URead(:, :, mu, ix, iy, iz, it) = TRANSPOSE(U_xd(:, :, nu, it, ix, iy, iz))
                     IF (fixSU3Set) THEN
                        CALL FixSU3Matrix(URead(:, :, mu, ix, iy, iz, it))
                     END IF

                  END DO
               END DO
            END DO
         END DO
      END DO

      ! Then write the gaugefield
      matrix_len = 16 * 3 * 3
      irecl = matrix_len * 4 * NX * NY * NZ
      ! write(*,*) matrix_len, irecl, 3, 3, 4, nx, ny, nz, nt
      OPEN (infl, file=TRIM(filename), form='unformatted', access='direct', &
            status='replace', action='write', recl=irecl, convert='BIG_ENDIAN')
      DO it = 1, NT
         WRITE (infl, rec=it) URead(:, :, :, :, :, :, it)
      END DO
      CLOSE (infl)

   END SUBROUTINE WriteGaugeField_ILDG

END MODULE FLUE_ILDG_bin
