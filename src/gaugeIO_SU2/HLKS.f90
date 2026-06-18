!! Functions to read & write SU2 HKLS format
MODULE FLUE_SU2_HKLS
   USE FLUE_constants, ONLY: WP, WC
   USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_DOUBLE_COMPLEX, C_LONG
   !use FLUE_SU3MatrixOps, only: FixSU3Matrix
   !use FLUE_wloops, only: genPlaquette
   IMPLICIT NONE(TYPE, EXTERNAL)
   PRIVATE
   PUBLIC :: ReadGaugeField_HKLS
   PUBLIC :: writeGaugeField_HKLS

CONTAINS

   SUBROUTINE ReadGaugeField_HKLS(filename, NX, NY, NZ, NT, U_xd, seed, old_nproc, bigEndian)
      CHARACTER(len=*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: NX, NY, NZ, NT
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(2, 2, 4, NT, NX, NY, NZ), INTENT(OUT) :: U_xd
      INTEGER(kind=C_INT), INTENT(OUT) :: old_nproc
      INTEGER(kind=C_LONG), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: seed
      LOGICAL, INTENT(IN), OPTIONAL :: bigEndian
      ! for reading
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(NX, NY, NZ, NT, 4, 2) :: UHold
      LOGICAL :: littleEndian
      !complex(kind=C_DOUBLE_COMPLEX), dimension(4) :: UTmp
      COMPLEX(kind=C_DOUBLE_COMPLEX) :: UTmp
      INTEGER, PARAMETER :: infl = 107
      ! counters
      INTEGER :: it, ix, iy, iz, mu, ab, aa, bb
      INTEGER :: ioStatus

      IF (present(bigEndian)) THEN
         littleEndian = .not. bigEndian
      ELSE
         littleEndian = .TRUE.
      END IF
      IF (littleEndian) THEN
         OPEN (infl, file=TRIM(filename), form="unformatted", access="stream", status="old", action="read", convert="little_endian")
      ELSE
         OPEN (infl, file=TRIM(filename), form="unformatted", access="stream", status="old", action="read", convert="big_endian")
      END IF
      READ (infl) old_nproc
      DO ab = 1, 2
         DO mu = 1, 4
            DO it = 1, NT
               DO iz = 1, NZ
                  DO iy = 1, NY
                     DO ix = 1, NX
                        READ (infl) UTmp
                        UHold(ix, iy, iz, it, mu, ab) = UTmp
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO

      ALLOCATE(seed(old_nproc))
      READ(infl, iostat=IOStatus) seed
      IF (IOStatus /= 0) THEN
         ! Fortran produced
         ! rewind
         REWIND(infl)
         ! make it length 1 instead
         DEALLOCATE(seed)
         ALLOCATE(seed(1))
         READ(infl) seed
      END IF
      CLOSE (infl)
      ! put it in 2x2
      DO it = 1, NT
         DO iz = 1, NZ
            DO iy = 1, NY
               DO ix = 1, NX
                  ! Unpack
                  U_xd(1, 1, :, it, ix, iy, iz) = UHold(ix, iy, iz, it, :, 1)
                  U_xd(2, 2, :, it, ix, iy, iz) = CONJG(UHold(ix, iy, iz, it, :, 1))
                  U_xd(1, 2, :, it, ix, iy, iz) = UHold(ix, iy, iz, it, :, 2)
                  U_xd(2, 1, :, it, ix, iy, iz) = -CONJG(UHold(ix, iy, iz, it, :, 2))

               END DO
            END DO
         END DO
      END DO
      ! Shift the mu ordering
      U_xd = CSHIFT(U_xd, -1, dim=3)
   END SUBROUTINE ReadGaugeField_HKLS

   SUBROUTINE WriteGaugeField_HKLS(filename, NX, NY, NZ, NT, U_xd, seed, old_nproc)
      CHARACTER(len=*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: NX, NY, NZ, NT
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(2, 2, 4, NT, NX, NY, NZ), INTENT(IN) :: U_xd
      INTEGER(kind=C_INT), INTENT(IN) :: old_nproc
      INTEGER(kind=C_LONG), DIMENSION(:), INTENT(IN) :: seed
      ! for writing
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(NX, NY, NZ, NT, 4, 2) :: UHold
      COMPLEX(kind=C_DOUBLE_COMPLEX) :: UTmp
      INTEGER, PARAMETER :: infl = 107
      ! counters
      INTEGER :: it, ix, iy, iz, mu, ab, aa, bb

      ! put it in 1x2
      DO it = 1, NT
         DO iz = 1, NZ
            DO iy = 1, NY
               DO ix = 1, NX
                  UHold(ix, iy, iz, it, :, 1) = U_xd(1, 1, :, it, ix, iy, iz)
                  UHold(ix, iy, iz, it, :, 2) = U_xd(1, 2, :, it, ix, iy, iz)
               END DO
            END DO
         END DO
      END DO
      ! Shift the mu ordering
      UHold = CSHIFT(UHold, 1, dim=5)
      OPEN (infl, file=TRIM(filename), form="unformatted", access="stream", &
            status="replace", action="write", convert="little_endian")

      WRITE (infl) old_nproc
      DO ab = 1, 2
         DO mu = 1, 4
            DO it = 1, NT
               DO iz = 1, NZ
                  DO iy = 1, NY
                     DO ix = 1, NX
                        UTmp = UHold(ix, iy, iz, it, mu, ab)
                        WRITE (infl) UTmp
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      WRITE (infl) seed
      CLOSE (infl)
   END SUBROUTINE WriteGaugeField_HKLS

END MODULE FLUE_SU2_HKLS
