MODULE FLUE_c
   !gp!use FLUE, only: Q_Average, get_qhat, cone_cut, MultiplyMatMat, calc_mom_space_scalarD, scalarGluonProp
   USE FLUE, ONLY: ReadGaugeField_ILDG, writeGaugeField_ILDG, &
                   ReadGaugeField_OpenQCD, writeGaugeField_OpenQCD, &
                   ReadGaugefield_CSSM, &
                   genplaquette, StoutSmearLinks
  !! Note that the c_wrapper and python uses order [NT, NX, NY, NZ, mu, colour, colour]
  !! This contrasts the internal ordering of [colour, colour, mu, NT, NX, NY, NZ]
  !! This is for compatibility with i.e. lyncs_io
   USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_DOUBLE, C_INT, C_DOUBLE_COMPLEX, C_CHAR
   IMPLICIT NONE(EXTERNAL)
   PUBLIC

CONTAINS

   !gaugeIO
   ! Writers
   SUBROUTINE writeGaugeField_ILDG_c(filename, U_xd, NX, NY, NZ, NT)
      CHARACTER(kind=C_CHAR, len=*), INTENT(IN) :: filename
      INTEGER(kind=C_INT), INTENT(IN) :: NX, NY, NZ, NT
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(NT, NX, NY, NZ, 4, 3, 3), INTENT(IN) :: U_xd
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(3, 3, 4, NT, NX, NY, NZ) :: uinternal
      INTEGER :: mu, ac, bc
      DO ac = 1, 3
         DO bc = 1, 3
            DO mu = 1, 4
               uinternal(ac, bc, mu, :, :, :, :) = U_xd(:, :, :, :, mu, ac, bc)
            END DO
         END DO
      END DO
      CALL writeGaugeField_ILDG(filename, uinternal, NX, NY, NZ, NT, .TRUE.)
   END SUBROUTINE writeGaugeField_ILDG_c

   SUBROUTINE writeGaugeField_OpenQCD_c(filename, U_xd, NX, NY, NZ, NT)
      CHARACTER(kind=C_CHAR, len=*), INTENT(IN) :: filename
      INTEGER(kind=C_INT), INTENT(IN) :: NX, NY, NZ, NT
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(NT, NX, NY, NZ, 4, 3, 3), INTENT(IN) :: U_xd
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(3, 3, 4, NT, NX, NY, NZ) :: uinternal
      INTEGER :: mu, ac, bc
      DO ac = 1, 3
         DO bc = 1, 3
            DO mu = 1, 4
               uinternal(ac, bc, mu, :, :, :, :) = U_xd(:, :, :, :, mu, ac, bc)
            END DO
         END DO
      END DO
      CALL writeGaugeField_OpenQCD(filename, uinternal, NX, NY, NZ, NT)
   END SUBROUTINE writeGaugeField_OpenQCD_c

   ! Readers
   FUNCTION ReadGaugeField_ILDG_c(filename, NX, NY, NZ, NT) RESULT(U_xd)
      CHARACTER(kind=C_CHAR, len=*), INTENT(IN) :: filename
      INTEGER(kind=C_INT), INTENT(IN) :: NX, NY, NZ, NT
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(NT, NX, NY, NZ, 4, 3, 3) :: U_xd
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(3, 3, 4, NT, NX, NY, NZ) :: uinternal
      INTEGER :: mu, ac, bc
      uinternal = ReadGaugeField_ILDG(filename, NX, NY, NZ, NT, .TRUE.)
      DO ac = 1, 3
         DO bc = 1, 3
            DO mu = 1, 4
               U_xd(:, :, :, :, mu, ac, bc) = uinternal(ac, bc, mu, :, :, :, :)
            END DO
         END DO
      END DO
   END FUNCTION ReadGaugeField_ILDG_c

   FUNCTION ReadGaugeField_OpenQCD_c(filename, NX, NY, NZ, NT) RESULT(U_xd)
      CHARACTER(kind=C_CHAR, len=*), INTENT(IN) :: filename
      INTEGER(kind=C_INT), INTENT(IN) :: NX, NY, NZ, NT
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(NT, NX, NY, NZ, 4, 3, 3) :: U_xd
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(3, 3, 4, NT, NX, NY, NZ) :: uinternal
      INTEGER :: mu, ac, bc
      uinternal = ReadGaugeField_OpenQCD(filename, NX, NY, NZ, NT, .TRUE.)
      DO ac = 1, 3
         DO bc = 1, 3
            DO mu = 1, 4
               U_xd(:, :, :, :, mu, ac, bc) = uinternal(ac, bc, mu, :, :, :, :)
            END DO
         END DO
      END DO
   END FUNCTION ReadGaugeField_OpenQCD_c

   FUNCTION ReadGaugeField_CSSM_c(filename, NX, NY, NZ, NT) RESULT(U_xd)
      CHARACTER(kind=C_CHAR, len=*), INTENT(IN) :: filename
      INTEGER(kind=C_INT), INTENT(IN) :: NX, NY, NZ, NT
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(NT, NX, NY, NZ, 4, 3, 3) :: U_xd
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(3, 3, 4, NT, NX, NY, NZ) :: uinternal
      INTEGER :: mu, ac, bc
      uinternal = ReadGaugeField_CSSM(filename, NX, NY, NZ, NT, .TRUE.)
      DO ac = 1, 3
         DO bc = 1, 3
            DO mu = 1, 4
               U_xd(:, :, :, :, mu, ac, bc) = uinternal(ac, bc, mu, :, :, :, :)
            END DO
         END DO
      END DO
   END FUNCTION ReadGaugeField_CSSM_c

   ! wloops
   SUBROUTINE genplaquette_c(data, nt, nx, ny, nz, mustart, muend, nuend, sumtrp, np, time)
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(nt, nx, ny, nz, 4, 3, 3), INTENT(IN) :: data
      INTEGER(kind=C_INT), INTENT(IN) :: mustart, muend, nuend
      INTEGER(kind=C_INT), INTENT(IN) :: nt, nx, ny, nz
      REAL(kind=C_DOUBLE), INTENT(OUT) :: sumtrp, time
      INTEGER(kind=C_INT), INTENT(OUT) :: np
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(3, 3, 4, NT, NX, NY, NZ) :: uinternal
      INTEGER :: mu, ac, bc
      DO ac = 1, 3
         DO bc = 1, 3
            DO mu = 1, 4
               uinternal(ac, bc, mu, :, :, :, :) = data(:, :, :, :, mu, ac, bc)
            END DO
         END DO
      END DO
      CALL genplaquette(uinternal, nt, nx, ny, nz, mustart, muend, nuend, sumtrp, np, time)
   END SUBROUTINE genplaquette_c

   SUBROUTINE stoutsmearlinks_c(data, rho, nSweeps, nt, nx, ny, nz, usmeared)
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(nt, nx, ny, nz, 4, 3, 3), INTENT(IN) :: data
      REAL(kind=C_DOUBLE), INTENT(IN) :: rho
      INTEGER(kind=C_INT), INTENT(IN) :: nSweeps
      INTEGER(kind=C_INT), INTENT(IN) :: nx, ny, nz, nt
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(nt, nx, ny, nz, 4, 3, 3), INTENT(OUT) :: usmeared
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(3, 3, 4, NT, NX, NY, NZ) :: uinternal, uinternalsmeared
      INTEGER :: mu, ac, bc
      DO ac = 1, 3
         DO bc = 1, 3
            DO mu = 1, 4
               uinternal(ac, bc, mu, :, :, :, :) = data(:, :, :, :, mu, ac, bc)
            END DO
         END DO
      END DO
      CALL StoutSmearLinks(uinternal, rho, nsweeps, uinternalsmeared)
      DO ac = 1, 3
         DO bc = 1, 3
            DO mu = 1, 4
               usmeared(:, :, :, :, mu, ac, bc) = uinternalsmeared(ac, bc, mu, :, :, :, :)
            END DO
         END DO
      END DO
   END SUBROUTINE stoutsmearlinks_c

END MODULE FLUE_c
