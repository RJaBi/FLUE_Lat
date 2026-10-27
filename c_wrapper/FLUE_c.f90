MODULE FLUE_c
   !gp!use FLUE, only: Q_Average, get_qhat, cone_cut, MultiplyMatMat, calc_mom_space_scalarD, scalarGluonProp
   USE FLUE, ONLY: ReadGaugeField_ILDG, writeGaugeField_ILDG, &
                   ReadGaugeField_OpenQCD, writeGaugeField_OpenQCD, &
                   ReadGaugefield_CSSM, &
                   ReadGaugeField_NERSC, writeGaugeField_NERSC, header_NERSC, writeGaugeField_NERSC_FullHeader, &
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

   SUBROUTINE writeGaugeField_NERSC_c(filename, U_xd, NX, NY, NZ, NT)
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
      CALL writeGaugeField_NERSC(filename, uinternal, NX, NY, NZ, NT)
    END SUBROUTINE writeGaugeField_NERSC_c

    SUBROUTINE writeGaugeField_NERSC_FullHeader_c( &
         filename, U_xd, NX, NY, NZ, NT, &
         storage_format, floating_point, datatype, checksum, &
         sequence_number, creator, creator_hardware, &
         creation_date, archive_date, ensemble_id, ensemble_label)
      CHARACTER(kind=C_CHAR, len=*), INTENT(IN) :: filename
      INTEGER(kind=C_INT), INTENT(IN) :: NX, NY, NZ, NT
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(NT, NX, NY, NZ, 4, 3, 3), INTENT(IN) :: U_xd
      CHARACTER(kind=C_CHAR, len=*), INTENT(IN) :: storage_format
      CHARACTER(kind=C_CHAR, len=*), INTENT(IN) :: floating_point
      CHARACTER(kind=C_CHAR, len=*), INTENT(IN) :: datatype
      CHARACTER(kind=C_CHAR, len=*), INTENT(IN) :: checksum
      INTEGER(kind=C_INT), INTENT(IN) :: sequence_number
      CHARACTER(kind=C_CHAR, len=*), INTENT(IN) :: creator
      CHARACTER(kind=C_CHAR, len=*), INTENT(IN) :: creator_hardware
      CHARACTER(kind=C_CHAR, len=*), INTENT(IN) :: creation_date
      CHARACTER(kind=C_CHAR, len=*), INTENT(IN) :: archive_date
      CHARACTER(kind=C_CHAR, len=*), INTENT(IN) :: ensemble_id
      CHARACTER(kind=C_CHAR, len=*), INTENT(IN) :: ensemble_label
      ! Will put all the header in here
      TYPE(header_NERSC) :: header
      COMPLEX(kind=C_DOUBLE_COMPLEX),  DIMENSION(3,3,4,NT,NX,NY,NZ) :: uinternal

      INTEGER :: mu, ac, bc

      DO ac = 1, 3
         DO bc = 1, 3
            DO mu = 1, 4
               uinternal(ac,bc,mu,:,:,:,:) = &
                    U_xd(:,:,:,:,mu,ac,bc)
            END DO
         END DO
      END DO
      ! Load the header
      header%storage_format   = storage_format
      header%floating_point   = floating_point
      header%datatype         = datatype
      header%checksum         = checksum

      header%sequence_number  = sequence_number

      header%creator          = creator
      header%creator_hardware = creator_hardware
      header%creation_date    = creation_date
      header%archive_date     = archive_date
      header%ensemble_id      = ensemble_id
      header%ensemble_label   = ensemble_label

      CALL WriteGaugeField_NERSC_FullHeader( &
           filename, uinternal, NX, NY, NZ, NT, header)

   END SUBROUTINE writeGaugeField_NERSC_FullHeader_c

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

   FUNCTION ReadGaugeField_NERSC_c(filename, NX, NY, NZ, NT) RESULT(U_xd)
      CHARACTER(kind=C_CHAR, len=*), INTENT(IN) :: filename
      INTEGER(kind=C_INT), INTENT(IN) :: NX, NY, NZ, NT
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(NT, NX, NY, NZ, 4, 3, 3) :: U_xd
      COMPLEX(kind=C_DOUBLE_COMPLEX), DIMENSION(3, 3, 4, NT, NX, NY, NZ) :: uinternal
      INTEGER :: mu, ac, bc
      uinternal = ReadGaugeField_NERSC(filename, NX, NY, NZ, NT, .TRUE.)
      DO ac = 1, 3
         DO bc = 1, 3
            DO mu = 1, 4
               U_xd(:, :, :, :, mu, ac, bc) = uinternal(ac, bc, mu, :, :, :, :)
            END DO
         END DO
      END DO
    END FUNCTION ReadGaugeField_NERSC_c

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
