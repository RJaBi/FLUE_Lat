!! Functions to read & write openQCD format gaugefields
MODULE FLUE_openQCDFileIO_SA
   USE FLUE_constants, ONLY: WP, WC
   USE FLUE_SU3MatrixOps, ONLY: FixSU3Matrix
   USE FLUE_wloops, ONLY: genPlaquette
   IMPLICIT NONE(TYPE, EXTERNAL)
   PRIVATE
   PUBLIC :: ReadGaugeField_OpenQCD
   PUBLIC :: writeGaugeField_OpenQCD

CONTAINS

  !! Maps an integer a to the set of integers [1,b] i.e. positive integers with cycle length b.
   ELEMENTAL FUNCTION modc(a, b) RESULT(c)
    !! Taken directly from COLA
      INTEGER, INTENT(IN) :: a, b
      INTEGER :: c
      !c = a - ((a-1)/b)*b
      c = MODULO(a - 1, b) + 1
   END FUNCTION modc

!  function determinant(matrix) result(det)
!    complex(c_double_complex), dimension(3,3), intent(in) :: matrix
!    complex(c_double_complex) :: det
!
!    det = matrix(1,1)*(matrix(2,2)*matrix(3,3) - matrix(2,3)*matrix(3,2)) - &
!         matrix(1,2)*(matrix(2,1)*matrix(3,3) - matrix(2,3)*matrix(3,1)) + &
!         matrix(1,3)*(matrix(2,1)*matrix(3,2) - matrix(2,2)*matrix(3,1))
!  end function determinant

   FUNCTION ReadGaugeField_OpenQCD(filename, NX, NY, NZ, NT, fixSU3) RESULT(U_out)
      CHARACTER(len=*), INTENT(IN) :: filename
      INTEGER, INTENT(IN) :: NX, NY, NZ, NT
      LOGICAL, OPTIONAL, INTENT(IN) :: fixSU3
      COMPLEX(kind=WC), DIMENSION(NT, NX, NY, NZ, 4, 3, 3) :: U_xd
      COMPLEX(kind=WC), DIMENSION(3, 3, 4, NT, NX, NY, NZ) :: U_out
      !complex(kind=WP), dimension(:,:,:,:,:,:,:), allocatable :: U_xd

      COMPLEX(kind=WC), DIMENSION(3, 3) :: UTmp
      INTEGER, PARAMETER :: infl = 107
      ! Header info
      REAL(kind=WP) :: plaq
      INTEGER :: ntdim, nxdim, nydim, nzdim
      ! counters
      INTEGER :: it, ix, iy, iz, mu, id
      INTEGER :: jx, jy, jz, jt
      INTEGER, DIMENSION(4) :: dmu
      LOGICAL :: fixSU3Set

      IF (PRESENT(fixSU3)) THEN
         fixSU3Set = fixSU3
      ELSE
         fixSU3Set = .TRUE.
      END IF

      WRITE (*, *) "here ", TRIM(filename)
      OPEN (infl, file=TRIM(filename), form="unformatted", access="stream", status="old", action="read", convert="little_endian")
      READ (infl) ntdim, nxdim, nydim, nzdim, plaq

      !allocate(U_xd(nxdim,nydim,nzdim,ntdim,4,3,3))

      ! z varies quickest, then y, then x, then t
      DO it = 1, ntdim
         DO ix = 1, nxdim
            DO iy = 1, nydim
               DO iz = 1, nzdim
                  IF (MODULO(ix + iy + iz + it - 4, 2) == 0) CYCLE  ! Format only considers odd points

                  DO id = 1, 4
                     mu = modc(id - 1, 4)  ! Time dimension first: mu = 4, 1, 2, 3

                     dmu(:) = 0
                     dmu(mu) = 1

                     ! Get the backward site under periodic boundary conditions
                     jx = modc(ix - dmu(1), nxdim)
                     jy = modc(iy - dmu(2), nydim)
                     jz = modc(iz - dmu(3), nzdim)
                     jt = modc(it - dmu(4), ntdim)
                     ! Read the forward and backward links in mu direction
                     !read(infl) U_g(mu,ix,iy,iz,it)%cl(:,:)
                     !read(infl) U_g(mu,jx,jy,jz,jt)%cl(:,:)
                     READ (infl) UTmp
                     U_xd(it, ix, iy, iz, mu, :, :) = UTmp
                     READ (infl) UTmp
                     U_xd(jt, jx, jy, jz, mu, :, :) = UTmp

                     U_xd(it, ix, iy, iz, mu, :, :) = TRANSPOSE(U_xd(it, ix, iy, iz, mu, :, :))
                     U_xd(jt, jx, jy, jz, mu, :, :) = TRANSPOSE(U_xd(jt, jx, jy, jz, mu, :, :))
                     IF (fixSU3Set) THEN
                        ! Welll FixSU3Matrix did nothing to the average plaquette value
                        UTmp = U_xd(it, ix, iy, iz, mu, :, :)
                        CALL FixSU3Matrix(UTmp)
                        U_xd(it, ix, iy, iz, mu, :, :) = UTmp
                        UTmp = U_xd(jt, jx, jy, jz, mu, :, :)
                        CALL FixSU3Matrix(UTmp)
                        U_xd(jt, jx, jy, jz, mu, :, :) = UTmp
                     END IF
                  END DO
               END DO
            END DO
         END DO
      END DO

      CLOSE (infl)

      U_xd = CSHIFT(U_xd, -1, dim=5)
      DO it = 1, 3
         DO iz = 1, 3
            DO mu = 1, 4
               U_out(it, iz, mu, :, :, :, :) = U_xd(:, :, :, :, mu, it, iz)
            END DO
         END DO
      END DO

   END FUNCTION ReadGaugeField_OpenQCD

   SUBROUTINE writeGaugeField_OpenQCD(filename, U_in, NX, NY, NZ, NT)
      CHARACTER(len=*), INTENT(IN) :: filename
      COMPLEX(kind=WC), DIMENSION(3, 3, 4, NT, NX, NY, NZ), INTENT(IN) :: U_in
      INTEGER, INTENT(IN) :: NX, NY, NZ, NT
      COMPLEX(kind=WC), DIMENSION(NT, NX, NY, NZ, 4, 3, 3) :: U

      COMPLEX(kind=WC), DIMENSION(3, 3) :: UTmp
      INTEGER, PARAMETER :: infl = 107
      ! Header info
      REAL(kind=WP) :: plaq, sumTrP, time
      INTEGER :: NP
      ! counters
      INTEGER :: it, ix, iy, iz, mu, id
      INTEGER :: jx, jy, jz, jt
      INTEGER, DIMENSION(4) :: dmu
      LOGICAL :: fixSU3Set

      DO it = 1, 3
         DO iz = 1, 3
            DO mu = 1, 4
               U(:, :, :, :, mu, it, iz) = U_in(it, iz, mu, :, :, :, :)
            END DO
         END DO
      END DO

      ! Calculate the plaquette as needed by oqcd header
      CALL genPlaquette(U_in, NT, NX, NY, NZ, 1, 4, 4, sumTrp, NP, time)
      plaq = sumTrp / real(NP, kind=WC)

      !write (*, *) 'here ', TRIM(filename)
      OPEN (infl, file=TRIM(filename), form="unformatted", access="stream", &
            status="replace", action="write", convert="little_endian")
      WRITE (infl) nt, nx, ny, nz, plaq

      !allocate(U(nxdim,nydim,nzdim,ntdim,4,3,3))

      U = CSHIFT(U, 1, dim=5)

      ! z varies quickest, then y, then x, then t
      DO it = 1, nt
         DO ix = 1, nx
            DO iy = 1, ny
               DO iz = 1, nz
                  IF (MODULO(ix + iy + iz + it - 4, 2) == 0) CYCLE  ! Format only considers odd points

                  DO id = 1, 4
                     mu = modc(id - 1, 4)  ! Time dimension first: mu = 4, 1, 2, 3

                     dmu(:) = 0
                     dmu(mu) = 1

                     ! Get the backward site under periodic boundary conditions
                     jx = modc(ix - dmu(1), nx)
                     jy = modc(iy - dmu(2), ny)
                     jz = modc(iz - dmu(3), nz)
                     jt = modc(it - dmu(4), nt)
                     ! Read the forward and backward links in mu direction
                     !read(infl) U_g(mu,ix,iy,iz,it)%cl(:,:)
                     !read(infl) U_g(mu,jx,jy,jz,jt)%cl(:,:)

                     U(it, ix, iy, iz, mu, :, :) = TRANSPOSE(U(it, ix, iy, iz, mu, :, :))
                     U(jt, jx, jy, jz, mu, :, :) = TRANSPOSE(U(jt, jx, jy, jz, mu, :, :))

                     UTmp = U(it, ix, iy, iz, mu, :, :)
                     WRITE (infl) UTmp
                     UTmp = U(jt, jx, jy, jz, mu, :, :)
                     WRITE (infl) UTmp

                  END DO
               END DO
            END DO
         END DO
      END DO

      CLOSE (infl)

   END SUBROUTINE WriteGaugeField_OpenQCD

END MODULE FLUE_openQCDFileIO_SA
