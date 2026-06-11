!! Functions to write NRQ2CD format

MODULE FLUE_SU2_NRQ2CD
  USE FLUE_constants, ONLY: WC
  IMPLICIT NONE(TYPE, EXTERNAL)
  PRIVATE
  PUBLIC :: writeGaugeField_NRQ2CD
  PUBLIC :: readGaugeField_NRQ2CD
CONTAINS
  SUBROUTINE writeGaugeField_NRQ2CD(filename, NX, NY, NZ, NT, U)
    CHARACTER(len=*), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: NX, NY, NZ, NT
    COMPLEX(kind=WC), DIMENSION(2,2,4, NT, NX, NY, NZ), INTENT(IN) :: U
    INTEGER :: ix, iy, iz, it, mu
    INTEGER :: ix2, ip, par
    INTEGER :: c1, c2, infl
    INTEGER :: ip5d2
   COMPLEX(kind=WC), DIMENSION(:,:,:,:,:), ALLOCATABLE :: s
   ip5d2 = int(NX*NY*NZ*NT*0.5)
   ALLOCATE(s(ip5d2, 2, 2, 4, 0:1))
   WRITE(*,*) 'a', U(:,:,2,1,1,3,1)
   !------------------------------------------------------------
   ! Rearrangement
   !------------------------------------------------------------
   DO it = 1, NT
      DO iz = 1, NZ
         DO iy = 1, NY
            DO ix = 1, NX
               !! get parity
               par = mod(ix + iy + iz + it, 2)
               !! get x loc
               ix2    = (ix - 1) / 2
               !! nx/2 is integer division
               !! Get linear mapping
               ip = ((( (it-1)*nz + (iz-1) )*ny + (iy-1) )*nx/2) &
                    + ix2 + 1
               DO mu = 1, 4
                  DO c2 = 1, 2
                     DO c1 = 1, 2
                        s(ip, c1, c2, mu, par) =  U(c1, c2, mu, it, ix, iy, iz)
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
   END DO
   ! Shift the mu-ordering so that xyzt
   s = cshift(s, 1, dim=4)
   WRITE(*,*) 'b25o', s(25,:,:,1,0)
   WRITE(*,*) 'b25e', s(25,:,:,1,1)
   !------------------------------------------------------------
   ! Write to disk (raw binary, stream I/O)
   !------------------------------------------------------------
   OPEN(newunit=infl, file=filename, form="unformatted", &
        access="stream", status="replace")
   WRITE(infl) s
   CLOSE(infl)
   DEALLOCATE(s)
 END SUBROUTINE writeGaugeField_NRQ2CD


 SUBROUTINE readGaugeField_NRQ2CD(filename, NX, NY, NZ, NT, U)
   CHARACTER(len=*), INTENT(IN) :: filename
   INTEGER, INTENT(IN) :: NX, NY, NZ, NT
   COMPLEX(kind=WC), DIMENSION(2,2,4,NT,NX,NY,NZ), INTENT(OUT) :: U
   INTEGER :: ix, iy, iz, it, mu
   INTEGER :: ix2, ip, par
   INTEGER :: infl, ip5d2
   COMPLEX(kind=WC), DIMENSION(:,:,:,:,:), ALLOCATABLE :: s
   ! Optional but sensible: the mapping assumes NX is even
   IF (MOD(NX, 2) /= 0) THEN
      WRITE(*,*) "readGaugeField_NRQ2CD: NX must be even.", NX
      STOP
   END IF
   ip5d2 = NX*NY*NZ*NT/2
   ALLOCATE(s(ip5d2, 2, 2, 4, 0:1))
   !------------------------------------------------------------
   ! Read raw binary block
   !------------------------------------------------------------
   OPEN(newunit=infl, file=filename, form="unformatted", &
        access="stream", status="old", action="read")
   READ(infl) s
   CLOSE(infl)
   !------------------------------------------------------------
   ! Inverse rearrangement
   !------------------------------------------------------------
   DO it = 1, NT
      DO iz = 1, NZ
         DO iy = 1, NY
            DO ix = 1, NX
               par = MOD(ix + iy + iz + it, 2)
               ix2 = (ix - 1) / 2
               ip  = ((((it-1)*NZ + (iz-1))*NY + (iy-1)) * (NX/2)) + ix2 + 1
               DO mu = 1, 4
                  U(:,:,mu,it,ix,iy,iz) = s(ip,:,:,mu,par)
               END DO
            END DO
         END DO
      END DO
   END DO
   ! Shift the mu-index so that txyz
   U = cshift(U, -1, dim=3)
   DEALLOCATE(s)
 END SUBROUTINE readGaugeField_NRQ2CD



END MODULE FLUE_SU2_NRQ2CD
