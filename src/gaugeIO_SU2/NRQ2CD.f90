!! Functions to write NRQ2CD format

MODULE FLUE_SU2_NRQ2CD
  USE FLUE_constants, ONLY: WC
  IMPLICIT NONE(TYPE, EXTERNAL)
  PRIVATE
  PUBLIC :: writeGaugeField_NRQ2CD
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
   !------------------------------------------------------------
   ! Write to disk (raw binary, stream I/O)
   !------------------------------------------------------------
   OPEN(newunit=infl, file=filename, form="unformatted", &
        access="stream", status="replace")
   WRITE(infl) s
   CLOSE(infl)
   DEALLOCATE(s)

 END SUBROUTINE writeGaugeField_NRQ2CD

END MODULE FLUE_SU2_NRQ2CD
