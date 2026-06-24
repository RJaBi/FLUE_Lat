!! Functions to write NRQ2CD format

module FLUE_SU2_NRQ2CD
   use FLUE_constants, only: WC
   implicit none(type, external)
   private
   public :: writeGaugeField_NRQ2CD
   public :: readGaugeField_NRQ2CD
contains
   subroutine writeGaugeField_NRQ2CD(filename, NX, NY, NZ, NT, U)
      character(len=*), intent(IN) :: filename
      integer, intent(IN) :: NX, NY, NZ, NT
      complex(kind=WC), dimension(2, 2, 4, NT, NX, NY, NZ), intent(IN) :: U
      integer :: ix, iy, iz, it, mu
      integer :: ix2, ip, par
      integer :: c1, c2, infl
      integer :: ip5d2
      complex(kind=WC), dimension(:, :, :, :, :), allocatable :: s
      ip5d2 = INT(NX * NY * NZ * NT * 0.5)
      allocate (s(ip5d2, 2, 2, 4, 0:1))
      write (*, *) 'a', U(:, :, 2, 1, 1, 3, 1)
      !------------------------------------------------------------
      ! Rearrangement
      !------------------------------------------------------------
      do it = 1, NT
         do iz = 1, NZ
            do iy = 1, NY
               do ix = 1, NX
               !! get parity
                  par = MOD(ix + iy + iz + it, 2)
               !! get x loc
                  ix2 = (ix - 1) / 2
               !! nx/2 is integer division
               !! Get linear mapping
                  ip = ((((it - 1) * nz + (iz - 1)) * ny + (iy - 1)) * nx / 2) &
                       + ix2 + 1
                  do mu = 1, 4
                     do c2 = 1, 2
                        do c1 = 1, 2
                           s(ip, c1, c2, mu, par) = U(c1, c2, mu, it, ix, iy, iz)
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
      ! Shift the mu-ordering so that xyzt
      s = CSHIFT(s, 1, dim=4)
      write (*, *) 'b25o', s(25, :, :, 1, 0)
      write (*, *) 'b25e', s(25, :, :, 1, 1)
      !------------------------------------------------------------
      ! Write to disk (raw binary, stream I/O)
      !------------------------------------------------------------
      open (newunit=infl, file=filename, form="unformatted", &
            access="stream", status="replace")
      write (infl) s
      close (infl)
      deallocate (s)
   end subroutine writeGaugeField_NRQ2CD

   subroutine readGaugeField_NRQ2CD(filename, NX, NY, NZ, NT, U)
      character(len=*), intent(IN) :: filename
      integer, intent(IN) :: NX, NY, NZ, NT
      complex(kind=WC), dimension(2, 2, 4, NT, NX, NY, NZ), intent(OUT) :: U
      integer :: ix, iy, iz, it, mu
      integer :: ix2, ip, par
      integer :: infl, ip5d2
      complex(kind=WC), dimension(:, :, :, :, :), allocatable :: s
      ! Optional but sensible: the mapping assumes NX is even
      if (MOD(NX, 2) /= 0) then
         write (*, *) "readGaugeField_NRQ2CD: NX must be even.", NX
         stop
      end if
      ip5d2 = NX * NY * NZ * NT / 2
      allocate (s(ip5d2, 2, 2, 4, 0:1))
      !------------------------------------------------------------
      ! Read raw binary block
      !------------------------------------------------------------
      open (newunit=infl, file=filename, form="unformatted", &
            access="stream", status="old", action="read")
      read (infl) s
      close (infl)
      !------------------------------------------------------------
      ! Inverse rearrangement
      !------------------------------------------------------------
      do it = 1, NT
         do iz = 1, NZ
            do iy = 1, NY
               do ix = 1, NX
                  par = MOD(ix + iy + iz + it, 2)
                  ix2 = (ix - 1) / 2
                  ip = ((((it - 1) * NZ + (iz - 1)) * NY + (iy - 1)) * (NX / 2)) + ix2 + 1
                  do mu = 1, 4
                     U(:, :, mu, it, ix, iy, iz) = s(ip, :, :, mu, par)
                  end do
               end do
            end do
         end do
      end do
      ! Shift the mu-index so that txyz
      U = CSHIFT(U, -1, dim=3)
      deallocate (s)
   end subroutine readGaugeField_NRQ2CD

end module FLUE_SU2_NRQ2CD
