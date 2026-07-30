module FLUE_openQCDFileIO_SA
  !< Functions to read and write openqcd format gaugefields
   use FLUE_constants, only: WP, WC
   use FLUE_SU3MatrixOps, only: FixSU3Matrix
   use FLUE_wloops, only: genPlaquette
   implicit none(type, external)
   private
   public :: ReadGaugeField_OpenQCD
   public :: writeGaugeField_OpenQCD

contains

  elemental function modc(a, b) result(c)
    !< Maps an integer a to the set of integers [1,b] i.e. positive integers with cycle length b.
    !! Taken directly from COLA
    integer, intent(IN) :: a, b
    !< input variables
    integer :: c
    !< result
    !c = a - ((a-1)/b)*b
    c = MODULO(a - 1, b) + 1
   end function modc

   function ReadGaugeField_OpenQCD(filename, NX, NY, NZ, NT, fixSU3) result(U_out)
     !< Reads an openqcd format gaugefield into internal form [3,3,4,NT,NX,NY,NZ]
     character(len=*), intent(IN) :: filename
     !< The file to be read
     integer, intent(IN) :: NX, NY, NZ, NT
     !< Dimensions of gaugefield
     logical, optional, intent(IN) :: fixSU3
     !< Optionally reproject each link to SU3
     complex(kind=WC), dimension(NT, NX, NY, NZ, 4, 3, 3) :: U_xd
     !< For reading
     complex(kind=WC), dimension(3, 3, 4, NT, NX, NY, NZ) :: U_out
     !< Output gaugefield
      !complex(kind=WP), dimension(:,:,:,:,:,:,:), allocatable :: U_xd

     complex(kind=WC), dimension(3, 3) :: UTmp
     !< SU3 matrix. Read each link individually
     integer, parameter :: infl = 107
     !< read from this file unit
      ! Header info
     real(kind=WP) :: plaq
     !< plaquette in the header
     integer :: ntdim, nxdim, nydim, nzdim
     !< lattice dimensions in header
      ! counters
     integer :: it, ix, iy, iz, mu, id
     !< counters 1
     integer :: jx, jy, jz, jt
     !< counters 2
     integer, dimension(4) :: dmu
     !< For backwards/forwards links
     logical :: fixSU3Set
     !< whether we are doing SU3 re-projection

      if (PRESENT(fixSU3)) then
         fixSU3Set = fixSU3
      else
         fixSU3Set = .TRUE.
      end if

      write (*, *) "here ", TRIM(filename)
      open (infl, file=TRIM(filename), form="unformatted", access="stream", status="old", action="read", convert="little_endian")
      read (infl) ntdim, nxdim, nydim, nzdim, plaq

      !allocate(U_xd(nxdim,nydim,nzdim,ntdim,4,3,3))

      ! z varies quickest, then y, then x, then t
      do it = 1, ntdim
         do ix = 1, nxdim
            do iy = 1, nydim
               do iz = 1, nzdim
                  if (MODULO(ix + iy + iz + it - 4, 2) == 0) cycle  ! Format only considers odd points

                  do id = 1, 4
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
                     read (infl) UTmp
                     U_xd(it, ix, iy, iz, mu, :, :) = UTmp
                     read (infl) UTmp
                     U_xd(jt, jx, jy, jz, mu, :, :) = UTmp

                     U_xd(it, ix, iy, iz, mu, :, :) = TRANSPOSE(U_xd(it, ix, iy, iz, mu, :, :))
                     U_xd(jt, jx, jy, jz, mu, :, :) = TRANSPOSE(U_xd(jt, jx, jy, jz, mu, :, :))
                     if (fixSU3Set) then
                        ! Welll FixSU3Matrix did nothing to the average plaquette value
                        UTmp = U_xd(it, ix, iy, iz, mu, :, :)
                        call FixSU3Matrix(UTmp)
                        U_xd(it, ix, iy, iz, mu, :, :) = UTmp
                        UTmp = U_xd(jt, jx, jy, jz, mu, :, :)
                        call FixSU3Matrix(UTmp)
                        U_xd(jt, jx, jy, jz, mu, :, :) = UTmp
                     end if
                  end do
               end do
            end do
         end do
      end do

      close (infl)

      U_xd = CSHIFT(U_xd, -1, dim=5)
      do it = 1, 3
         do iz = 1, 3
            do mu = 1, 4
               U_out(it, iz, mu, :, :, :, :) = U_xd(:, :, :, :, mu, it, iz)
            end do
         end do
      end do

   end function ReadGaugeField_OpenQCD

   subroutine writeGaugeField_OpenQCD(filename, U_in, NX, NY, NZ, NT)
     !< Writes an openqcd gaugefield to file
     character(len=*), intent(IN) :: filename
     !< file to write to
     complex(kind=WC), dimension(3, 3, 4, NT, NX, NY, NZ), intent(IN) :: U_in
     !< gaugefield to write
     integer, intent(IN) :: NX, NY, NZ, NT
     !< lattice dimensions
     complex(kind=WC), dimension(NT, NX, NY, NZ, 4, 3, 3) :: U
     !< reordered lattice
     complex(kind=WC), dimension(3, 3) :: UTmp
     !< single link writes
     integer, parameter :: infl = 107
     !< use this file unit
      ! Header info
     real(kind=WP) :: plaq, sumTrP, time
     !< plaquette calc variables
     integer :: NP
     !< number of plaquettes
     integer :: it, ix, iy, iz, mu, id
     !< counters 1
     integer :: jx, jy, jz, jt
     !< counters 2
     integer, dimension(4) :: dmu
     !< backwards/forwards
     logical :: fixSU3Set
     !< whether we are doing SU3 re-projection
      do it = 1, 3
         do iz = 1, 3
            do mu = 1, 4
               U(:, :, :, :, mu, it, iz) = U_in(it, iz, mu, :, :, :, :)
            end do
         end do
      end do

      ! Calculate the plaquette as needed by oqcd header
      call genPlaquette(U_in, NT, NX, NY, NZ, 1, 4, 4, sumTrp, NP, time)
      plaq = sumTrp / real(NP, kind=WC)

      !write (*, *) 'here ', TRIM(filename)
      open (infl, file=TRIM(filename), form="unformatted", access="stream", &
            status="replace", action="write", convert="little_endian")
      write (infl) nt, nx, ny, nz, plaq

      !allocate(U(nxdim,nydim,nzdim,ntdim,4,3,3))

      U = CSHIFT(U, 1, dim=5)

      ! z varies quickest, then y, then x, then t
      do it = 1, nt
         do ix = 1, nx
            do iy = 1, ny
               do iz = 1, nz
                  if (MODULO(ix + iy + iz + it - 4, 2) == 0) cycle  ! Format only considers odd points

                  do id = 1, 4
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
                     write (infl) UTmp
                     UTmp = U(jt, jx, jy, jz, mu, :, :)
                     write (infl) UTmp

                  end do
               end do
            end do
         end do
      end do

      close (infl)

   end subroutine WriteGaugeField_OpenQCD

end module FLUE_openQCDFileIO_SA
