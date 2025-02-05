!! Functions to read (maybe write later) ILDG binary data formats
module FLUE_ILDG_bin
   use FLUE_constants, only: WP, WC
   !use stdlib_linalg, only: det
   use ISO_FORTRAN_ENV, only: compiler_version, compiler_options, OUTPUT_UNIT
   implicit none(external)
   private
   public :: ReadGaugeField_ILDG
   public :: ReadGaugeTransformation_cola
   public :: FixSU3Matrix

contains

   function ReadGaugeField_ILDG(filename, NX, NY, NZ, NT) result(U_xd)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: NX, NY, NZ, NT
      complex(kind=WC), dimension(NT, NX, NY, NZ, 4, 3, 3) :: U_xd
      complex(kind=WC), dimension(3, 3, 4, NX, NY, NZ, NT) :: URead
      integer, parameter :: infl = 101
      integer :: matrix_len, irecl
      ! counters
      integer :: it, ix, iy, iz, mu, nu
      ! First read the gaugefield
      write (OUTPUT_UNIT, *) TRIM(filename)

      matrix_len = 16 * 3 * 3
      irecl = matrix_len * 4 * NX * NY * NZ
      ! write(*,*) matrix_len, irecl, 3, 3, 4, nx, ny, nz, nt
      open (infl, file=TRIM(filename), form='unformatted', access='direct', &
            status='old', action='read', recl=irecl, convert='BIG_ENDIAN')
      do it = 1, NT
         read (infl, rec=it) URead(:, :, :, :, :, :, it)
      end do
      close (infl)
      ! THen really dumbly re-order the indices
      do it = 1, NT
         do ix = 1, NX
            do iy = 1, NY
               do iz = 1, NZ
                  do mu = 1, 4
                     ! Mu is the openqcd index
                     ! i.e. starts t,x,y,z
                     ! nu is the ILDG index
                     ! i.e. starts x,y,z,t
                     select case (mu)
                     case (1)
                        nu = 2
                     case (2)
                        nu = 3
                     case (3)
                        nu = 4
                     case (4)
                        nu = 1
                     end select
                     U_xd(it, ix, iy, iz, nu, :, :) = TRANSPOSE(URead(:, :, mu, ix, iy, iz, it))
                     call FixSU3Matrix(U_xd(it, ix, iy, iz, nu, :, :))
                  end do
               end do
            end do
         end do
      end do

   end function ReadGaugeField_ILDG

   function ReadGaugeTransformation_cola(filename, NX, NY, NZ, NT) result(G_x)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: NX, NY, NZ, NT
      !complex(kind=WC), dimension(:,:,:,:,:,:) :: G_x
      complex(kind=WC), dimension(NT, NX, NY, NZ, 3, 3) :: G_x

      real(WP), dimension(:, :, :, :, :, :), allocatable :: ReG, ImG
      integer :: ix, iy, iz, it, ic, irank
      integer :: jx, jy, jz, jt

      complex(WC), dimension(:, :, :, :, :, :), allocatable :: G_tr
      !type(colour_vector) :: v1,v2,v3
      complex(kind=WC), dimension(3) :: v1, v2, v3

      integer, parameter :: nc = 3

      allocate (G_tr(nx, ny, nz, nt, nc, nc - 1))
      allocate (ReG(nx, ny, nz, nt, nc, nc - 1))
      allocate (ImG(nx, ny, nz, nt, nc, nc - 1))
      open (101, file=filename, form='unformatted', status='old', action='read', convert='BIG_ENDIAN')
      ! File format means we read the first two rows of G_x.
      ! Instead we read the first two columns of G_tr (the transpose of G_x) for better memory alignment.
      do ic = 1, nc - 1
         read (101) ReG(:, :, :, :, :, ic)
         read (101) ImG(:, :, :, :, :, ic)
      end do
      close (101)
      G_tr(1:nx, 1:ny, 1:nz, 1:nt, :, :) = CMPLX(ReG(1:nx, 1:ny, 1:nz, 1:nt, :, :), ImG(1:nx, 1:ny, 1:nz, 1:nt, :, :), kind=WC)
      do it = 1, nt; do iz = 1, nz; do iy = 1, ny; do ix = 1, nx
                  v1 = G_tr(ix, iy, iz, it, :, 1)
                  v2 = G_tr(ix, iy, iz, it, :, 2)
                  call orthogonalise_vectors(v2, v1)
                  call vector_product(v3, v1, v2)
                  ! G_x is the transpose of G_tr
                  G_x(it, ix, iy, iz, 1, :) = v1
                  G_x(it, ix, iy, iz, 2, :) = v2
                  G_x(it, ix, iy, iz, 3, :) = v3
               end do; end do; end do; end do
      deallocate (G_tr, ReG, ImG)

   end function ReadGaugeTransformation_cola

   ! stripped from cola and de-colour vectored
   pure subroutine orthogonalise_vectors(w, v)

      complex(kind=WC), dimension(3), intent(inout) :: w
      complex(kind=WC), dimension(3), intent(in) :: v
      complex(kind=WC) :: vdotw

      vdotw = SUM(CONJG(v) * w)

      w = w - v * vdotw

   end subroutine orthogonalise_vectors
   ! stripped from cola and de-colour vectored
   pure subroutine vector_product(x, v, w)

      complex(kind=WC), dimension(3), intent(out) :: x
      complex(kind=WC), dimension(3), intent(in) :: v, w
      integer :: ic, jc, kc

      integer, parameter :: nc = 3

      do ic = 1, nc
         jc = MODULO(ic, 3) + 1
         kc = MODULO(jc, 3) + 1

         x(ic) = CONJG(v(jc) * w(kc) - v(kc) * w(jc))
      end do

   end subroutine vector_product

   ! stripped from cola and de-colour vectored
   subroutine FixSU3Matrix(U_x)
      complex(kind=WC), dimension(3, 3), intent(inout) :: U_x
      complex(kind=WC), dimension(3) :: v1, v2, v3
      v1 = U_x(1, :)
      call normalise_vector(v1)
      v2(:) = U_x(2, :)
      call orthogonalise_vectors(v2, v1)
      call normalise_vector(v2)
      call vector_product(v3, v1, v2)
      call normalise_vector(v3)
      U_x(1, :) = v1(:)
      U_x(2, :) = v2(:)
      U_x(3, :) = v3(:)
   end subroutine FixSU3Matrix

   ! stripped from cola and de-colour vectored
   subroutine normalise_vector(v)
      complex(kind=WC), dimension(3), intent(inout) :: v
      real(WP) :: norm
      norm = SQRT(SUM(real(v)**2 + AIMAG(v)**2))
      v = v / norm
   end subroutine normalise_vector

end module FLUE_ILDG_bin
