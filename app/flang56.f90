program flang56
  !! Entirely self-contained example showing how the amdflang-new
  !! is broken for offload for O0 (but works for O3)
  !! it compiles but gives wrong answer
   use, intrinsic :: ISO_C_BINDING, only: C_DOUBLE
   implicit none !(external)
   integer, parameter :: WP = C_DOUBLE
   ! set up size
   integer :: NX, NY, NZ, NT
   ! counters
   integer :: nnt, nnx, nny, nnz, ll
   ! sum
   real(kind=WP) :: sumTrP
   integer :: nP
   NX = 32
   NY = 32
   NZ = 32
   NT = 8
   nP = 0
   sumTrP = 0.0_WP
   !$omp target teams distribute parallel do collapse(2) reduction(+:sumTrp,nP) private(ll,nnz,nnt) shared(nz,nt)
   do nnx = 1, nx
      do nny = 1, ny
         do nnz = 1, nz
            do nnt = 1, nt
               do ll=1, 10
                  sumTrP = sumTrP + 1.0_WP
                  nP = nP + 1
               end do
            end do
         end do
      end do
   end do

   write (*, *) 'Done', sumTrP, nP

 end program flang56

