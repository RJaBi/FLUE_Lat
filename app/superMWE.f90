program superMWE
  !! Entirely self-contained example showing how the AMD compiler breaks for reductions
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
#ifdef AMD
   do concurrent(nnx=1:nx, nny=1:ny, nnz=1:nz, nnt=1:nt) &
      reduce(+:sumTrP, nP)
#elif OMP
      !$omp parallel do collapse(4) reduction(+:sumTrP,nP) private(ll)
      do nnx = 1, nx
         do nny = 1, ny
            do nnz = 1, nz
               do nnt = 1, nt
#elif OMPOff
      !$omp target teams distribute parallel do collapse(4) reduction(+:sumTrp,nP) private(ll)
       do nnx = 1, nx
          do nny = 1, ny
             do nnz = 1, nz
                do nnt = 1, nt
#elif OMPOffLoop
      !$omp target teams loop collapse(4) reduction(+:sumTrp,nP) private(ll)
       do nnx = 1, nx
          do nny = 1, ny
             do nnz = 1, nz
                do nnt = 1, nt
#else
                  do concurrent(nnx=1:nx, nny=1:ny, nnz=1:nz, nnt=1:nt)
#endif
                     do ll=1, 1000
                        sumTrP = sumTrP + 1.0_WP
                        nP = nP + 1
                     end do
#ifdef OMP
                  end do
               end do
            end do
         end do
#elif OMPOff
                  end do
               end do
            end do
         end do
#elif OMPOffLoop
                  end do
               end do
            end do
         end do 
#else
      end do
#endif

      write (*, *) 'Done', sumTrP, nP

      end program superMWE
