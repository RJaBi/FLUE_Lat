program superMWE
  !! Entirely self-contained example showing how the AMD compiler breaks for reductions
  !! it compiles but gives wrong answer
   use ISO_C_BINDING, only: C_DOUBLE
   implicit none(external)
   integer, parameter :: WP = C_DOUBLE
   ! set up size
   integer :: NX, NY, NZ, NT
   ! counters
   integer :: nnt, nnx, nny, nnz
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
#else
      do concurrent(nnx=1:nx, nny=1:ny, nnz=1:nz, nnt=1:nt)
#endif
         sumTrP = sumTrP + 1.0_WP
         nP = nP + 1
      end do

      write (*, *) 'Done', sumTrP, nP

      end program superMWE
