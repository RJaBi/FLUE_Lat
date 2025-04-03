program flangArrays
  !! Entirely self-contained example showing how the amdflang-new
  !! is broken for offload for O0 (but works for O3)
  !! It does not compile
  !! See bottom of file
   use, intrinsic :: ISO_C_BINDING, only: C_DOUBLE
   implicit none !(external)
   integer, parameter :: WP = C_DOUBLE
   ! set up size
   integer :: NX
   ! counters
   integer :: nnx, ll
   ! arrays
   integer, dimension(4) :: ar1
   integer, dimension(4) :: ar2
   integer, dimension(4,4) :: arT
   ! sum
   real(kind=WP) :: sumTrP
   integer :: nP
   NX = 32
   nP = 0
   ar2 = (/5, 20, 10, 1/)
   arT = 1.0_WP
   sumTrP = 0.0_WP
   !$omp target teams distribute parallel do reduction(+:sumTrp,nP) private(ll,ar1) shared(ar2)
   do nnx = 1, nx
      do ll=1, 10
         ar1 = ar2
         nP = nP + sum(ar1)
         arT(1,2) = nP
         arT = transpose(arT)
         sumTrP = sumTrP + arT(2,1)
      end do
   end do

   write (*, *) 'Done', sumTrP, nP

 end program flangArrays



!
! ld.lld: error: undefined symbol: _FortranAAssign
!>>> referenced by flangAssignA.f90:28
!>>>               /tmp/a.out.amdgcn.gfx942-d6effa.img.lto.o:(__omp_offloading_32_73374__QQmain_l25..omp_par)
!>>> referenced by flangAssignA.f90:28
!>>>               /tmp/a.out.amdgcn.gfx942-d6effa.img.lto.o:(__omp_offloading_32_73374__QQmain_l25..omp_par)
!
!ld.lld: error: undefined symbol: _FortranASumInteger4
!>>> referenced by flangAssignA.f90:29
!>>>               /tmp/a.out.amdgcn.gfx942-d6effa.img.lto.o:(__omp_offloading_32_73374__QQmain_l25..omp_par)
!>>> referenced by flangAssignA.f90:29
!>>>               /tmp/a.out.amdgcn.gfx942-d6effa.img.lto.o:(__omp_offloading_32_73374__QQmain_l25..omp_par)
!
!ld.lld: error: undefined symbol: _FortranATranspose
!>>> referenced by flangAssignA.f90:31
!>>>               /tmp/a.out.amdgcn.gfx942-d6effa.img.lto.o:(__omp_offloading_32_73374__QQmain_l25..omp_par)
!>>> referenced by flangAssignA.f90:31
!>>>               /tmp/a.out.amdgcn.gfx942-d6effa.img.lto.o:(__omp_offloading_32_73374__QQmain_l25..omp_par)
!clang: error: ld.lld command failed with exit code 1 (use -v to see invocation)
!/opt/rocmplus-6.3.3/rocm-afar-7110-drop-5.3.0/lib/llvm/bin/clang-linker-wrapper: error: 'clang' failed
!flang-new: error: linker command failed with exit code 1 (use -v to see invocation)
