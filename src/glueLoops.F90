!! per coordinate/plane wilson loops
module FLUE_glueLoops

   use FLUE_constants, only: WP, WC
   use FLUE_wloops, only: genericPath

   implicit none(external)

   public

   abstract interface
      pure function coordMuNuInterface(data, NT, NX, NY, NZ, mu, nu) result(U_xd)
         import :: WC
         implicit none(external)
         complex(kind=WC), dimension(NT, NX, NY, NZ, 4, 3, 3), intent(in) :: data
         integer, intent(in) :: NT, NX, NY, NZ, mu, nu
         complex(kind=WC), dimension(NT, NX, NY, NZ, 3, 3) :: U_xd
      end function coordMuNuInterface
   end interface

contains

   pure function genSpaceAverageTr(coordMuNuFunc, data, NT, NX, NY, NZ, mu, nu) result(genTau)
      ! Traces over the function coordMuNuFunc for each coord
      ! Then averages over space
      complex(kind=WC), dimension(NT, NX, NY, NZ, 4, 3, 3), intent(in) :: data
      integer, intent(in) :: NT, NX, NY, NZ, mu, nu
      procedure(coordMuNuInterface), pointer :: coordMuNuFunc
      complex(kind=WC), dimension(NT) :: genTau
      ! working storage
      complex(kind=WC), dimension(NT, NX, NY, NZ, 3, 3) :: U_xd
      integer :: nnx, nny, nnz, nnt
      ! Calculate for each coordinate
      U_xd = coordMuNuFunc(data, NT, NX, NY, NZ, mu, nu)
      ! now take trace
      genTau = 0.0_WP
      do nnx = 1, nx
         do nny = 1, ny
            do nnz = 1, nz
               do nnt = 1, nt
                  genTau(nnt) = genTau(nnt) &
                                + U_xd(nnt, nnx, nnz, nnz, 1, 1) &
                                + U_xd(nnt, nnx, nnz, nnz, 2, 2) &
                                + U_xd(nnt, nnx, nnz, nnz, 3, 3)
               end do
            end do
         end do
      end do
      ! average over space
      genTau = genTau / real(nx * ny * nz, kind=WP)
   end function genSpaceAverageTr

   pure function genPlaquetteMuNu(data, NT, NX, NY, NZ, mu, nu) result(U_xd)
      complex(kind=WC), dimension(NT, NX, NY, NZ, 4, 3, 3), intent(in) :: data
      integer, intent(in) :: NT, NX, NY, NZ, mu, nu
      integer, dimension(4) :: coord
      complex(kind=WC), dimension(NT, NX, NY, NZ, 3, 3) :: U_xd
      complex(kind=WC), dimension(3, 3) :: clovLeaf
      integer, dimension(4) :: plaqPath
      integer, dimension(4) :: coordBase
      ! counters
      integer :: nnx, nny, nnz, nnt
      ! top right
      plaqPath(:) = (/mu, nu, -mu, -nu/)
#ifdef LOCALITYSUPPORT
      do concurrent(nnx=1:nx, nny=1:ny, nnz=1:nz, nnt=1:nt) &
         default(none) local_init(plaqPath) &
         local(coordBase, clovLeaf) shared(U_xd, data)
#else
         do nnx = 1, nx
            do nny = 1, ny
               do nnz = 1, nz
                  do nnt = 1, nt
#endif
                     ! Calculate the plaquette for this coordinate
                     coordBase = (/nnt, nnx, nny, nnz/)
                     ! 1
                     clovLeaf = genericPath(data, coordBase, plaqPath(:))
                     U_xd(nnt, nnx, nny, nnz, :, :) = clovLeaf
#ifdef LOCALITYSUPPORT
                  end do
#else
               end do
            end do
         end do
      end do
#endif
   end function genPlaquetteMuNu

end module FLUE_glueLoops
