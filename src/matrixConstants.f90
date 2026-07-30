module FLUE_matrixConstants
  !< Module: FLUE_matrixConstants
  !< Purpose:
  !<   This module provides a collection of predefined constant matrices used
  !<   throughout FLUE. These include identity matrices for SU(2) and SU(3)
  !<   operations together with the three Pauli matrices that form a basis of the
  !<   SU(2) Lie algebra.
  !<
  !< All matrices are defined as compile-time constants with complex-valued
  !< elements of kind ``WC``. Matrix entries are written in row-major order for
  !< human readability and converted to Fortran's native column-major storage
  !< layout through the ``RESHAPE`` intrinsic using ``order=[2,1]``.
  !<
  !<
  !< The Pauli matrices satisfy
  !< $
  !< [\sigma_i,\sigma_j] = 2 i \epsilon_{ijk}\sigma_k
  !< $
  !< and
  !< $
  !< \{\sigma_i,\sigma_j\} = 2\delta_{ij}I_2.
  !< $
   use FLUE_constants, only: WC, WP
   implicit none(type, external)
   complex(kind=WC), dimension(3, 3), parameter :: Ident3x3 = RESHAPE(source=[ &
                                                                      (1.0_WP, 0.0_WP), (0.0_WP, 0.0_WP), (0.0_WP, 0.0_WP), &
                                                                      (0.0_WP, 0.0_WP), (1.0_WP, 0.0_WP), (0.0_WP, 0.0_WP), &
                                                                      (0.0_WP, 0.0_WP), (0.0_WP, 0.0_WP), (1.0_WP, 0.0_WP)], &
                                                                      shape=[3, 3], order=[2, 1])
   !< 3x3 Identity matrix in complex variable
   complex(kind=WC), dimension(2, 2), parameter :: Ident2x2 = RESHAPE(source=[ &
                                                                      (1.0_WP, 0.0_WP), (0.0_WP, 0.0_WP), &
                                                                      (0.0_WP, 0.0_WP), (1.0_WP, 0.0_WP)], &
                                                                      shape=[2, 2], order=[2, 1])
   !< 2x2 Identity matrix in complex variable
   complex(kind=WC), dimension(2, 2), parameter :: sigma1 = RESHAPE(source=[ &
                                                                    (0.0_WP, 0.0_WP), (1.0_WP, 0.0_WP), &
                                                                    (1.0_WP, 0.0_WP), (0.0_WP, 0.0_WP)], &
                                                                    shape=[2, 2], order=[2, 1])
   !< First Pauli matrix (sigma_x)
   complex(kind=WC), dimension(2, 2), parameter :: sigma2 = RESHAPE(source=[ &
                                                                    (0.0_WP, 0.0_WP), (0.0_WP, -1.0_WP), &
                                                                    (0.0_WP, 1.0_WP), (0.0_WP, 0.0_WP)], &
                                                                    shape=[2, 2], order=[2, 1])
   !< Second Pauli matrix (sigma_y)
   complex(kind=WC), dimension(2, 2), parameter :: sigma3 = RESHAPE(source=[ &
                                                                    (1.0_WP, 0.0_WP), (0.0_WP, 0.0_WP), &
                                                                    (0.0_WP, 0.0_WP), (-1.0_WP, 0.0_WP)], &
                                                                    shape=[2, 2], order=[2, 1])
   !< Third Pauli matrix (sigma_z)
   private

   public :: Ident3x3, Ident2x2
   public :: sigma1, sigma2, sigma3

end module FLUE_matrixConstants
