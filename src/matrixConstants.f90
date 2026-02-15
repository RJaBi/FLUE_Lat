!! Matrix constants module
!!
!! Provides predefined constant matrices used throughout FLUE, including
!! the 3x3 identity matrix for SU(3) operations and the 2x2 Pauli/Sigma matrices
!! for SU(2) operations.
!!
!! Note: Matrix elements are stored in row-major order for visual clarity.
!! Fortran uses column-major order, so the optional order=[2,1] argument
!! is used in RESHAPE for correct representation.
!!
!! @author FLUE Developers
!! @version 1.0

module FLUE_matrixConstants
  use FLUE_constants, only: WC, WP
  implicit none(external)

  !! 3x3 Identity matrix for SU(3) operations
  !! Used in group theory calculations and matrix operations
  complex(kind=WC), dimension(3, 3), parameter :: Ident3x3 = RESHAPE(source=[ &
       (1.0_WP, 0.0_WP), (0.0_WP, 0.0_WP), (0.0_WP, 0.0_WP), &
       (0.0_WP, 0.0_WP), (1.0_WP, 0.0_WP), (0.0_WP, 0.0_WP), &
       (0.0_WP, 0.0_WP), (0.0_WP, 0.0_WP), (1.0_WP, 0.0_WP)], &
       shape=[3, 3], order=[2,1])

  !! Pauli/Sigma matrices for SU(2) operations
  !! First Pauli matrix (sigma_x)
  complex(kind=WC), dimension(2, 2), parameter :: sigma1 = RESHAPE(source=[ &
       (0.0_WP, 0.0_WP), (1.0_WP, 0.0_WP), &
       (1.0_WP, 0.0_WP), (0.0_WP, 0.0_WP)], &
       shape=[2, 2], order=[2,1])

  !! Second Pauli matrix (sigma_y)
  complex(kind=WC), dimension(2, 2), parameter :: sigma2 = RESHAPE(source=[ &
       (0.0_WP, 0.0_WP), (0.0_WP, -1.0_WP), &
       (0.0_WP, 1.0_WP), (0.0_WP, 0.0_WP)], &
       shape=[2, 2], order=[2,1])

  !! Third Pauli matrix (sigma_z)
  complex(kind=WC), dimension(2, 2), parameter :: sigma3 = RESHAPE(source=[ &
       (1.0_WP, 0.0_WP), (0.0_WP, 0.0_WP), &
       (0.0_WP, 0.0_WP), (-1.0_WP, 0.0_WP)], &
       shape=[2, 2], order=[2,1])

  private

  public :: Ident3x3
  public :: sigma1, sigma2, sigma3

end module FLUE_matrixConstants
