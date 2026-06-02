module FLUE_SU2_random
   use FLUE_constants, only: WP, WC
   use FLUE_matrixConstants, only: Ident2x2, sigma1, sigma2, sigma3
   ! use stdlib_random, only: random_seed
   use stdlib_stats_distribution_uniform, only: rvs_uniform
   implicit none(type, external)
   private
   public :: constructSU2Matrix
   public :: randomNumbers
contains

   pure function constructSU2Matrix(r) result(U)
      real(kind=WP), dimension(0:3), intent(in) :: r
      complex(kind=WC), dimension(2, 2) :: U
      real(kind=WP), dimension(0:3) :: x
      real(kind=WP) :: normr
      ! U = x0 I + i xvec * sigmavec
      ! eqn 4.24
      !
      normr = SQRT(SUM(r**2))
      if (normr <= TINY(1.0_WP)) then
         x = [1.0_WP, 0.0_WP, 0.0_WP, 0.0_WP]
      else
         x = r / normr
      end if

      U = x(0) * Ident2x2
      U = U + CMPLX(0.0_WP, 1.0_WP) * x(1) * sigma1
      U = U + CMPLX(0.0_WP, 1.0_WP) * x(2) * sigma2
      U = U + CMPLX(0.0_WP, 1.0_WP) * x(3) * sigma3
   end function constructSU2Matrix

  !! get 4 random numbers uniformly distributed in (-0.5, 0.5)
   function randomNumbers() result(r)
      real(kind=WP), dimension(4) :: r
      r = rvs_uniform(loc=-0.5_WP, scale=1.0_WP, array_size=4)
   end function randomNumbers

end module FLUE_SU2_random
