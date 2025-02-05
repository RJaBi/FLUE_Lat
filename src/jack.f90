module FLUE_jack

   use FLUE_constants, only: WP
   implicit none(external)
   private
   public :: Complement
   public :: jackknife_wp

   interface Complement
      module procedure Complement0
      module procedure Complement1
      module procedure Complement2
      module procedure Complement3
   end interface Complement

contains

   pure subroutine Complement0(cset, data)
      real(kind=WP), intent(out) :: cset
      real(kind=WP), dimension(:), intent(in) :: data
      cset = SUM(data) / real(SIZE(data), kind=WP)
   end subroutine Complement0

   pure subroutine Complement1(ncon, cset, data)
      integer, intent(in) :: ncon
      real(kind=WP), dimension(ncon), intent(out) :: cset
      real(kind=WP), dimension(ncon), intent(in) :: data
      cset = (SUM(data) - data) / real(SIZE(data) - 1, kind=WP)
   end subroutine Complement1

   pure subroutine Complement2(ncon, cset, data, icon)
      integer, intent(in) :: ncon
      real(kind=WP), dimension(ncon), intent(out) :: cset
      real(kind=WP), dimension(ncon), intent(in) :: data
      integer, intent(in) :: icon
      cset = (SUM(data) - data(icon) - data) / real(SIZE(data) - 2, kind=WP)
   end subroutine Complement2

   pure subroutine Complement3(ncon, cset, data, icon, jcon)
      integer, intent(in) :: ncon
      real(kind=WP), dimension(ncon), intent(out) :: cset
      real(kind=WP), dimension(ncon), intent(in) :: data
      integer, intent(in) :: icon
      integer, intent(in) :: jcon
      cset = (SUM(data) - data(icon) - data(jcon) - data) / real(SIZE(data) - 3, kind=WP)
   end subroutine Complement3

   subroutine JackKnife_wp(ncon, c, err, bias)
      ! This subroutine works
      ! but just rescale np.cov by (ncon - 1)
      integer, intent(in) :: ncon
      real(kind=WP), dimension(0:ncon), intent(in) :: c
      real(kind=WP), intent(out) :: err
      real(kind=WP), optional, intent(out) :: bias
      real(kind=WP) :: avg
      avg = SUM(c(1:)) / real(SIZE(c) - 1, kind=WP)
      err = SQRT(SUM((c(1:) - avg)**2) * real(SIZE(c) - 2, kind=WP) / real(SIZE(c) - 1, kind=WP))
      if (PRESENT(bias)) bias = c(0) - avg
   end subroutine JackKnife_wp

end module FLUE_jack
