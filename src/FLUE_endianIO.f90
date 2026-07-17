module FLUE_endianIO
   use, intrinsic :: iso_fortran_env, only: int8, int16, int32
   use, intrinsic :: iso_c_binding, only: c_int, c_long, c_double, c_double_complex
   use FLUE_constants, only: WP, WC
   implicit none(type, external)
   private

   !public :: host_is_little_endian
   !public :: host_is_big_endian
   public :: need_endian_swap


   interface endianSwap
      ! FLUE_constants does not define a type for ints
      ! so can keep them here
      module procedure endian_swap_c_int
      module procedure endian_swap_c_long
      ! It does for WP, WC
      ! so put them here
      module procedure endian_swap_wp
      module procedure endian_swap_wc
   end interface endianSwap

   public :: endianSwap
   public :: endian_swap_i32
   !public :: endian_swap_c_int
   !public :: endian_swap_c_long
   !public :: endian_swap_wp
   !public :: endian_swap_wc
   ! Explicit c_double are separate
   ! in case WP is not c_double, etc
   public :: endian_swap_c_double
   public :: endian_swap_c_double_complex

contains

   logical function host_is_little_endian()
      host_is_little_endian = &
         1_int16 == transfer([1_int8, 0_int8], 0_int16)
   end function host_is_little_endian


   logical function host_is_big_endian()
      host_is_big_endian = .not. host_is_little_endian()
   end function host_is_big_endian


   logical function need_endian_swap(file_big_endian)
      logical, intent(in) :: file_big_endian

      need_endian_swap = host_is_big_endian() .neqv. file_big_endian
   end function need_endian_swap


   elemental function endian_swap_i32(x) result(y)
      integer(int32), intent(in) :: x
      integer(int32) :: y
      integer(int8) :: bytes(storage_size(x) / 8)

      bytes = transfer(x, bytes)
      bytes = bytes(size(bytes):1:-1)
      y = transfer(bytes, y)
   end function endian_swap_i32


   elemental function endian_swap_c_int(x) result(y)
      integer(c_int), intent(in) :: x
      integer(c_int) :: y
      integer(int8) :: bytes(storage_size(x) / 8)

      bytes = transfer(x, bytes)
      bytes = bytes(size(bytes):1:-1)
      y = transfer(bytes, y)
   end function endian_swap_c_int


   elemental function endian_swap_c_long(x) result(y)
      integer(c_long), intent(in) :: x
      integer(c_long) :: y
      integer(int8) :: bytes(storage_size(x) / 8)

      bytes = transfer(x, bytes)
      bytes = bytes(size(bytes):1:-1)
      y = transfer(bytes, y)
   end function endian_swap_c_long


   elemental function endian_swap_wp(x) result(y)
      real(WP), intent(in) :: x
      real(WP) :: y
      integer(int8) :: bytes(storage_size(x) / 8)

      bytes = transfer(x, bytes)
      bytes = bytes(size(bytes):1:-1)
      y = transfer(bytes, y)
   end function endian_swap_wp


   elemental function endian_swap_c_double(x) result(y)
      real(c_double), intent(in) :: x
      real(c_double) :: y
      integer(int8) :: bytes(storage_size(x) / 8)

      bytes = transfer(x, bytes)
      bytes = bytes(size(bytes):1:-1)
      y = transfer(bytes, y)
   end function endian_swap_c_double


   elemental function endian_swap_wc(z) result(w)
      complex(WC), intent(in) :: z
      complex(WC) :: w

      w = cmplx( &
         endian_swap_wp(real(z, kind=WP)), &
         endian_swap_wp(aimag(z)), &
         kind=WC &
      )
   end function endian_swap_wc


   elemental function endian_swap_c_double_complex(z) result(w)
      complex(c_double_complex), intent(in) :: z
      complex(c_double_complex) :: w

      w = cmplx( &
         endian_swap_c_double(real(z, kind=c_double)), &
         endian_swap_c_double(aimag(z)), &
         kind=c_double_complex &
      )
   end function endian_swap_c_double_complex

 end module FLUE_endianIO
