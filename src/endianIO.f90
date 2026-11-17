module FLUE_endianIO
  !< A module to handle endian swapping
  use, intrinsic :: iso_fortran_env, only: int8, int16, int32, int64
  use, intrinsic :: iso_c_binding, only: c_int, c_long, c_double, c_double_complex
  use FLUE_constants, only: WP, WC
  implicit none(type, external)
  private
  public :: need_endian_swap


  interface endianSwap
     !< An interface to handle swapping different kinds
     !< Because kind numbers are sometimes reusable
     !< can not do i.e. real64 and c_long together in this interface
     module procedure endian_swap_c_int
     !< For integer c_int
     module procedure endian_swap_c_long
     !< for integer c_long
     module procedure endian_swap_wp
     !< for real WP
     module procedure endian_swap_wc
     !< for complex WC
  end interface endianSwap

   public :: endianSwap
   public :: endian_swap_i32
   public :: endian_swap_i64
   ! Explicit c_double are separate
   ! in case WP is not c_double, etc
   public :: endian_swap_c_double
   public :: endian_swap_c_double_complex

contains

  function host_is_little_endian() result(host_little)
    !< Determines if the host is little endian
    logical :: host_little
    !< true if host is little endian
    host_little = &
         1_int16 == transfer([1_int8, 0_int8], 0_int16)
  end function host_is_little_endian

  function host_is_big_endian() result(host_big)
    !< Determines if host is big endian
    logical :: host_big
    !< true if host is big endian
    host_big = .not. host_is_little_endian()
   end function host_is_big_endian

   function need_endian_swap(file_big_endian) result(need_swap)
     !< Determines if we need to swap endian before writing to the file
     !< i.e. if desired endianess of file is same as host
     logical, intent(in) :: file_big_endian
     !< true if file is big endian, false if little
     logical :: need_swap
     !< true if host and file are different endian, false otherwise
     need_swap = host_is_big_endian() .neqv. file_big_endian
   end function need_endian_swap

   elemental function endian_swap_i32(x) result(y)
     !< Swaps the endian of input int32 x
     !< does not do in place
     integer(int32), intent(in) :: x
     !< input x to be swapped
     integer(int32) :: y
     !< endian swapped output
     integer(int8), dimension(storage_size(x) / 8) :: bytes
     !< int8 used to reverse
     bytes = transfer(x, bytes)
     bytes = bytes(size(bytes):1:-1)
     y = transfer(bytes, y)
   end function endian_swap_i32

   elemental function endian_swap_i64(x) result(y)
     !< Swaps the endian of input int64 x
     !< does not do in place
     integer(int64), intent(in) :: x
     !< input x to be swapped
     integer(int64) :: y
     !< endian swapped output
     integer(int8), dimension(storage_size(x) / 8) :: bytes
     !< int8 used to reverse
     bytes = transfer(x, bytes)
     bytes = bytes(size(bytes):1:-1)
     y = transfer(bytes, y)
   end function endian_swap_i64

   elemental function endian_swap_c_int(x) result(y)
     !< Swaps the endian of input c_int x
     !< does not do in place
     integer(c_int), intent(in) :: x
     !< input x to be swapped
     integer(c_int) :: y
     !< endian swapped output
     integer(int8), dimension(storage_size(x) / 8) :: bytes
     !< int8 used to reverse
     bytes = transfer(x, bytes)
     bytes = bytes(size(bytes):1:-1)
     y = transfer(bytes, y)
   end function endian_swap_c_int

   elemental function endian_swap_c_long(x) result(y)
     !< Swaps the endian of input c_long x
     !< does not do in place
     integer(c_long), intent(in) :: x
     !< input x to be swapped
     integer(c_long) :: y
     !< endian swapped output
     integer(int8), dimension(storage_size(x) / 8) :: bytes
     !< int8 used to reverse
     bytes = transfer(x, bytes)
     bytes = bytes(size(bytes):1:-1)
     y = transfer(bytes, y)
   end function endian_swap_c_long

   elemental function endian_swap_wp(x) result(y)
     !< Swaps the endian of input real WP
     !< does not do in place
     real(WP), intent(in) :: x
     !< input x to be swapped
     real(WP) :: y
     !< endian swapped output
     integer(int8), dimension(storage_size(x) / 8) :: bytes
     !< int8 used to reverse
     bytes = transfer(x, bytes)
     bytes = bytes(size(bytes):1:-1)
     y = transfer(bytes, y)
   end function endian_swap_wp

   elemental function endian_swap_c_double(x) result(y)
     !< Swaps the endian of input real c_double
     !< does not do in place
     real(c_double), intent(in) :: x
     !< input x to be swapped
     real(c_double) :: y
     !< endian swapped output
     integer(int8), dimension(storage_size(x) / 8) :: bytes
     !< int8 used to reverse
     bytes = transfer(x, bytes)
     bytes = bytes(size(bytes):1:-1)
     y = transfer(bytes, y)
   end function endian_swap_c_double

   elemental function endian_swap_wc(z) result(w)
     !< Swaps the endian of input complex WC
     !< does not do in place
     complex(WC), intent(in) :: z
     !< input z to be swapped
     complex(WC) :: w
     !< endian swapped output
     w = cmplx( &
          endian_swap_wp(real(z, kind=WP)), &
          endian_swap_wp(aimag(z)), &
          kind=WC &
          )
   end function endian_swap_wc

   elemental function endian_swap_c_double_complex(z) result(w)
     !< Swaps the endian of input complex c_double_complex
     complex(c_double_complex), intent(in) :: z
     !< input z to be swapped
     complex(c_double_complex) :: w
     !< endian swapped output
     w = cmplx( &
          endian_swap_c_double(real(z, kind=c_double)), &
          endian_swap_c_double(aimag(z)), &
          kind=c_double_complex &
          )
   end function endian_swap_c_double_complex

 end module FLUE_endianIO
