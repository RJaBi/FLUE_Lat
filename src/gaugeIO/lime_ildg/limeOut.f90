module FLUE_LIMEOut
  !< Functions to handle LIME files
  !< LIME = Lattice QCD Interchange Message Encapsulation
  !< LIME = Large Internet Message Encapsulation
  !< See the [SCIDAC documentation](https://usqcd-software.github.io/c-LIME/LIME_1p2.pdf)
  !< A LIME file consists of any number of concatenated LIME records
  !< Each record consists of a header followed by the data with trailing padding to an integer multiple of 8 bytes
  !< | header       | 144 bytes |
  !< | data         | up to $2^63$ bytes |
  !< | null padding | 0 to 7 bytes as needed |
  !<
  !< The header is documented below
  !< All is written in big-endian
  use, intrinsic :: iso_c_binding, only: c_null_char
  use, intrinsic :: iso_fortran_env, only: int64, int8
  use FLUE_constants, only: WP, WC
  use FLUE_endianIO, only: need_endian_swap, endian_swap_i64
  implicit none(type, external)
  private

  public :: writeLIMEHeader
  public :: writeLIMEPadding

contains
  pure function getControlWord(mb, me) result(word)
    !< Forms the int64 word to go in the header
    !< From the magic number, version and message end/begin bits
    integer, intent(in) :: mb
    !< 1 for message begin, 0 for continuinig
    integer, intent(in) :: me
    !< 1 for message end, 0 for continuing
    integer(int64) :: word
    !< output word of magic number, version, mb, me
    word = 0_int64
    ! Magic number
    word = ior(word, ishft(int(z'456789AB', int64), 32))
    ! Version = 1
    word = ior(word, ishft(1_int64, 16))
    ! MB/ME bits
    word = ior(word, ishft(int(iand(mb,1), int64), 15))
    word = ior(word, ishft(int(iand(me,1), int64), 14))
  end function getControlWord

  subroutine writeLIMEHeader(myUnit, mb, me, nbytes, recordType)
    !< This subroutine writes the LIME header to the open unit
    !< The Lime header is
    !< | bits    | contents |
    !< | 0 - 31  | LIME Magic Number |
    !< | 32 - 47 | LIME file version number |
    !< | 48 - 48 | LIME begin bit |
    !< | 49 - 49 | LIME end bit |
    !< | 50 - 63 | reserved |
    !< The LIME magic number 456789a (base16) is used to identify a record. This is a `long`
    !< The LIME file version number is a `short` integer
    !< The begin/end bit specifies where this record fits in the message
    !< i.e. 1 to begin, 1 to end
    !< So if you had 3 records, it would be
    !< 1, 0  for the first
    !< 0, 0  for the second
    !< 0, 1  for the third
    integer, intent(in) :: myUnit
    !< unit to write to
    integer, intent(in) :: mb
    !< mb=1 if start of message, 0 otherwise
    integer, intent(in) :: me
    !< me=1 if end of message, 0 otherwise
    integer(int64), intent(in) :: nbytes
    !< number of bytes of the payload
    character(*), intent(in) :: recordType
    !< a descriptive character string
    logical :: needEndianSwap
    !< LIME is big endian so figure out if we need to swap
    character(len=128) :: recordWrite
    !< the record which is written
    integer(int64) :: controlWord
    !< the control heaader: magic, mb, me, nbytes
    integer(int64) :: nBytesSwap
    !< nbytes, now in big endian
    needEndianSwap = need_endian_swap(.true.)
    if (len(trim(recordType)) > 127) then
       write(*,*) 'trim(recordType) = '//trim(recordType)
       error stop 'trim(recordType) too long'
    end if
    ! Fill with null char
    recordWrite = repeat(c_null_char,128)
    recordWrite(1:len(trim(recordType))) = recordType
    controlWord = getControlWord(mb, me)
    if (needEndianSwap) then
       nBytesSwap = endian_swap_i64(nBytes)
       controlWord = endian_swap_i64(controlWord)
    end if
    write(myUnit) controlWord
    write(myUnit) nbytesSwap
    write(myUnit) recordWrite
  end subroutine writeLIMEHeader

  subroutine writeLIMEPadding(myUnit, payLoadBytes)
    !< LIME records must end on multipe of 8 bytes
    !< Pad with zeros
    integer, intent(in) :: myUnit
    !< unit to write to
    integer(kind=int64), intent(in) :: payLoadBytes
    !< size of the data in the LIME record
    integer(kind=int8), dimension(7) :: padding
    !< max padding required
    integer :: paddingBytes
    !< calculated number of bytes required to be padded
    padding = 0_int8
    ! Calculate the padding required
    paddingBytes = modulo(8_int64 - modulo(payloadBytes, 8_int64), 8_int64)
    write(myUnit) padding(1:paddingBytes)
  end subroutine writeLIMEPadding

end module FLUE_LIMEOut
