module FLUE_ILDG_lime
  !< Functions to read and write ILDG (lime) format gaugefields
  use, intrinsic :: ISO_FORTRAN_ENV, only: OUTPUT_UNIT, int8, int32, int64
  use FLUE_constants, only: WP, WC
  use FLUE_ILDG_Bin, only: writeGaugeField_ILDG_Unit
  use FLUE_SU3MatrixOps, only: FixSU3Matrix
   use FLUE_limeOut, only: writeLimeHeader, writeLimePadding
   use foxy, only: xml_file, xml_tag
   use stdlib_strings, only: to_string
  implicit none(type, external)
  private
  public :: writeGaugeField_ILDG_lime

contains


  subroutine writeGaugeField_ILDG_lime(filename, U_xd, NX, NY, NZ, NT, lfn, xml, fixSU3)
    !< Writes the gaugefield to an ILDG compliant lime (xml) file
    !< xml are additionally messages/records
    character(*), intent(in) :: filename
    !< filename to write to
    integer, intent(in) :: NX, NY, NZ, NT
    !< Lattice dimensions
    complex(kind=WC), dimension(3, 3, 4, NT, NX, NY, NZ), intent(in) :: U_xd
    !< input gaugefield
    character(len=*), intent(in) :: lfn
    !< A string to go in message/record ildg-data-lfn
    type(xml_file), optional, dimension(:), intent(in) :: xml
    !< Optional additional XML to write
    !< Note this is not the ildg-required data
    logical, optional, intent(in) :: fixSU3
    !< optionall re-project to SU3
    integer :: myUnit
    !< where we write to
    type(xml_file) :: xmlWrite
    !< xml for required ILDG parts
    integer :: root
    !< updating parent for xml tag
    integer :: id
    !< updating id for xml tag
    integer(kind=int64) :: payloadBytes
    !< the payload bytes for the header
    integer :: ix
    !< counter
    logical :: fixSU3Set
    fixSU3Set = .false.
    if (present(fixSU3)) fixSU3Set = fixSU3

    open(newunit=myUnit, file=trim(filename), access='stream', &
         form='unformatted', status='replace', action='write')

    ! Write each of the additional records
    if (present(xml)) then
       do ix=1, size(xml)
          payloadBytes = len(xml(ix)%stringify(), kind=int64)
          ! First need to write the header
          call writeLimeHeader(myUnit, 1, 0, payloadBytes, 'ExtraXML_'//to_string(ix))
          ! Now write the xml
          write(myUnit) xml(ix)%stringify()
          call writeLimePadding(myUnit, payLoadBytes)
       end do
    end if

    ! Write the ildg-data-lfn message
    payloadBytes = len(trim(lfn), kind=int64)
    call writeLimeHeader(myUnit, 1, 1, payloadBytes, 'ildg-data-lfn')
    write(myUnit) trim(lfn)
    call writeLimePadding(myUnit, payLoadBytes)

    ! Assemble the ildgFormat record
    call xmlWrite%add_tag( &
         xml_tag(name='ildgFormat'), &
         id=root)
    call xmlWrite%add_tag(xml_tag(name='version'), parent=root, id=id)
    call xmlWrite%set_content(id, '1.2')
    call xmlWrite%add_tag(xml_tag(name='field'), parent=root, id=id)
    call xmlWrite%set_content(id, 'su3gauge')
    call xmlWrite%add_tag(xml_tag(name='precision'), parent=root, id=id)
    call xmlWrite%set_content(id, 64)
    call xmlWrite%add_tag(xml_tag(name='lx'), parent=root, id=id)
    call xmlWrite%set_content(id, NX)
    call xmlWrite%add_tag(xml_tag(name='ly'), parent=root, id=id)
    call xmlWrite%set_content(id, NY)
    call xmlWrite%add_tag(xml_tag(name='lz'), parent=root, id=id)
    call xmlWrite%set_content(id, NZ)
    call xmlWrite%add_tag(xml_tag(name='lt'), parent=root, id=id)
    call xmlWrite%set_content(id, NT)
    payloadBytes = len(xmlWrite%stringify(), kind=int64)
    write(*,*) payLoadBytes
    write(*,*) xmlWrite%stringify()
    ! First need to write the header
    call writeLimeHeader(myUnit, 1, 0, payloadBytes, 'ildg-format')
    ! Now write the xml
    write(myUnit) xmlWrite%stringify()
    call writeLimePadding(myUnit, payLoadBytes)
    ! Now need to write the 'ildg-binary-data' record
    payloadBytes = int( &
         storage_size(U_xd(1,1,1,1,1,1,1)) &
         * size(U_xd) / 8, &
         kind=int64)
    call writeLimeHeader(myUnit, 0, 1, payloadBytes, 'ildg-binary-data')
    call writeGaugeField_ILDG_Unit(myUnit, U_xd, NX, NY, NZ, NT, fixSU3Set)
    call writeLimePadding(myUnit, payLoadBytes)

    close(myUnit)

  end subroutine writeGaugeField_ILDG_lime

 end module FLUE_ILDG_lime
