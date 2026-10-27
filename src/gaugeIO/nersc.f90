module FLUE_NERSCFileIO
  !< Subroutines to read and write NERSC format gaugefields
  use FLUE_constants, only : WP, WC
  use FLUE_SU3MatrixOps, only: FixSU3Matrix
  use FLUE_wloops, only : genPlaquette
  implicit none(type, external)
  private
  public :: WriteGaugeField_NERSC
  public :: ReadGaugeField_NERSC
  public :: header_NERSC
  public :: writeGaugeField_NERSC_FullHeader

  type :: header_NERSC
     character(len=3) :: storage_format = '1.0'
     !< I don't know what other versions are...
     integer :: dimension_1, dimension_2, dimension_3, dimension_4
     !< Lattice dimensions, nx, ny, nz, nt
     character(len=12) :: floating_point = 'IEEE64LITTLE'
     !< Technically this should be IEEE64 in whatever compiled endian & precision (WC)
     character(len=16) :: datatype = '4D_SU3_GAUGE_3X3'
     !< Only support 4D_SU3_GAUGE_3X3
     character(len=7) :: checksum = '1234567'
     !< Proper checksum not implemented yet
     real(kind=WP) :: plaquette, link_trace
     !< Average plaquette and link_trace
     integer :: sequence_number = -1
     !< Probably the config number?. negative numbers wont be written
     character(len=256) :: creator = ''
     !< The creator of the nersc file
     character(len=256) :: creator_hardware = ''
     !< what hardware was used
     character(len=256) :: creation_date = ''
     !< When the gaugefield was generated in format for C++ std::get_time(&tm, "%A %B %d %H:%M:%S %Y")
     character(len=256) :: archive_date = ''
     !< When it was made into NERSC format in format for C++ std::get_time(&tm, "%A %B %d %H:%M:%S %Y")
     ! "%A %B %d %H:%M:%S %Y" Means
     ! Full Weekday | Full Month | Day of the month 2 digits |
     ! hour range [0, 23], leading zero allowed but not required |
     ! minute range [0, 59], leading zero allowed but not required |
     ! second range [00, 59], leading zero allowed but not required |
     ! year as 4 digit number, leading zero allowed but not required
     ! i.e. an example is
     ! Tuesday October 27 14:48:30 2026
     character(len=256) :: ensemble_id = ''
     !< some form of ID for the ensemble?
     character(len=256) :: ensemble_label = ''
     !< some form of label for the ensemble
  end type header_NERSC


contains

  subroutine WriteGaugeField_NERSC(filename, U_in, NX, NY, NZ, NT)
    !< Writes a gaugefield to file in NESRC format
    character(len=*), intent(in) :: filename
    !< File to write to
    integer, intent(in) :: NX, NY, NZ, NT
    !< Lattice dimensions
    complex(kind=WC), dimension(3,3,4,NT,NX,NY,NZ), intent(in) :: U_in

    type(header_NERSC) :: header
    call writeGaugeField_NERSC_FullHeader(filename, U_in, NX, NY, NZ, NT, header)
  end subroutine WriteGaugeField_NERSC

  subroutine WriteGaugeField_NERSC_FullHeader(filename, U_in, NX, NY, NZ, NT,header)
    !< Writes a gaugefield to file in NESRC format
    character(len=*), intent(in) :: filename
    !< File to write to
    integer, intent(in) :: NX, NY, NZ, NT
    !< Lattice dimensions
    complex(kind=WC), dimension(3,3,4,NT,NX,NY,NZ), intent(in) :: U_in
    !< gaugefield to write
    type(header_NERSC), intent(in) :: header
    !< Header to write
    complex(kind=WC), dimension(3,3,4,NT,NX,NY,NZ) :: UShift
    !< gaugefield but cshift mu so x is first
    integer :: fileUnit
    !< file units
    real(kind=WP) :: plaq
    !< plaquette value
    real(kind=WP) :: sumTrP
    !< sum of all plaquettes
    real(kind=WP) :: timing
    !< time taken to calculate plaquette
    integer :: NP
    !< number of plaquette
    real(kind=WP) :: link_trace, singleTrace
    !< mean of real trace (colour) of al links
    integer :: it, ix, iy, iz, mu, pos
    !< counters
    ! Calculate plaquette to go in header
    call genPlaquette(U_in, NT, NX, NY, NZ, 1, 4, 4, &
         sumTrP, NP, timing)
    plaq = sumTrP / real(NP, kind=WP)
    plaq = plaq / 3.0_WP  ! Extra factor to match Jonivar
    !< Calculate trace of links to go in header
    link_trace = 0.0_WP
    do concurrent(it=1:nt, ix=1:nx, iy=1:ny, iz=1:nz, mu=1:4) &
         default(none) reduce(+:link_trace) shared(U_in) local(singleTrace)
       singleTrace = real( 0.0_WP &
            + U_in(1,1,mu,it,ix,iy,iz) &
            + U_in(2,2,mu,it,ix,iy,iz) &
            + U_in(3,3,mu,it,ix,iy,iz) &
            )
       link_trace = link_trace + singleTrace
    end do
    link_trace = link_trace / &
         real(4*NT*NX*NY*NZ, kind=WP)
    link_trace = link_trace / 3.0_WP  ! Extra factor to match Jonivar
    ! Write header in ascii, i.e. 'a'
    open(newunit=fileUnit, &
         file=trim(filename), &
         access='stream', &
         form='formatted', &
         status='replace', &
         action='write')

    write(fileUnit,'(a)') 'BEGIN_HEADER'
    write(fileUnit,'(a)') 'STORAGE_FORMAT = '// trim(header%storage_format)
    write(fileUnit,'(a,i0)') 'DIMENSION_1 = ', NX
    write(fileUnit,'(a,i0)') 'DIMENSION_2 = ', NY
    write(fileUnit,'(a,i0)') 'DIMENSION_3 = ', NZ
    write(fileUnit,'(a,i0)') 'DIMENSION_4 = ', NT
    write(fileUnit,'(a)') 'FLOATING_POINT = '// trim(header%floating_point)
    write(fileUnit,'(a)') 'DATATYPE = '// trim(header%datatype)
    write(fileUnit,'(a)') 'CHECKSUM = '// trim(header%checksum)
    write(fileUnit,'(a,eS24.16e3)') 'PLAQUETTE = ', plaq
    write(fileUnit,'(a,eS24.16e3)') 'LINK_TRACE = ', link_trace
    ! Ones below this are not in Jonivar's convertor
    if (header%sequence_number .ge. 0) then
       write(fileUnit,'(a,i0)') 'SEQUENCE_NUMBER = ', header%sequence_number
    end if
    if (len(trim(header%creator)) .gt. 0 ) then
       write(fileUnit, '(a)') 'CREATOR = '// trim(header%creator)
    end if
    if (len(trim(header%creator_hardware)) .gt. 0 ) then
       write(fileUnit, '(a)') 'CREATOR_HARDWARE = '// trim(header%creator_hardware)
    end if
    if (len(trim(header%creation_date)) .gt. 0 ) then
       write(fileUnit, '(a)') 'CREATION_DATE = '// trim(header%creation_date)
    end if
    if (len(trim(header%archive_date)) .gt. 0 ) then
       write(fileUnit, '(a)') 'ARCHIVE_DATE = '// trim(header%archive_date)
    end if
    write(fileUnit,'(a)') 'END_HEADER'
    close(fileUnit)
    ! cshift to mu is xyzt
    UShift = cshift(U_in, 1, dim=3)
    ! Re-open in append mode and write binary gaugefield
    open(newunit=fileUnit, &
         file=trim(filename), &
         access='stream', &
         form='unformatted', &
         status='old', &
         action='write', &
         position='append', &
         convert='little_endian')
    !
    ! NERSC ordering:
    ! x fastest, then y, z, t
    ! write four directions for every site
    !
    do it = 1, NT
       do iz = 1, NZ
          do iy = 1, NY
             do ix = 1, NX
                do mu = 1, 4
                   write(fileUnit) transpose(UShift(:,:,mu,it,ix,iy,iz))
                end do
             end do
          end do
       end do
    end do

    close(fileUnit)

  end subroutine WriteGaugeField_NERSC_FullHeader

  function ReadGaugeField_NERSC(filename, NX, NY, NZ, NT, fixSU3) result(U_out)
    !< Reads a NERSC format gaugefield into internal form [3,3,4,NT,NX,NY,NZ]
    character(len=*), intent(in) :: filename
    !< File to be read
    integer, intent(in) :: NX, NY, NZ, NT
    !< Lattice dimensions
    logical, optional, intent(IN) :: fixSU3
    !< Optionally reproject each link to SU3
    complex(kind=WC), dimension(3,3,4,NT,NX,NY,NZ) :: U_out
    !< Output gaugefield
    integer :: fileUnit
    !< File unit
    integer :: it, ix, iy, iz, mu, pos
    !< counters
    logical :: fixSU3Set
    !< whether we are doing SU3 re-projection
    character(len=256) :: stringRead
    !< To read each line of the header into a string
    integer, parameter :: recordLength = 8*2*3*3
    !< double precsion complex * SU3 link
    complex(kind=WC), dimension(3,3) :: ULink
    !< Holds single SU3 link for reading
    if (PRESENT(fixSU3)) then
       fixSU3Set = fixSU3
    else
       fixSU3Set = .TRUE.
    end if

    ! read header in ascii, i.e. 'a'
    open(newunit=fileUnit, &
         file=trim(filename), &
         access='stream', &
         form='formatted', &
         status='replace', &
         action='read')

    headerLoop: do
       read(fileUnit,'(a)') stringRead  ! 'BEGIN_HEADER'
       if (trim(stringRead) == 'END_HEADER') then
          exit headerLoop
       end if
    end do headerLoop
    ! write(*,*) trim(stringRead)
    !read(fileUnit,'(a)') stringRead  ! 'STORAGE_FORMAT = 1.0'
    ! write(*,*) trim(stringRead)
    !read(fileUnit,'(a)') stringRead  ! 'DIMENSION_1 = ', NX
    ! write(*,*) trim(stringRead)
    !read(fileUnit,'(a)') stringRead  ! 'DIMENSION_2 = ', NY
    ! write(*,*) trim(stringRead)
    !read(fileUnit,'(a)') stringRead  ! 'DIMENSION_3 = ', NZ
    ! write(*,*) trim(stringRead)
    !read(fileUnit,'(a)') stringRead  ! 'DIMENSION_4 = ', NT
    ! write(*,*) trim(stringRead)
    !read(fileUnit,'(a)') stringRead  ! 'FLOATING_POINT = IEEE64LITTLE'
    ! write(*,*) trim(stringRead)
    !read(fileUnit,'(a)') stringRead  ! 'DATATYPE = 4D_SU3_GAUGE_3X3'
    ! write(*,*) trim(stringRead)
    !read(fileUnit,'(a)') stringRead  ! 'CHECKSUM = 1234567'
    ! write(*,*) trim(stringRead)
    !read(fileUnit,'(a)') stringRead  ! 'PLAQUETTE = ', plaq
    ! write(*,*) trim(stringRead)
    !read(fileUnit,'(a)') stringRead  ! 'LINK_TRACE = ', link_trace
    ! write(*,*) trim(stringRead)
    !read(fileUnit,'(a)') stringRead  ! 'END_HEADER'
    ! write(*,*) trim(stringRead)
    inquire(fileUnit, pos=pos)
    close(fileUnit)
    ! write(*,*) 'pos after header is ', pos
    ! Re-open in append mode and write binary gaugefield
    open(newunit=fileUnit, &
         file=trim(filename), &
         access='stream', &
         form='unformatted', &
         status='old', &
         action='read')
    !
    ! NERSC ordering:
    ! x fastest, then y, z, t
    ! write four directions for every site (xyzt)
    !
    do it = 1, NT
       do iz = 1, NZ
          do iy = 1, NY
             do ix = 1, NX
                do mu = 1, 4
                   read(fileUnit, pos=pos) ULink
                   if (fixSU3Set) then
                      call FixSU3Matrix(ULink)
                   end if
                   pos = pos + recordLength
                   U_out(:,:,mu,it,ix,iy,iz) = transpose(ULink)
                end do
             end do
          end do
       end do
    end do
    ! cshift to match txyz mu ordering
    U_out = cshift(U_out, -1, dim=3)

    close(fileUnit)

  end function ReadGaugeField_NERSC

end module FLUE_NERSCFileIO
