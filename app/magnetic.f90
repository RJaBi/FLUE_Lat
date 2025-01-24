!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the chromomagnetic portion of the field strength tensor
!! Ryan Bignell 2025
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program calcB
   use FLUE, only: WP, WC, writeCompiler, writeGit, &
                   ReadGaugeField_ILDG, ReadGaugeField_OpenQCD, genPlaquette, plaquette, &
                   magnetic
   use tomlf, only: toml_table, toml_load, toml_array, get_value, toml_path
   implicit none(external)
   ! IO vars
   character(len=128) :: tomlName
   character(len=:), allocatable :: strRead
   ! toml vars
   type(toml_table), allocatable :: table
   type(toml_array), pointer :: top_array
   ! What are we doing?
   integer :: nFixes, nTrans
   character(len=128), dimension(:), allocatable :: fixLabels
   character(len=:), allocatable :: gaugePath, gaugeFormat, transFormat
   character(len=512), dimension(:), allocatable :: transFiles
   procedure(ReadGaugeInterface), pointer :: readGauge
   complex(kind=WC), dimension(:, :, :, :, :, :, :), allocatable :: U, UT
   complex(kind=WC), dimension(:, :, :, :, :, :), allocatable :: G
   ! lattice dimensions
   integer :: NT, NS
   ! counters
   integer :: ii, nn, mu
   ! plaq
   real(kind=WP) :: plaqFactor, aplaq, splaq, tplaq
   real(kind=WP) :: time, sumTrP
   integer :: nP

   abstract interface
      function ReadGaugeInterface(filename, NX, NY, NZ, NT) result(U_xd)
         import :: WC
         implicit none(external)
         character(len=*), intent(in) :: filename
         integer, intent(in) :: NX, NY, NZ, NT
         complex(kind=WC), dimension(NT, NX, NY, NZ, 4, 3, 3) :: U_xd
      end function ReadGaugeInterface
   end interface

   call writeCompiler()
   call writeGit()
   write (*, *) ''

   if (COMMAND_ARGUMENT_COUNT() > 0) then
      call GET_COMMAND_ARGUMENT(1, tomlName)
      write (*, *) 'Reading toml file from ', TRIM(tomlName)
   else
      write (*, *) 'Pass the full path to the input toml on the command line'
      write (*, *) 'i.e. fpm run uzerobar -- mydir/input.toml'
      stop
   end if

   call toml_load(table, tomlName)
   call get_value(table, 'fixNum', nFixes)
   write (*, *) 'We will analyse ', nFixes, ' different cases'
   call get_value(table, 'fixLabels', top_array)
   allocate (fixLabels(nFixes))
   write (*, *) 'These are'
   do ii = 1, nFixes
      call get_value(top_array, ii, strRead)
      fixLabels(ii) = strRead
      write (*, *) TRIM(fixLabels(ii))
   end do
   ! Now let's do the load and do the analysis
   ! Rather than pre-load and do all analysis at once
   !allocate (aplaq(nFixes), splaq(nFixes), tplaq(nFixes))
   do ii = 1, nFixes
      write (*, *) ''
      write (*, *) ''
      write (*, *) ii, fixLabels(ii)
      ! Get singleton variables
      call get_value(table, toml_path('fix', TRIM(fixLabels(ii)), 'gaugePath'), gaugePath)
      call get_value(table, toml_path('fix', TRIM(fixLabels(ii)), 'gaugeFormat'), gaugeFormat)
      call get_value(table, toml_path('fix', TRIM(fixLabels(ii)), 'NT'), NT)
      call get_value(table, toml_path('fix', TRIM(fixLabels(ii)), 'NS'), NS)
      ! Now do stuff
      select case (gaugeFormat)
      case ('cssmILDG')
         ReadGauge => ReadGaugeField_ILDG
         plaqFactor = 9.0_WP
      case ('openqcd')
         readGauge => ReadGaugeField_OpenQCD
         plaqFactor = 1.0_WP
      case default
         write (*, *) 'gaugeFormat ', TRIM(gaugeFormat), ' not supported. Exiting'
         stop
      end select
      ! allocate space for the gaugefield
      allocate (U(NT, NS, NS, NS, 4, 3, 3))
      U = ReadGauge(gaugePath, NS, NS, NS, NT)

      write (*, *) 'For gauge format ', TRIM(gaugeFormat)
      write (*, *) 'Factor on plaquette is ', plaqFactor

      write (*, *) 'Using plaquette'
      call plaquette(U, 1, 4, 4, sumTrP, nP, time)
      !call plaquette(U, 1, 2, 2, sumTrP, nP, time)
      aplaq = plaqFactor * sumTrp / nP
      call plaquette(U, 1, 2, 4, sumTrP, nP, time)
      splaq = plaqFactor * sumTrp / nP
      call plaquette(U, 1, 1, 4, sumTrP, nP, time)
      tplaq = plaqFactor * sumTrp / nP
      write (*, *) 'U  Plaquette for', ' is ', aplaq, 'and took', time, 'seconds'
      write (*, *) 'U SPlaquette for', ' is ', splaq, 'and took', time, 'seconds'
      write (*, *) 'U TPlaquette for', ' is ', tplaq, 'and took', time, 'seconds'
      write (*, *) 'Using genPlaquette'
      call genPlaquette(U, NT, NS, NS, NS, 1, 4, 4, sumTrP, nP, time)
      !call genPlaquette(U, 1, 2, 2, sumTrP, nP, time)
      aplaq = plaqFactor * sumTrp / nP
      call genPlaquette(U, NT, NS, NS, NS, 1, 2, 4, sumTrP, nP, time)
      splaq = plaqFactor * sumTrp / nP
      call genPlaquette(U, NT, NS, NS, NS, 1, 1, 4, sumTrP, nP, time)
      tplaq = plaqFactor * sumTrp / nP
      write (*, *) 'U  Plaquette for', ' is ', aplaq, 'and took', time, 'seconds'
      write (*, *) 'U SPlaquette for', ' is ', splaq, 'and took', time, 'seconds'
      write (*, *) 'U TPlaquette for', ' is ', tplaq, 'and took', time, 'seconds'

      write (*, *) 'magnetic is ', magnetic(U, NT, NS, NS, NS)

      deallocate (U)
   end do
contains

end program calcB
