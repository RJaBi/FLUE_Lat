!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate the chromomagnetic portion of the field strength tensor
!! Ryan Bignell 2025
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program calcB
   use FLUE, only: WP, WC, writeCompiler, writeGit, &
                   ReadGaugeField_ILDG, ReadGaugeField_OpenQCD, genPlaquette, plaquette, &
                   magnetic, complement, jackknife_wp
   use tomlf, only: toml_table, toml_load, toml_array, get_value, toml_path
   use stdlib_io_npy, only: load_npy
   use csv_module, only: csv_file
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
   character(len=:), allocatable :: gaugePath, gaugeFormat, cfgListFile
   procedure(ReadGaugeInterface), pointer :: readGauge
   complex(kind=WC), dimension(:, :, :, :, :, :, :), allocatable :: U
   type(csv_file) :: csvf
   logical :: status_ok
   character(len=128), dimension(:), allocatable :: cfgList
   character(len=512) :: gaugeFile
   ! lattice dimensions
   integer :: NT, NS
   ! counters
   integer :: ii, nn, mu, icon
   ! plaq
   real(kind=WP) :: plaqFactor
   real(kind=WP) :: time, sumTrP
   integer :: nP
   ! jack variables
   real(kind=WP), dimension(:), allocatable :: B, aplaq, splaq, tplaq
   real(kind=WP), dimension(:), allocatable :: BJ, aplaqJ, splaqJ, tplaqJ
   real(kind=WP) :: BErr, aplaqErr, splaqErr, tplaqErr
   integer :: ncon

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
      call get_value(table, toml_path('fix', TRIM(fixLabels(ii)), 'cfgList'), cfglistFile)
      ! Get list of configurations
      write (*, *) 'cfgListFile is ', TRIM(cfgListFile)
      call csvf%read(TRIM(cfgListFile), status_ok=status_ok)
      call csvf%get(1, cfgList, status_ok)
      call csvf%destroy()
      ncon = SIZE(cfgList)
      ! config data
      allocate (B(ncon), aplaq(ncon), splaq(ncon), tplaq(ncon))
      ! jackknifes
      allocate (BJ(0:ncon), aplaqJ(0:ncon), splaqJ(0:ncon), tplaqJ(0:ncon))
      ! loop over them
      do icon = 1, ncon
         write (*, *) icon, TRIM(cfgList(icon))
         gaugeFile = TRIM(gaugePath)//'/'//TRIM(cfgList(icon))
         ! Now do stuff
         select case (gaugeFormat)
         case ('cssmILDG')
            ReadGauge => ReadGaugeField_ILDG
            plaqFactor = 1.0_WP
         case ('openqcd')
            readGauge => ReadGaugeField_OpenQCD
            plaqFactor = 1.0_WP
         case ('littlenpy')
            readGauge => npyFunctionWrapper
            plaqFactor = 1.0_WP
         case default
            write (*, *) 'gaugeFormat ', TRIM(gaugeFormat), ' not supported. Exiting'
            stop
         end select
         ! allocate space for the gaugefield
         allocate (U(NT, NS, NS, NS, 4, 3, 3))
         U = ReadGauge(TRIM(gaugeFile), NS, NS, NS, NT)
         call genPlaquette(U, NT, NS, NS, NS, 1, 4, 4, sumTrP, nP, time)
         aplaq(icon) = plaqFactor * sumTrp / real(nP, kind=WP)
         call genPlaquette(U, NT, NS, NS, NS, 1, 2, 4, sumTrP, nP, time)
         splaq(icon) = plaqFactor * sumTrp / real(nP, kind=WP)
         call genPlaquette(U, NT, NS, NS, NS, 1, 1, 4, sumTrP, nP, time)
         tplaq(icon) = plaqFactor * sumTrp / real(nP, kind=WP)
         !write (*, *) 'U  Plaquette for', ' is ', aplaq(icon), 'and took', time, 'seconds'
         !write (*, *) 'U SPlaquette for', ' is ', splaq(icon), 'and took', time, 'seconds'
         !write (*, *) 'U TPlaquette for', ' is ', tplaq(icon), 'and took', time, 'seconds'
         !write (*, *) 'magnetic is'
         B(icon) = magnetic(U, NT, NS, NS, NS)
         !write (*, *) 'magnetic', B(icon)
         deallocate (U)
      end do
      ! Take jackknifes
      ! Do mean first
      call Complement(BJ(0), B)
      call Complement(aplaqJ(0), aplaq)
      call Complement(splaqJ(0), splaq)
      call Complement(tplaqJ(0), tplaq)
      ! now 1st order jackknifes
      call Complement(ncon, BJ(1:), B)
      call Complement(ncon, aplaqJ(1:), aplaq)
      call Complement(ncon, splaqJ(1:), splaq)
      call Complement(ncon, tplaqJ(1:), tplaq)
      ! Now calculate uncertainties
      call Jackknife_wp(ncon, BJ, BErr)
      call Jackknife_wp(ncon, aplaqJ, aplaqErr)
      call Jackknife_wp(ncon, splaqJ, splaqErr)
      call Jackknife_wp(ncon, tplaqJ, tplaqErr)
      write (*, *) 'B is ', BJ(0), ' +- ', BErr
      write (*, *) 'aplaq is ', aplaqJ(0), ' +- ', aplaqErr
      write (*, *) 'splaq is ', splaqJ(0), ' +- ', splaqErr
      write (*, *) 'tplaq is ', tplaqJ(0), ' +- ', tplaqErr

      deallocate (B, aplaq, splaq, tplaq)
      deallocate (BJ, aplaqJ, splaqJ, tplaqJ)
      deallocate (cfgList)
   end do
contains

   function npyFunctionWrapper(filename, NX, NY, NZ, NT) result(G_x)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: NX, NY, NZ, NT
      complex(kind=WC), dimension(NT, NX, NY, NZ, 4, 3, 3) :: G_x
      complex(kind=WC), dimension(:, :, :, :, :, :, :), allocatable :: GNPY
      !counters
      integer :: xx, yy, zz, tt, aa, bb, mu
      integer, dimension(7) :: PYShape
      call load_npy(filename, GNPY)
      pyShape = SHAPE(GNPY)
      if (pyShape(1) < pyShape(7)) then
         ! For some reason, load_npy sometimes loads the array in opposite order..
         do concurrent(tt=1:NT, xx=1:NS, yy=1:NS, zz=1:NS, mu=1:4, aa=1:3, bb=1:3)
            G_x(tt, xx, yy, zz, mu, aa, bb) = GNPY(bb, aa, mu, zz, yy, xx, tt)
         end do
      else
         G_x = GNPY
      end if
   end function npyFunctionWrapper

end program calcB
