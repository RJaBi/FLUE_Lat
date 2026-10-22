!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Read one OpenQCD gauge field and fix it to Coulomb gauge
!!
!! Usage:
!!   fpm run OQCD_coulombGauge -- gauge_file NT NS [tolerance] [max_iterations]
!!
!! Example:
!!   fpm run OQCD_coulombGauge -- config.openqcd 64 32 1.0e-12 5000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program OQCD_gaugeFix

   use FLUE, only: WP, WC,                         &
                   writeCompiler, writeGit,         &
                   readGaugeField_OpenQCD,          &
                   genPlaquette, gaugeFix,    &
                                 gaugeFix_options,   &
                                 gaugeFix_info, &
                                 applyFiniteGaugeTransformField

   implicit none(type, external)

   ! Command-line values
   character(len=512) :: input_file
   character(len=64)  :: argument

   ! Gauge field
   complex(kind=WC), dimension(:, :, :, :, :, :, :), allocatable :: gauge_field, UFixed

   complex(kind=WC), dimension(:,:,:,:,:,:), allocatable :: GTrans

   ! Coulomb-gauge controls and diagnostics
   type(gaugeFix_options) :: options
   type(gaugeFix_info)    :: info

   ! Lattice dimensions
   integer :: nt
   integer :: ns

   ! Plaquette measurements
   real(kind=WP) :: plaquette_before
   real(kind=WP) :: plaquette_after
   real(kind=WP) :: plaquette_difference
   real(kind=WP) :: plaquette_time
   integer :: number_plaquettes

   integer :: argument_status

   call writeCompiler()
   call writeGit()
   write (*, *)

   !---------------------------------------------------------------------------
   ! Parse command line
   !---------------------------------------------------------------------------

   if (command_argument_count() < 3) then
      call print_usage()
      error stop "Insufficient command-line arguments"
   end if

   call get_command_argument(1, input_file, status=argument_status)
   if (argument_status /= 0) then
      error stop "Could not read gauge-field filename"
   end if

   call get_integer_argument(2, nt)
   call get_integer_argument(3, ns)

   options = gaugeFix_options()

   if (command_argument_count() >= 4) then
      call get_real_argument(4, options%outerTolerance)
   end if

   if (command_argument_count() >= 5) then
      call get_integer_argument(5, options%maxOuterIterations)
   end if

   ! Inner CG controls
   !options% = 1.0e-10_WP
   !options%cg_atol = 1.0e-14_WP
   !options%max_cg_iterations = 1000

   ! Nonlinear update and line search
   !options%initial_step = 1.0_WP
   !options%minimum_step = 1.0e-8_WP
   !options%max_linesearch = 28

   !options%verbose = .true.
   !options%reunitarise = .false.

   options%fixType = 'coulomb'

   if (nt <= 0 .or. ns <= 0) then
      error stop "NT and NS must be positive"
   end if

   write (*, '(a)') "Gauge fixing"
   write (*, '(a,a)') "  Input file:         ", trim(input_file)
   write (*, '(a,i0)') "  NT:                 ", nt
   write (*, '(a,i0)') "  NS:                 ", ns
   write (*, '(a,es12.4)') "  Gauge tolerance:    ", options%outerTolerance
   write (*, '(a,es12.4)') "  Inner relative tol.:   ", options%innerRelTol
   write (*, '(a,i0)') "  Maximum outer iterations: ", options%maxOuterIterations
   write(*,*) " Gauge Fixing Type: "//trim(options%fixType)
   write (*, *)

   !---------------------------------------------------------------------------
   ! Read the gauge configuration
   !---------------------------------------------------------------------------

   allocate(gauge_field(3, 3, 4, nt, ns, ns, ns))
   allocate(GTrans(3,3,nt,ns,ns,ns))
   allocate(UFixed, mold=gauge_field)

   write (*, '(a)') "Reading OpenQCD gauge field..."

   gauge_field = readGaugeField_OpenQCD( &
      trim(input_file),                  &
      ns, ns, ns, nt,                    &
      fixSU3=.true.)

   write (*, '(a)') "Gauge field read successfully."
   write (*, *)

   !---------------------------------------------------------------------------
   ! Measurements before gauge fixing
   !---------------------------------------------------------------------------

   call genPlaquette( &
      gauge_field,                          &
      nt, ns, ns, ns,                      &
      1, 4, 4,                             &
      plaquette_before,                    &
      number_plaquettes,                   &
      plaquette_time)

   if (number_plaquettes <= 0) then
      error stop "genPlaquette returned no plaquettes"
   end if

   plaquette_before = plaquette_before &
      / real(number_plaquettes, kind=WP)

   write (*, '(a)') "Before gauge fixing"
   write (*, '(a,es24.16)') "  Average plaquette:   ", plaquette_before
   write (*, *)
   write(*,*) 'Before U(1,1,1,1,1,1,1)', gauge_field(1,1,1,1,1,1,1)

   ! Gauge fix
   call gaugeFix(gauge_field, options, info, GTrans)
   call applyFiniteGaugeTransformField(gauge_field, GTrans, UFixed)


   !---------------------------------------------------------------------------
   ! Measurements after gauge fixing
   !---------------------------------------------------------------------------

   call genPlaquette( &
      UFixed,                          &
      nt, ns, ns, ns,                      &
      1, 4, 4,                             &
      plaquette_after,                     &
      number_plaquettes,                   &
      plaquette_time)

   plaquette_after = plaquette_after &
      / real(number_plaquettes, kind=WP)


   plaquette_difference = abs(plaquette_after - plaquette_before)

   !---------------------------------------------------------------------------
   ! Report
   !---------------------------------------------------------------------------

   write (*, *)
   write (*, '(a)') repeat("=", 72)
   write (*, '(a)') "Coulomb gauge-fixing result"
   write (*, '(a)') repeat("=", 72)

   write (*, '(a,l1)')       "  Converged:              ", info%converged
   write (*, '(a,i0)')       "  Outer iterations:       ", info%OuterIterations
   write (*, '(a,i0)')       "  Inner iterations:    ", &
      info%InnerIterations

   write (*, '(a,es24.16)') "  Initial theta:          ", &
      info%initialTheta
   write (*, '(a,es24.16)') "  Final theta:            ", &
      info%finalTheta
   write (*, '(a,es24.16)') "  Maximum timeslice theta:", &
      info%maxFinalThetaPart

   write (*, '(a,es24.16)') "  Initial functional:     ", &
      info%initialFunctional
   write (*, '(a,es24.16)') "  Final functional:       ", &
      info%finalFunctional
   write (*, *)
   write (*, '(a,es24.16)') "  Plaquette before:       ", &
      plaquette_before
   write (*, '(a,es24.16)') "  Plaquette after:        ", &
      plaquette_after
   write (*, '(a,es24.16)') "  Absolute difference:    ", &
      plaquette_difference
   write(*,*) 'After UFixed(1,1,1,1,1,1,1)', UFixed(1,1,1,1,1,1,1)



   !---------------------------------------------------------------------------
   ! Sanity checks
   !---------------------------------------------------------------------------


   if (plaquette_difference > 1.0e3_WP * epsilon(1.0_WP)) then
      write (*, *)
      write (*, '(a)') &
         "WARNING: the plaquette changed by more than the roundoff target."
      write (*, '(a)') &
         "Check ExpIQ unitarity and the gauge-transformation indexing."
   end if

   if (.not. info%converged) then
      write (*, *)
      write (*, '(a)') &
         "WARNING: Coulomb gauge fixing did not reach the requested tolerance."
   end if

   deallocate(gauge_field)

contains

   !---------------------------------------------------------------------------
   ! Read an integer command-line argument.
   !---------------------------------------------------------------------------

   subroutine get_integer_argument(position, value)
      integer, intent(in) :: position
      integer, intent(out) :: value

      character(len=64) :: string
      integer :: io_status
      integer :: command_status

      call get_command_argument( &
         position, string, status=command_status)

      if (command_status /= 0) then
         write (*, '(a,i0)') &
            "Could not obtain command-line argument ", position
         error stop
      end if

      read (string, *, iostat=io_status) value

      if (io_status /= 0) then
         write (*, '(a,i0,a,a)') &
            "Argument ", position,                 &
            " is not a valid integer: ", trim(string)
         error stop
      end if

   end subroutine get_integer_argument


   !---------------------------------------------------------------------------
   ! Read a real command-line argument.
   !---------------------------------------------------------------------------

   subroutine get_real_argument(position, value)
      integer, intent(in) :: position
      real(kind=WP), intent(out) :: value

      character(len=64) :: string
      integer :: io_status
      integer :: command_status

      call get_command_argument( &
         position, string, status=command_status)

      if (command_status /= 0) then
         write (*, '(a,i0)') &
            "Could not obtain command-line argument ", position
         error stop
      end if

      read (string, *, iostat=io_status) value

      if (io_status /= 0) then
         write (*, '(a,i0,a,a)') &
            "Argument ", position,              &
            " is not a valid real: ", trim(string)
         error stop
      end if

   end subroutine get_real_argument


   !---------------------------------------------------------------------------
   ! Print command usage.
   !---------------------------------------------------------------------------

   subroutine print_usage()

      write (*, '(a)') &
         "Usage:"
      write (*, '(a)') &
         "  fpm run OQCD_coulombGauge -- gauge_file NT NS " // &
         "[tolerance] [max_iterations]"
      write (*, *)
      write (*, '(a)') &
         "Example:"
      write (*, '(a)') &
         "  fpm run OQCD_coulombGauge -- config.openqcd " // &
         "64 32 1.0e-12 5000"

   end subroutine print_usage

 end program OQCD_gaugeFix
