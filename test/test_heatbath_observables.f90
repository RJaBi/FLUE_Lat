module test_heatbath_observables
   use, intrinsic :: IEEE_ARITHMETIC, only: ieee_is_finite
   use FLUE_constants, only: WP, WC
   use FLUE_heatbath, only: updateLinks, build_colour_sites
   use FLUE_wloops, only: genPlaquette
   use philox, only: C64
   use stdlib_ascii, only: to_lower
   use stdlib_strings, only: to_string
   use test_helpers, only: seed_rng_fixed, fill_identity_su3
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use tomlf, only: toml_table, toml_error, toml_load, get_value
   implicit none(type, external)
   private
   public :: collect_heatbath_observables
contains
   subroutine collect_heatbath_observables(testsuite)
      type(unittest_type), allocatable, intent(OUT) :: testsuite(:)
      testsuite = [ &
                  new_unittest("heatbath_observables_wilson_xi1", test_heatbath_observables_wilson_xi1), &
                  new_unittest("heatbath_observables_wilson_xi1_beta1", test_heatbath_observables_wilson_xi1_beta1), &
                  new_unittest("heatbath_observables_wilson_xi10", test_heatbath_observables_wilson_xi10), &
                  new_unittest("heatbath_observables_symanzik_xi1", test_heatbath_observables_symanzik_xi1), &
                  new_unittest("heatbath_observables_symanzik_xi10", test_heatbath_observables_symanzik_xi10) &
                  ]
   end subroutine collect_heatbath_observables
   subroutine test_heatbath_observables_wilson_xi1(error)
      type(error_type), allocatable, intent(OUT) :: error
      call run_heatbath_observable_case(error, "wilson_xi1", "wilson", 1.0_WP, 6.0_WP)
   end subroutine test_heatbath_observables_wilson_xi1
   subroutine test_heatbath_observables_wilson_xi1_beta1(error)
      type(error_type), allocatable, intent(OUT) :: error
      call run_heatbath_observable_case(error, "wilson_xi1_beta1", "wilson", 1.0_WP, 1.0_WP)
   end subroutine test_heatbath_observables_wilson_xi1_beta1
   subroutine test_heatbath_observables_wilson_xi10(error)
      type(error_type), allocatable, intent(OUT) :: error
      call run_heatbath_observable_case(error, "wilson_xi10", "wilson", 10.0_WP, 6.0_WP)
   end subroutine test_heatbath_observables_wilson_xi10
   subroutine test_heatbath_observables_symanzik_xi1(error)
      type(error_type), allocatable, intent(OUT) :: error
      call run_heatbath_observable_case(error, "symanzik_xi1", "symanzik", 1.0_WP, 6.0_WP)
   end subroutine test_heatbath_observables_symanzik_xi1
   subroutine test_heatbath_observables_symanzik_xi10(error)
      type(error_type), allocatable, intent(OUT) :: error
      call run_heatbath_observable_case(error, "symanzik_xi10", "symanzik", 10.0_WP, 6.0_WP)
   end subroutine test_heatbath_observables_symanzik_xi10
   subroutine run_heatbath_observable_case(error, case_name, action_tag, xi, beta)
      type(error_type), allocatable, intent(OUT) :: error
      character(len=*), intent(IN) :: case_name, action_tag
      real(kind=WP), intent(IN) :: xi, beta
      integer, parameter :: NS = 4
      integer, parameter :: NT = 4
      integer, parameter :: nTraj = 20
      complex(kind=WC) :: U(3, 3, 4, NT, NS, NS, NS)
      complex(kind=WC) :: UNew(3, 3, 4, NT, NS, NS, NS)
      ! Sites linearisation
      integer, allocatable, dimension(:, :, :) :: sites_t, sites_x, sites_y, sites_z
      integer, allocatable, dimension(:, :) :: counts
      real(kind=WP) :: aplaq, splaq, tplaq
      real(kind=WP) :: ref_aplaq, ref_splaq, ref_tplaq
      real(kind=WP) :: tol_aplaq, tol_splaq, tol_tplaq
      real(kind=WP) :: sumTrP, time
      integer :: nPlaq, iTraj
      logical :: found_ref
      logical :: use_symanzik_sites
      integer(kind=C64), dimension(2), parameter :: key = (/1_C64, 2_C64/)
      use_symanzik_sites = .FALSE.
      use_symanzik_sites = (TRIM(to_lower(action_tag)) == 'symanzik')
      call fill_identity_su3(U)
      UNew = U
      call build_colour_sites(NT, NS, NS, NS, use_symanzik_sites, sites_t, sites_x, sites_y, sites_z, counts)
      do iTraj = 1, nTraj
         call updateLinks(U, beta, UNew, key, iTraj, sites_t, sites_x, sites_y, sites_z, counts, TRIM(action_tag), xi=xi)
         U = UNew
      end do
      call genPlaquette(U, NT, NS, NS, NS, 1, 4, 4, sumTrP, nPlaq, time)
      aplaq = sumTrP / (3.0_WP * real(nPlaq, kind=WP))
      if (xi /= 1.0_WP) then
         call genPlaquette(U, NT, NS, NS, NS, 2, 4, 4, sumTrP, nPlaq, time)
         splaq = sumTrP / (3.0_WP * real(nPlaq, kind=WP))
         call genPlaquette(U, NT, NS, NS, NS, 1, 1, 4, sumTrP, nPlaq, time)
         tplaq = sumTrP / (3.0_WP * real(nPlaq, kind=WP))
      else
         splaq = 0.0_WP
         tplaq = 0.0_WP
      end if
      call read_heatbath_observable_reference( &
         "testdata/reference_values.toml", case_name, xi /= 1.0_WP, &
         ref_aplaq, tol_aplaq, ref_splaq, tol_splaq, ref_tplaq, tol_tplaq, found_ref)
      call check(error, found_ref, &
                 "heatbath observable reference [observables.heatbath."//TRIM(case_name)//"] must exist")
      if (ALLOCATED(error)) return
      ! Basic sanity first
      call check(error, ieee_is_finite(aplaq) .AND. ABS(aplaq) <= 1.0_WP + 1.0E-12_WP, &
                 "average plaquette "//to_string(aplaq)//" should be finite and bounded")
      if (ALLOCATED(error)) return
      if (xi /= 1.0_WP) then
         call check(error, ieee_is_finite(splaq) .AND. ABS(splaq) <= 1.0_WP + 1.0E-12_WP, &
                    "spatial plaquette "//to_string(splaq)//"  should be finite and bounded")
         if (ALLOCATED(error)) return
         call check(error, ieee_is_finite(tplaq) .AND. ABS(tplaq) <= 1.0_WP + 1.0E-12_WP, &
                    "temporal plaquette "//to_string(tplaq)//" should be finite and bounded")
         if (ALLOCATED(error)) return
      end if
      ! Reference-value regression
      call check(error, ABS(aplaq - ref_aplaq) < tol_aplaq, &
                 "average plaquette "//to_string(aplaq)//" should match the stored heatbath reference "//to_string(ref_aplaq))
      if (ALLOCATED(error)) return
      if (xi /= 1.0_WP) then
         call check(error, ABS(splaq - ref_splaq) < tol_splaq, &
              "spatial plaquette "//to_string(splaq)//" should match the stored heatbath reference "//to_string(ref_splaq))
         if (ALLOCATED(error)) return
         call check(error, ABS(tplaq - ref_tplaq) < tol_tplaq, &
              "temporal plaquette "//to_string(tplaq)//" should match the stored heatbath reference "//to_string(ref_tplaq))
      end if
   end subroutine run_heatbath_observable_case

   subroutine read_heatbath_observable_reference(filename, case_name, anisotropic, &
                                                 ref_aplaq, tol_aplaq, &
                                                 ref_splaq, tol_splaq, &
                                                 ref_tplaq, tol_tplaq, found)
      character(len=*), intent(IN) :: filename
      character(len=*), intent(IN) :: case_name
      logical, intent(IN) :: anisotropic
      real(kind=WP), intent(OUT) :: ref_aplaq, tol_aplaq
      real(kind=WP), intent(OUT) :: ref_splaq, tol_splaq
      real(kind=WP), intent(OUT) :: ref_tplaq, tol_tplaq
      logical, intent(OUT) :: found

      type(toml_table), allocatable :: root
      type(toml_table), pointer :: obs
      type(toml_table), pointer :: hb
      type(toml_table), pointer :: case_tbl
      type(toml_error), allocatable :: err
      integer :: stat
      ref_aplaq = -HUGE(1.0_WP)
      tol_aplaq = -1.0_WP
      ref_splaq = -HUGE(1.0_WP)
      tol_splaq = -1.0_WP
      ref_tplaq = -HUGE(1.0_WP)
      tol_tplaq = -1.0_WP
      found = .FALSE.
      call toml_load(root, TRIM(filename), error=err)
      if (ALLOCATED(err)) return
      call get_value(root, "observables", obs, stat=stat)
      if (stat /= 0 .OR. .NOT. ASSOCIATED(obs)) return
      call get_value(obs, "heatbath", hb, stat=stat)
      if (stat /= 0 .OR. .NOT. ASSOCIATED(hb)) return
      call get_value(hb, TRIM(case_name), case_tbl, stat=stat)
      if (stat /= 0 .OR. .NOT. ASSOCIATED(case_tbl)) return
      call get_value(case_tbl, "aplaq", ref_aplaq, stat=stat)
      if (stat /= 0) return
      call get_value(case_tbl, "aplaq_tol", tol_aplaq, stat=stat)
      if (stat /= 0) return
      if (anisotropic) then
         call get_value(case_tbl, "splaq", ref_splaq, stat=stat)
         if (stat /= 0) return
         call get_value(case_tbl, "splaq_tol", tol_splaq, stat=stat)
         if (stat /= 0) return
         call get_value(case_tbl, "tplaq", ref_tplaq, stat=stat)
         if (stat /= 0) return
         call get_value(case_tbl, "tplaq_tol", tol_tplaq, stat=stat)
         if (stat /= 0) return
      end if
      found = .TRUE.
   end subroutine read_heatbath_observable_reference

end module test_heatbath_observables
