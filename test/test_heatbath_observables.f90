MODULE test_heatbath_observables
  USE FLUE_constants, ONLY : WP, WC
   USE FLUE_heatbath, ONLY : updateLinks
   USE FLUE_wloops, ONLY : genPlaquette
   USE, INTRINSIC :: ieee_arithmetic, ONLY : ieee_is_finite
   USE test_helpers, ONLY : seed_rng_fixed, fill_identity_su3
   USE testdrive, ONLY : new_unittest, unittest_type, error_type, check
   USE tomlf, ONLY : toml_table, toml_error, toml_load, get_value
   IMPLICIT NONE(TYPE, EXTERNAL)
   PRIVATE
   PUBLIC :: collect_heatbath_observables
CONTAINS
   SUBROUTINE collect_heatbath_observables(testsuite)
      TYPE(unittest_type), ALLOCATABLE, INTENT(OUT) :: testsuite(:)
      testsuite = [ &
         new_unittest("heatbath_observables_wilson_xi1",    test_heatbath_observables_wilson_xi1), &
         new_unittest("heatbath_observables_wilson_xi10",   test_heatbath_observables_wilson_xi10), &
         new_unittest("heatbath_observables_symanzik_xi1",  test_heatbath_observables_symanzik_xi1), &
         new_unittest("heatbath_observables_symanzik_xi10", test_heatbath_observables_symanzik_xi10) &
      ]
   END SUBROUTINE collect_heatbath_observables
   SUBROUTINE test_heatbath_observables_wilson_xi1(error)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      CALL run_heatbath_observable_case(error, "wilson_xi1", "wilson", 1.0_WP, 6.0_WP)
   END SUBROUTINE test_heatbath_observables_wilson_xi1
   SUBROUTINE test_heatbath_observables_wilson_xi10(error)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      CALL run_heatbath_observable_case(error, "wilson_xi10", "wilson", 10.0_WP, 6.8_WP)
   END SUBROUTINE test_heatbath_observables_wilson_xi10
   SUBROUTINE test_heatbath_observables_symanzik_xi1(error)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      CALL run_heatbath_observable_case(error, "symanzik_xi1", "symanzik", 1.0_WP, 6.0_WP)
   END SUBROUTINE test_heatbath_observables_symanzik_xi1
   SUBROUTINE test_heatbath_observables_symanzik_xi10(error)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      CALL run_heatbath_observable_case(error, "symanzik_xi10", "symanzik", 10.0_WP, 6.8_WP)
   END SUBROUTINE test_heatbath_observables_symanzik_xi10
   SUBROUTINE run_heatbath_observable_case(error, case_name, action_tag, xi, beta)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      CHARACTER(len=*), INTENT(IN) :: case_name, action_tag
      REAL(WP), INTENT(IN) :: xi, beta
      INTEGER, PARAMETER :: NS = 4
      INTEGER, PARAMETER :: NT = 4
      INTEGER, PARAMETER :: nTraj = 10
      INTEGER, PARAMETER :: seed_value = 24681357
      COMPLEX(WC) :: U(3,3,4,NT,NS,NS,NS)
      COMPLEX(WC) :: UNew(3,3,4,NT,NS,NS,NS)
      REAL(WP) :: aplaq, splaq, tplaq
      REAL(WP) :: ref_aplaq, ref_splaq, ref_tplaq
      REAL(WP) :: tol_aplaq, tol_splaq, tol_tplaq
      REAL(WP) :: sumTrP, time
      INTEGER  :: nPlaq, iTraj
      LOGICAL  :: found_ref
      CALL fill_identity_su3(U)
      UNew = U
      CALL seed_rng_fixed(seed_value)
      DO iTraj = 1, nTraj
         CALL updateLinks(U, beta, xi, UNew, trim(action_tag))
         U = UNew
      END DO
      CALL genPlaquette(U, NT, NS, NS, NS, 1, 4, 4, sumTrP, nPlaq, time)
      aplaq = sumTrP / (3.0_WP * real(nPlaq, kind=WP))
      IF (xi /= 1.0_WP) THEN
         CALL genPlaquette(U, NT, NS, NS, NS, 2, 4, 4, sumTrP, nPlaq, time)
         splaq = sumTrP / (3.0_WP * real(nPlaq, kind=WP))
         CALL genPlaquette(U, NT, NS, NS, NS, 1, 1, 4, sumTrP, nPlaq, time)
         tplaq = sumTrP / (3.0_WP * real(nPlaq, kind=WP))
      ELSE
         splaq = 0.0_WP
         tplaq = 0.0_WP
      END IF
      CALL read_heatbath_observable_reference( &
         "testdata/reference_values.toml", case_name, xi /= 1.0_WP, &
         ref_aplaq, tol_aplaq, ref_splaq, tol_splaq, ref_tplaq, tol_tplaq, found_ref)
      ! write(*,*) aplaq, splaq, tplaq
      CALL check(error, found_ref, &
         "heatbath observable reference [observables.heatbath."//trim(case_name)//"] must exist")
      IF (allocated(error)) RETURN
      ! Basic sanity first
      CALL check(error, ieee_is_finite(aplaq) .and. abs(aplaq) <= 1.0_WP + 1.0e-12_WP, &
         "average plaquette should be finite and bounded")
      IF (allocated(error)) RETURN
      IF (xi /= 1.0_WP) THEN
         CALL check(error, ieee_is_finite(splaq) .and. abs(splaq) <= 1.0_WP + 1.0e-12_WP, &
            "spatial plaquette should be finite and bounded")
         IF (allocated(error)) RETURN
         CALL check(error, ieee_is_finite(tplaq) .and. abs(tplaq) <= 1.0_WP + 1.0e-12_WP, &
            "temporal plaquette should be finite and bounded")
         IF (allocated(error)) RETURN
      END IF
      ! Reference-value regression
      CALL check(error, abs(aplaq - ref_aplaq) < tol_aplaq, &
         "average plaquette should match the stored heatbath reference")
      IF (allocated(error)) RETURN
      IF (xi /= 1.0_WP) THEN
         CALL check(error, abs(splaq - ref_splaq) < tol_splaq, &
            "spatial plaquette should match the stored heatbath reference")
         IF (allocated(error)) RETURN

         CALL check(error, abs(tplaq - ref_tplaq) < tol_tplaq, &
            "temporal plaquette should match the stored heatbath reference")
      END IF
      WRITE(*,'(A,1X,A)')  'Observable case:', trim(case_name)
      WRITE(*,'(A,1X,A)')  'Action tag:    ', trim(action_tag)
      WRITE(*,'(A,F10.6)') 'xi:            ', xi
      WRITE(*,'(A,F10.6)') 'beta:          ', beta
      WRITE(*,'(A,I0)')    'nTraj:         ', nTraj
      WRITE(*,'(A,F10.6)') 'aplaq:         ', aplaq
      IF (xi /= 1.0_WP) THEN
         WRITE(*,'(A,F10.6)') 'splaq:         ', splaq
         WRITE(*,'(A,F10.6)') 'tplaq:         ', tplaq
      END IF
    END SUBROUTINE run_heatbath_observable_case

   SUBROUTINE read_heatbath_observable_reference(filename, case_name, anisotropic, &
                                                 ref_aplaq, tol_aplaq, &
                                                 ref_splaq, tol_splaq, &
                                                 ref_tplaq, tol_tplaq, found)
      CHARACTER(len=*), INTENT(IN)  :: filename
      CHARACTER(len=*), INTENT(IN)  :: case_name
      LOGICAL,          INTENT(IN)  :: anisotropic
      REAL(WP),         INTENT(OUT) :: ref_aplaq, tol_aplaq
      REAL(WP),         INTENT(OUT) :: ref_splaq, tol_splaq
      REAL(WP),         INTENT(OUT) :: ref_tplaq, tol_tplaq
      LOGICAL,          INTENT(OUT) :: found

      TYPE(toml_table), ALLOCATABLE :: root
      TYPE(toml_table), POINTER     :: obs
      TYPE(toml_table), POINTER     :: hb
      TYPE(toml_table), POINTER     :: case_tbl
      TYPE(toml_error), ALLOCATABLE :: err
      INTEGER :: stat
      ref_aplaq = -huge(1.0_WP)
      tol_aplaq = -1.0_WP
      ref_splaq = -huge(1.0_WP)
      tol_splaq = -1.0_WP
      ref_tplaq = -huge(1.0_WP)
      tol_tplaq = -1.0_WP
      found = .FALSE.
      CALL toml_load(root, trim(filename), error=err)
      IF (allocated(err)) RETURN
      CALL get_value(root, "observables", obs, stat=stat)
      IF (stat /= 0 .or. .not. associated(obs)) RETURN
      CALL get_value(obs, "heatbath", hb, stat=stat)
      IF (stat /= 0 .or. .not. associated(hb)) RETURN
      CALL get_value(hb, trim(case_name), case_tbl, stat=stat)
      IF (stat /= 0 .or. .not. associated(case_tbl)) RETURN
      CALL get_value(case_tbl, "aplaq", ref_aplaq, stat=stat)
      IF (stat /= 0) RETURN
      CALL get_value(case_tbl, "aplaq_tol", tol_aplaq, stat=stat)
      IF (stat /= 0) RETURN
      IF (anisotropic) THEN
         CALL get_value(case_tbl, "splaq", ref_splaq, stat=stat)
         IF (stat /= 0) RETURN
         CALL get_value(case_tbl, "splaq_tol", tol_splaq, stat=stat)
         IF (stat /= 0) RETURN
         CALL get_value(case_tbl, "tplaq", ref_tplaq, stat=stat)
         IF (stat /= 0) RETURN
         CALL get_value(case_tbl, "tplaq_tol", tol_tplaq, stat=stat)
         IF (stat /= 0) RETURN
      END IF
      found = .TRUE.
   END SUBROUTINE read_heatbath_observable_reference

END MODULE test_heatbath_observables
