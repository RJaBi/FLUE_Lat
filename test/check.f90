program check
  use FLUE, only: WP, WC, &
       constructSU2Matrix, randomNumbers,  &
       constructSU3Matrix
  use stdlib_random, only: random_seed
  use stdlib_linalg, only: det
  use stdlib_math, only: is_close
  implicit none(external)
  ! Set up the random seed
  integer, parameter :: seedPut = 12345
  integer :: seedGet
  ! SU2 variables
  real(kind=WP), dimension(4) :: SU2_r
  complex(kind=WC), dimension(2,2) :: SU2_U
  real(kind=WP) :: SU2_det
  ! SU3 variables
  real(kind=WP), dimension(4,3) :: SU3_r
  complex(kind=WC), dimension(2,2) :: SU3_SU2_r, SU3_SU2_s, SU3_SU2_t
  complex(kind=WC), dimension(3,3) :: SU3_U
  real(kind=WP) :: SU3_det
  ! Do the seed
  call random_seed(put=seedPut, get=seedGet)
  ! now do some tests
  SU2_r = randomNumbers()
  SU2_U = constructSU2Matrix(SU2_r)
  SU2_det = det(SU2_U)
  if (.not. is_close(SU2_det, 1.0_WP)) then
     write(*,*) 'SU2_det', SU2_det, 'not close to 1.0'
     write(*,*) 'SU2_U is ', SU2_U
     stop
  end if
  !print *, "Put some tests in here!"
  write(*,*) 'SU2 checks done'
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Now do some SU3 checks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Get random numbers
  SU3_r(:, 1) = randomNumbers()
  SU3_r(:, 2) = randomNumbers()
  SU3_r(:, 3) = randomNumbers()
  ! Get some SU2 matrices
  SU3_SU2_r = constructSU2Matrix(SU3_r(:, 1))
  SU3_SU2_s = constructSU2Matrix(SU3_r(:, 2))
  SU3_SU2_t = constructSU2Matrix(SU3_r(:, 3))
  ! Do the SU3 matrix
  SU3_U = constructSU3Matrix(SU3_SU2_r, SU3_SU2_s, SU3_SU2_t)
  SU3_det = det(SU3_U)
  if (.not. is_close(SU3_det, 1.0_WP)) then
     write(*,*) 'SU3_det', SU3_det, 'not close to 1.0'
     write(*,*) 'SU3_U is ', SU3_U
     stop
  end if
  write(*,*) 'SU3 checks done'
end program check
