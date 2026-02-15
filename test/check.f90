program check
  use FLUE, only: WP, WC, &
       constructSU2Matrix, randomNumbers
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
  write(*,*) 'check done'
end program check
