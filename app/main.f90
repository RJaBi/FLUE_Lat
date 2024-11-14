program main
  use FLUE, only: WP, WC, scalarGluonProp, ReadGaugeField_OpenQCD, writeCompiler, calc_mom_space_scalarD!, !calculate_area, calculate_perimeter
  use csv_module
  implicit none
  integer, parameter :: NS=32, NT=56
  complex(kind=WC), dimension(NT, NS, NS, NS, 4, 3, 3) :: U
  real(kind=WP), dimension(NT, NS, NS, NS) :: D
  complex(kind=WC), dimension(0:NT/2,0:NS/2,0:NS/2,0:NS/2) :: D2
  ! output
  type(csv_file) :: csvf
  integer :: tt, xx, yy, zz
  logical :: csvStatus

  call writeCompiler()
  
  !U = 1.0_WP
  U = ReadGaugeField_OpenQCD('/home/ryan/Documents/2024/conf/Gen2L/56x32/Gen2l_56x32n1', &
       NS, NS, NS, NT)
  write(*,*) 'Loaded gauge '
  !call scalarGluonProp(NS, NT, U, 3.5_WP, D)
  call calc_mom_space_scalarD(NS, NT, U, D2, 0, 3.5_WP)
  write(*,*) 'Done Prop'

  ! Set optional flags
  call csvf%initialize(verbose = .true.)
  ! open the file
  call csvf%open('Btest.csv', n_cols=6, status_ok=csvStatus)
  ! add header
  call csvf%add(['qt', 'qx', 'qy', 'qz'])
  call csvf%add(['D(q)_s'])
  call csvf%add(['D4(q)_s'])
  call csvf%next_row()
  do tt = 1, NT
     do xx = 1, NS
        do yy = 1, NS
           do zz = 1, NS
              call csvf%add((/tt, xx, yy, zz/))
              !call csvf%add((/D(tt,xx,yy,zz), D(tt,xx,yy,zz)/))
              call csvf%add((/D2(tt,xx,yy,zz), D2(tt,xx,yy,zz)/))
              call csvf%next_row()
           end do
        end do
     end do
  end do
  call csvf%close(status_ok=csvStatus)
  write(*,*) 'finished writing'
  
  
  !real(WP) :: area, perimeter
  !real(WP) :: radius = 3.0_WP
  
  !call calculate_area(radius, area)
  !call  calculate_perimeter(radius, perimeter)
  
  !print *, "Circle radius: ", radius
  !print *, "Circle area: ", area
  !print *, "Circle perimeter: ", perimeter
end program main
