program uzerobar
  use FLUE, only: WP, WC, writeCompiler, writeGit, &
       ReadGaugeField_ILDG, ReadGaugeField_OpenQCD, ReadGaugeTransformation_cola, &
       plaquette, RealTraceMat
  use tomlf, only : toml_table, toml_load, toml_array, get_value, toml_path
  use stdlib_io_npy, only: load_npy
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
  procedure(ReadTransInterface), pointer :: readTrans
  complex(kind=WC), dimension(:,:,:,:,:,:,:), allocatable :: U, UT
  complex(kind=WC), dimension(:,:,:,:,:,:), allocatable :: G
  ! lattice dimensions
  integer :: NT, NS
  ! counters
  integer :: ii, nn, mu
  ! plaquette params
  real(kind=WP) :: sumTrP, time
  real(kind=WP), dimension(:), allocatable :: aplaq, splaq, tplaq
  integer :: nP
  real(kind=WP) :: plaqFactor
  ! u0 Landau params
  real(kind=WP) :: uzFactor
  abstract interface
     function ReadGaugeInterface(filename, NX, NY, NZ, NT) result(U_xd)
       import :: WC
       character(len=*), intent(in) :: filename
       integer, intent(in) :: NX, NY, NZ, NT
       complex(kind=WC), dimension(NT, NX, NY, NZ, 4, 3, 3) :: U_xd
     end function ReadGaugeInterface
  end interface

  abstract interface
     function ReadTransInterface(filename, NX, NY, NZ, NT) result(G_xd)
       import :: WC
       character(len=*), intent(in) :: filename
       integer, intent(in) :: NX, NY, NZ, NT
       complex(kind=WC), dimension(NT, NX, NY, NZ, 3, 3) :: G_xd
     end function ReadTransInterface
  end interface

  

  call writeCompiler()
  call writeGit()
  write(*,*) ''

  
  if (command_argument_count() > 0) then
     call get_command_argument(1, tomlName)
     write(*,*) 'Reading toml file from ', trim(tomlName)
   else
      write (*, *) 'Pass the full path to the input toml on the command line'
      write (*, *) 'i.e. fpm run uzerobar -- mydir/input.toml'
      stop
   end if

   call toml_load(table, tomlName)
   call get_value(table, 'fixNum', nFixes)
   write(*,*) 'We will analyse ', nFixes, ' different cases'
   call get_value(table, 'fixLabels', top_array)
   allocate(fixLabels(nFixes))
   write(*,*) 'These are'
   do ii=1, nFixes
      call get_value(top_array, ii, strRead)
      fixLabels(ii) = strRead
      write(*,*) trim(fixLabels(ii))
   end do
   ! Now let's do the load and do the analysis
   ! Rather than pre-load and do all analysis at once
   allocate(aplaq(nFixes), splaq(nFixes), tplaq(nFixes))
   do ii=1, nFixes
      write(*,*) ''
      write(*,*) ''
      write(*,*) ii, fixLabels(ii)
      ! Get singleton variables
      call get_value(table, toml_path('fix', trim(fixLabels(ii)), 'gaugePath'), gaugePath)
      call get_value(table, toml_path('fix', trim(fixLabels(ii)), 'gaugeFormat'), gaugeFormat)
      call get_value(table, toml_path('fix', trim(fixLabels(ii)), 'transFormat'), transFormat)
      call get_value(table, toml_path('fix', trim(fixLabels(ii)), 'nTrans'), nTrans)
      call get_value(table, toml_path('fix', trim(fixLabels(ii)), 'NT'), NT)
      call get_value(table, toml_path('fix', trim(fixLabels(ii)), 'NS'), NS)
      ! Get the list of transform files
      call get_value(table, toml_path('fix', trim(fixLabels(ii)), 'transFiles'), top_array)
      allocate(transFiles(nTrans))
      do nn=1, nTrans
         call get_value(top_array, nn, strRead)
         transFiles(nn) = strRead
      end do
      ! Now do stuff
      select case (gaugeFormat)
      case ('cssmILDG')
         ReadGauge => ReadGaugeField_ILDG
      case ('openqcd')
         readGauge => ReadGaugeField_OpenQCD
      case default
         write(*,*) 'gaugeFormat ', trim(gaugeFormat), ' not supported. Exiting'
         stop
      end select
      ! allocate space for the gaugefield
      allocate(U(NT,NS,NS,NS,4,3,3), UT(NT,NS,NS,NS,4,3,3))
      U = ReadGauge(gaugePath, NS, NS, NS, NT)

      ! Now for the gauge transform
      select case (transFormat)
      case ('cssm')
         readTrans => ReadGaugeTransformation_cola
         plaqFactor = 9.0_WP
         uzFactor = 1.0_WP
      case ('npy')
         readTrans => npyFunctionWrapper
         plaqFactor = 1.0_WP
         uzFactor = 3.0_WP
      case default
         write(*,*) 'transFormat ', trim(transFormat), ' not supported. Exiting'
         stop
      end select
      allocate(G(NT,NS,NS,NS,3,3))

      write(*,*) 'For gauge format ', trim(gaugeFormat)
      write(*,*) 'Factor on plaquette is ', plaqFactor
      write(*,*) 'Factor on uzerobar is  ', uzFactor

      call plaquette(U, 1, 4, 4, sumTrP, nP, time)
      aplaq(ii) = plaqFactor * sumTrp / nP
      call plaquette(U, 1, 2, 4, sumTrP, nP, time)
      splaq(ii) = plaqFactor * sumTrp / nP
      call plaquette(U, 1, 1, 4, sumTrP, nP, time)
      tplaq(ii) = plaqFactor * sumTrp / nP
      write (*, *) 'U  Plaquette for', ' is ', aplaq(ii), 'and took', time, 'seconds'
      write (*, *) 'U SPlaquette for', ' is ', splaq(ii), 'and took', time, 'seconds'
      write (*, *) 'U TPlaquette for', ' is ', tplaq(ii), 'and took', time, 'seconds'
      write (*, *) 'u0bar is (1/3*aplaq)**0.25', ((1.0_WP / 3.0_WP) * aplaq(ii))**0.25_WP
      write (*, *) 'u0bar^s is (1/3*splaq)**0.25', ((1.0_WP / 3.0_WP) * splaq(ii))**0.25_WP
      write (*, *) 'u0bar^t is (1/3*tplaq/u0s)**0.5', &
           ((1.0_WP / 3.0_WP) * tplaq(ii) / (((1.0_WP / 3.0_WP) * splaq(ii))**0.25_WP)**2.0_WP)**0.5_WP
      do mu=1, 4
         write(*,*) mu, uzFactor * uzeroLandau(U, NS, NS, NS, NT, mu)
      end do
      ! Read transform, apply,
      do nn=1, nTrans
         G = ReadTrans(transFiles(nn), NS, NS, NS, NT)
         write(*,*) 'Apply transform', trim(transFiles(nn))
         UT = applyGauge(U, G, NS, NS, NS, NT)
         U = UT
         call plaquette(U, 1, 4, 4, sumTrP, nP, time)
         aplaq(ii) = plaqFactor * sumTrp / nP
         !write (*, *) 'U Plaquette for', ' is ', plaq, 'and took', time, 'seconds'
         write (*, *) 'u0bar is (1/3*aplaq)**0.25', ((1.0_WP / 3.0_WP) * aplaq(ii))**0.25_WP
         do mu=1, 4
            write(*,*) mu, uzFactor * uzeroLandau(U, NS, NS, NS, NT, mu)
         end do
      end do

      deallocate(transFiles)
      deallocate(U, UT)
      deallocate(G)            
   end do
 contains

   pure function uzeroLandau(U, NX, NY, NZ, NT, mu) result(uz)
     complex(kind=WC), dimension(NT, NX, NY, NZ, 4, 3, 3), intent(in) :: U
     integer, intent(in) :: NX, NY, NZ, NT, mu
     real(kind=WP) :: uz
     ! counters
     integer :: tt, xx, yy, zz
     real(kind=WP) :: uzLoop

     uz = 0.0_WP
     do concurrent(tt=1:NT, xx=1:NX, yy=1:NY, zz=1:NZ)  !reduction here but gfortran doesnt support
        call RealTraceMat(uzLoop, U(tt,xx,yy,zz,mu,:,:))
        uz = uz + uzLoop
     end do
     uz = uz / real(NT*NX*NY*NZ, kind=WP)
   end function uzeroLandau

   pure function applyGauge(U, G, NX, NY, NZ, NT) result(UT)
     complex(kind=WC), dimension(NT, NX, NY, NZ, 4, 3, 3), intent(in) :: U
     complex(kind=WC), dimension(NT, NX, NY, NZ, 3, 3), intent(in) :: G
     integer, intent(in) :: NX, NY, NZ, NT
     complex(kind=WC), dimension(NT, NX, NY, NZ, 4, 3, 3) :: UT
     ! counters
     integer :: m, i, j, k, l, cc
     integer, dimension(4) :: muCoord, coord
     integer, dimension(4) :: dataShape
     dataShape = (/NT, NX, NY, NZ/)
     ! Apply gauge transform
     do m = 1, 4
        muCoord(:) = 0
        muCoord(m) = 1
        do i = 1, NT
           do j = 1, NX
              do k = 1, NY
                 do l = 1, NZ
                    coord = (/i, j, k, l/) + muCoord
                    !# respect periodic boundary conditions
                    do cc = 1, SIZE(coord)
                       if (coord(cc) > dataShape(cc)) then
                          coord(cc) = 1
                       end if
                    end do
                    
                    UT(i, j, k, l, m, :, :) = MATMUL(G(i, j, k, l, :, :), &
                         MATMUL(U(i, j, k, l, m, :, :), &
                         CONJG(TRANSPOSE(G(coord(1), coord(2), coord(3), coord(4), :, :)))))
                 end do
              end do
           end do
        end do
     end do
   end function applyGauge
   

   function npyFunctionWrapper(filename, NX, NY, NZ, NT) result(G_x)
      character(len=*), intent(in) :: filename
      integer, intent(in) :: NX, NY, NZ, NT
      complex(kind=WC), dimension(NT, NX, NY, NZ, 3, 3) :: G_x
      complex(kind=WC), dimension(:, :, :, :, :, :), allocatable :: GNPY
      !counters
      integer :: xx, yy, zz, tt, aa, bb
      call load_npy(filename, GNPY)
      ! For some reason, load_npy loads the array in opposite order..
      do concurrent(tt=1:NT, xx=1:NS, yy=1:NS, zz=1:NS, aa=1:3, bb=1:3)
         G_x(tt,xx,yy,zz,aa,bb) = GNPY(aa,bb,zz,yy,xx,tt)
      end do
      !G_x = GNPY
    end function npyFunctionWrapper
      
     

end program uzerobar
