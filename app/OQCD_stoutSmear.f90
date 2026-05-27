!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Read openqcd gaugefield, do some stout smearing
!! Ryan Bignell 2025
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program OQCD_to_ILDG
  use FLUE, only: WP, WC, writeCompiler, writeGit, &
      readGaugeField_OpenQCD, &
      genPlaquette, &
      StoutSmearLinks
  implicit none(type, external)
  ! IO vars
   character(len=256) :: inputFile, NTS, NSS, rhoS, nSmearS
  complex(kind=WC), dimension(:, :, :, :, :, :, :), allocatable :: U1, USmeared
  integer :: nSmear
  real(kind=WP) :: rho
real(kind=WP) :: plaq, time, plaqSmear
integer :: nplaq
  integer :: NT, NS

  call writeCompiler()
  call writeGit()
  write (*, *) ""

  if (COMMAND_ARGUMENT_COUNT() > 0) then
     call GET_COMMAND_ARGUMENT(1, inputFile)
     write (*, *) "Reading gauge file from ", TRIM(inputFile)
     call GET_COMMAND_ARGUMENT(2, rhoS)
     write (*, *) "Smearing with rho =  ", trim(rhoS)
      read(rhoS, *) rho
     call GET_COMMAND_ARGUMENT(3, nSmearS)
     write(*,*) "Smearing sweeps = ", trim(nSmearS)
      call str2int(nSmearS, nSmear)
     call GET_COMMAND_ARGUMENT(4, NTS)
     call str2int(NTS, NT)
     write (*, *) "NT IS", NT
     call GET_COMMAND_ARGUMENT(5, NSS)
     call str2int(NSS, NS)
     write (*, *) "NS IS", NS
  else
     write (*, *) "Pass the full path to the input toml on the command line"
     write (*, *) "i.e. fpm run OQCD_stoutSmear -- inputFile rho nSmear NT NS"
     stop
  end if

  allocate(U1(3, 3, 4, NT, NS, NS, NS))
  allocate(USmeared, mold=U1)
  U1 = ReadGaugeField_OpenQCD(trim(inputFile), NS, NS, NS, NT, fixSU3=.true.)

  ! First calculate the unsmeared plaquette
  call genPlaquette(U1, NT, NS, NS, NS, 1, 4, 4, plaq, nplaq, time)
  ! Now do some smearing
  call StoutSmearLinks(U1, rho, nSmear, USmeared)
  ! Now calculate smeared plaquette
  call genPlaquette(USmeared, NT, NS, NS, NS, 1, 4, 4, plaqSmear, nplaq, time)
  plaq = plaq / real(nplaq, kind=WP)
  plaqSmear = plaqSmear / real(nplaq, kind=WP)
  write(*,*) 'unsmeared'
  write(*,*) plaq
  write(*,*) 'smeared'
  write(*,*) plaqSmear
contains

  elemental subroutine str2int(str,int,stat)
    ! Modified from https://stackoverflow.com/a/24077338
    implicit none(external)
    ! Arguments
    character(len=*),intent(in) :: str
    integer,intent(out)         :: int
    integer,optional, intent(out)         :: stat
    character(len=25) :: mystr
    if (present(stat) )then
       read(str,*,iostat=stat)  int
       if (stat /= 0) then
          mystr = 'iostat was not equal to 0'
!          write(*,*) "iostat is ", stat
       end if
    else
       read(str,*)  int
    end if
  end subroutine str2int

end program OQCD_To_ILDG
