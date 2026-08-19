!!!
!!! Run some SU3 heatbath
!!! Measure and print plaquette after each update
!!! Does not Save
!!!

PROGRAM SU3_heatbath
   USE FLUE, ONLY: WP, WC, writeCompiler, writeGit, &
                   Ident3x3, updateLinks, genPlaquette, build_colour_sites
   USE philox, ONLY: c64
   use stdlib_ascii, only: to_lower
   IMPLICIT NONE(TYPE, EXTERNAL)

   ! Lattice geometry
   INTEGER, PARAMETER :: NS = 12
   INTEGER, PARAMETER :: NT = 12
   COMPLEX(kind=WC), DIMENSION(3, 3, 4, NT, NS, NS, NS) :: U, UNew
   ! Sites linearisation
   INTEGER, ALLOCATABLE, DIMENSION(:, :, :) :: sites_t, sites_x, sites_y, sites_z
   INTEGER, ALLOCATABLE, DIMENSION(:, :) :: counts
   ! steps
   INTEGER, PARAMETER :: nTherm = 0
   INTEGER, PARAMETER :: nTraj = 20
   INTEGER, PARAMETER :: nSkip = 1
   ! simulation params
   REAL(kind=WP), PARAMETER :: beta = 2.8_WP
   real(kind=WP), parameter :: xi = 1.0_WP
   CHARACTER(len=20), PARAMETER :: actionTag = 'Iwasaki'
   logical :: useSymanzik
  !! counters
   INTEGER :: it, ix, iy, iz, mu
   INTEGER :: iTraj
  !! plaq
   REAL(kind=WP) :: sumTrP, time
   REAL(kind=WP) :: aplaq, splaq, tplaq
   INTEGER :: nPlaq
   !! Random
   integer(kind=C64), dimension(2), parameter :: key = (/1_C64, 2_C64/)

   CALL writeCompiler()
   CALL writeGit()

   WRITE (*, *) ' beta is ', beta
   WRITE(*,*) 'Lattice is ', NT , 'x', NS, '^3'
   WRITE(*,*) 'Action is ', trim(actionTag)
   WRITE(*,*) 'Doing', nTherm, 'Thermallisation sweeps'
   WRITE(*,*) 'Doing', nTraj, 'Production sweeps'
   !! Initialise to identity
   DO CONCURRENT(it=1:nt, ix=1:ns, iy=1:ns, iz=1:ns, mu=1:4)
      U(:, :, mu, it, ix, iy, iz) = Ident3x3
   END DO
   select case (to_lower(TRIM(actionTag)))
   case ('symanzik')
      useSymanzik = .TRUE.
   case ('wilson')
      useSymanzik = .FALSE.
   case ('iwasaki')
      useSymanzik = .true.
   end select
   write(*,*) '1'
   CALL build_colour_sites(NT, NS, NS, NS, useSymanzik, sites_t, sites_x, sites_y, sites_z, counts)
   write(*,*) '2'
   ! Thermallise
   DO iTraj = 1, nTherm
      CALL updateLinks(U, beta, UNew, key, iTraj, sites_t, sites_x, sites_y, sites_z, counts, TRIM(actionTag), xi=xi)
      U = UNew
   END DO
   write(*,*) '3'
   !shortcircuit!CALL genPlaquette(U, NT, NS, NS, NS, 1, 4, 4, sumTrp, nPlaq, time)
   !shortcircuit!aPlaq = sumTrP / (3.0_WP * real(nPlaq, kind=WP))
   !shortcircuit!IF (xi /= 1.0_WP) THEN
   !shortcircuit!   CALL genPlaquette(U, NT, NS, NS, NS, 2, 4, 4, sumTrp, nPlaq, time)
   !shortcircuit!   sPlaq = sumTrP / (3.0_WP * real(nPlaq, kind=WP))
   !shortcircuit!   CALL genPlaquette(U, NT, NS, NS, NS, 1, 1, 4, sumTrp, nPlaq, time)
   !shortcircuit!   tPlaq = sumTrP / (3.0_WP * real(nPlaq, kind=WP))
   !shortcircuit!   WRITE (*, *) 'after therm', aplaq, splaq, tplaq
   !shortcircuit!ELSE
   !shortcircuit!   WRITE (*, *) 'after therm', aplaq
   !shortcircuit!END IF
   DO iTraj = 1, nTraj
      write(*,*) 'about to updatelinks'
      CALL updateLinks(U, beta, UNew, key, iTraj + nTherm, sites_t, sites_x, sites_y, sites_z, counts, TRIM(actionTag), xi=xi)
      U = UNew
      IF (MOD(iTraj, nSkip) == 0) THEN
         CALL genPlaquette(U, NT, NS, NS, NS, 1, 4, 4, sumTrp, nPlaq, time)
         aPlaq = sumTrP / (3.0_WP * real(nPlaq, kind=WP))
         IF (xi /= 1.0_WP) THEN
            CALL genPlaquette(U, NT, NS, NS, NS, 2, 4, 4, sumTrp, nPlaq, time)
            sPlaq = sumTrP / (3.0_WP * real(nPlaq, kind=WP))
            CALL genPlaquette(U, NT, NS, NS, NS, 1, 1, 4, sumTrp, nPlaq, time)
            tPlaq = sumTrP / (3.0_WP * real(nPlaq, kind=WP))
            WRITE (*, *) 'traj', iTraj, 'plaq', aplaq, splaq, tplaq
         ELSE
            WRITE (*, *) 'traj', iTraj, 'plaq' , aplaq
         END IF
      END IF
   END DO
END PROGRAM SU3_heatbath
