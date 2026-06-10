!!!
!!! Run some SU3 heatbath
!!! Measure and print plaquette after each update
!!! Does not Save
!!!

PROGRAM SU3_heatbath
   USE FLUE, ONLY: WP, WC, writeCompiler, writeGit, &
                   Ident3x3, updateLinks, genPlaquette
   USE stdlib_random, ONLY: random_seed
   IMPLICIT NONE(TYPE, EXTERNAL)

   ! Lattice geometry
   INTEGER, PARAMETER :: NS = 4
   INTEGER, PARAMETER :: NT = 4
   COMPLEX(kind=WC), DIMENSION(3, 3, 4, NT, NS, NS, NS) :: U, UNew
   ! steps
   INTEGER, PARAMETER :: nTherm = 0
   INTEGER, PARAMETER :: nTraj = 200
   INTEGER, PARAMETER :: nSkip = 1
   ! simulation params
   REAL(kind=WP), PARAMETER :: beta = 5.6_WP
   CHARACTER(len=20), PARAMETER :: actionTag='Symanzik'
  !! counters
   INTEGER :: it, ix, iy, iz, mu
   INTEGER :: iTraj
  !! plaq
   REAL(kind=WP) :: sumTrP, time
   INTEGER :: nPlaq

   CALL writeCompiler()
   CALL writeGit()

   WRITE (*, *) ' beta is ', beta
  !! Initialise to identity
   DO CONCURRENT(it=1:nt, ix=1:ns, iy=1:ns, iz=1:ns, mu=1:4)
      U(:, :, mu, it, ix, iy, iz) = Ident3x3
   END DO

   ! Thermallise
   DO iTraj = 1, nTherm
      CALL updateLinks(U, beta, UNew, trim(actionTag))
      U = UNew
   END DO

   CALL genPlaquette(U, NT, NS, NS, NS, 1, 4, 4, sumTrp, nPlaq, time)
   WRITE (*, *) 'after therm', sumTrP / (3.0_WP * real(nPlaq, kind=WP))
   DO iTraj = 1, nTraj
      !write(*,*) iTraj
      CALL updateLinks(U, beta, UNew, trim(actionTag))
      !write(*,*) 'updated'
      U = UNew
      IF (MOD(iTraj, nSkip) == 0) THEN
         CALL genPlaquette(U, NT, NS, NS, NS, 1, 4, 4, sumTrp, nPlaq, time)
         WRITE (*, *) 'traj', iTraj, 'plaq', sumTrP / (3.0_WP * real(nPlaq, kind=WP))
      END IF
   END DO
END PROGRAM SU3_heatbath
