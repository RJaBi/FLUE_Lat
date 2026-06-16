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
   INTEGER, PARAMETER :: NS = 16
   INTEGER, PARAMETER :: NT = 256
   COMPLEX(kind=WC), DIMENSION(3, 3, 4, NT, NS, NS, NS) :: U, UNew
   ! steps
   INTEGER, PARAMETER :: nTherm = 0
   INTEGER, PARAMETER :: nTraj = 1000
   INTEGER, PARAMETER :: nSkip = 1
   ! simulation params
   REAL(kind=WP), PARAMETER :: beta = 6.8_WP
   REAL(kind=WP), PARAMETER :: xi = 10_WP
   CHARACTER(len=20), PARAMETER :: actionTag='wilson'
  !! counters
   INTEGER :: it, ix, iy, iz, mu
   INTEGER :: iTraj
  !! plaq
   REAL(kind=WP) :: sumTrP, time
   REAL(kind=WP) :: aplaq, splaq, tplaq
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
      CALL updateLinks(U, beta, xi, UNew, trim(actionTag))
      U = UNew
   END DO

   CALL genPlaquette(U, NT, NS, NS, NS, 1, 4, 4, sumTrp, nPlaq, time)
   aPlaq = sumTrP / (3.0_WP * real(nPlaq, kind=WP))
   IF (xi /= 1.0_WP) THEN
      CALL genPlaquette(U, NT, NS, NS, NS, 2, 4, 4, sumTrp, nPlaq, time)
      sPlaq = sumTrP / (3.0_WP * real(nPlaq, kind=WP))
      CALL genPlaquette(U, NT, NS, NS, NS, 1, 1, 4, sumTrp, nPlaq, time)
      tPlaq = sumTrP / (3.0_WP * real(nPlaq, kind=WP))
      WRITE (*, *) 'after therm', aplaq, splaq, tplaq
   ELSE
      WRITE (*, *) 'after therm', aplaq
   END IF
   DO iTraj = 1, nTraj
      !write(*,*) iTraj
      CALL updateLinks(U, beta, xi, UNew, trim(actionTag))
      !write(*,*) 'updated'
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
