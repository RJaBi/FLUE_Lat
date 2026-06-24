!!!
!!! Run some SU2 heatbath
!!! Measure and print plaquette after each update
!!! Does not Save
!!!

PROGRAM SU2_heatbath
   USE FLUE, ONLY: WP, WC, writeCompiler, writeGit, &
                   Ident2x2, SU2_updateLinks, SU2_genPlaquette
   USE philox, ONLY: c64
   IMPLICIT NONE(TYPE, EXTERNAL)

   ! Lattice geometry
   INTEGER, PARAMETER :: NS = 4
   INTEGER, PARAMETER :: NT = 4
   COMPLEX(kind=WC), DIMENSION(2, 2, 4, NT, NS, NS, NS) :: U, UNew
   ! steps
   INTEGER, PARAMETER :: nTherm = 0
   INTEGER, PARAMETER :: nTraj = 2000
   INTEGER, PARAMETER :: nSkip = 1
   ! simulation params
   REAL(kind=WP), PARAMETER :: beta = 1.60_WP
  !! counters
   INTEGER :: it, ix, iy, iz, mu
   INTEGER :: iTraj
  !! plaq
   REAL(kind=WP) :: sumTrP, time
   INTEGER :: nPlaq

   CALL writeCompiler()
   CALL writeGit()

  !! Initialise to identity
   DO CONCURRENT(it=1:nt, ix=1:ns, iy=1:ns, iz=1:ns, mu=1:4)
      U(:, :, mu, it, ix, iy, iz) = Ident2x2
   END DO

   ! Thermallise
   DO iTraj = 1, nTherm
      CALL SU2_updateLinks(U, beta, UNew, (/1_C64, 2_C64/), iTraj)
      U = UNew
   END DO

   CALL SU2_genPlaquette(U, NT, NS, NS, NS, 1, 4, 4, sumTrp, nPlaq, time)
   WRITE (*, *) 'after therm', sumTrP / real(nPlaq, kind=WP)
   DO iTraj = 1, nTraj
      !write(*,*) iTraj
      CALL SU2_updateLinks(U, beta, UNew, (/1_C64, 2_C64/), iTraj + nTherm)
      !write(*,*) 'updated'
      U = UNew
      !if (mod(iTraj, nSkip) == 0) then
      CALL SU2_genPlaquette(U, NT, NS, NS, NS, 1, 4, 4, sumTrp, nPlaq, time)
      WRITE (*, *) 'traj', iTraj, 'plaq', sumTrP / real(nPlaq, kind=WP)
      !end if
   END DO
END PROGRAM SU2_heatbath
