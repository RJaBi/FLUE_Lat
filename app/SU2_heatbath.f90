!!!
!!! Run some SU2 heatbath
!!! Measure and print plaquette after each update
!!! Does not Save
!!!

program SU2_heatbath
   use FLUE, only: WP, WC, writeCompiler, writeGit, &
                   Ident2x2, SU2_updateLinks, SU2_genPlaquette
   use stdlib_random, only: random_seed
   implicit none(type, external)

   ! Lattice geometry
   integer, parameter :: NS = 4
   integer, parameter :: NT = 4
   complex(kind=WC), dimension(2, 2, 4, NT, NS, NS, NS) :: U, UNew
   ! steps
   integer, parameter :: nTherm = 0
   integer, parameter :: nTraj = 2000
   integer, parameter :: nSkip = 1
   ! simulation params
   real(kind=WP), parameter :: beta = 1.60_WP
  !! counters
   integer :: it, ix, iy, iz, mu
   integer :: iTraj
  !! plaq
   real(kind=WP) :: sumTrP, time
   integer :: nPlaq

   call writeCompiler()
   call writeGit()

  !! Initialise to identity
   do concurrent(it=1:nt, ix=1:ns, iy=1:ns, iz=1:ns, mu=1:4)
      U(:, :, mu, it, ix, iy, iz) = Ident2x2
   end do

   ! Thermallise
   do iTraj = 1, nTherm
      call SU2_updateLinks(U, beta, UNew)
      U = UNew
   end do

   call SU2_genPlaquette(U, NT, NS, NS, NS, 1, 4, 4, sumTrp, nPlaq, time)
   write (*, *) 'after therm', sumTrP / real(nPlaq, kind=WP)
   do iTraj = 1, nTraj
      !write(*,*) iTraj
      call SU2_updateLinks(U, beta, UNew)
      !write(*,*) 'updated'
      U = UNew
      !if (mod(iTraj, nSkip) == 0) then
      call SU2_genPlaquette(U, NT, NS, NS, NS, 1, 4, 4, sumTrp, nPlaq, time)
      write (*, *) 'traj', iTraj, 'plaq', sumTrP / real(nPlaq, kind=WP)
      !end if
   end do
end program SU2_heatbath
