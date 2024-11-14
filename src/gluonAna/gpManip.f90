!! Module to do gluon propagator manipulation
module FLUE_gpManip
   use FLUE_constants, only: WP, pi
   !use stdlib_stats, only: median
   implicit none

contains

   subroutine Q_Average(NQ, q, D, q_slices, DOUT, QOUT, Qcount)
    !! Averages over adjacent q slices
    !! Averages both q and D
      ! Arguments
      integer, intent(in) :: NQ
    !! The total number of Q and D values
      real(kind=WP), dimension(NQ), intent(in) :: q, D
    !! The Q and D values to be averaged
      integer, intent(in) :: q_slices
    !! The number of adjacent to average over
      integer, intent(out) :: Qcount
    !! The final number of averaged values
      real(kind=WP), dimension(1:NQ + 1), intent(out) :: DOut, QOUT
    !! The output, averaged values. Note only the first QCount are used
      ! Some working variables
      real(kind=WP) :: dq, q_pos, maxq
      logical, dimension(NQ) :: mask
      QOUt = 0.0_wp
      DOUT = 0.0_wp
      ! setup the loop
      maxq = maxval(q)
      dq = maxq / real(q_slices, kind=WP)
      q_pos = 0.0_wp
      qcount = 1
      do
         ! Get a mask of 'adjacent' q values
         mask = q > q_pos .and. q < q_pos + dq
         if (any(mask)) then
            ! write(*,*) qcount, 'inside mask', q_pos
            ! If there is anything in the mask, then keep the output
            QOUT(qcount) = sum(q, mask=mask) / count(mask, kind=WP)
            !DOUT(qcount) = median(D, mask=mask)
            DOUT(qcount) = D(1)
            !write(*,*) D(merge(1,0,mask))
            !stop
            !DOUT(qcount) = median(D(mask), len(D(mask)))!sum(D, mask=mask) / count(mask, kind=WP)
            qcount = qcount + 1
         end if
         q_pos = q_pos + dq
         ! Check the loop condition
         if (q_pos + dq >= maxq) exit
      end do
      !write(*,*) 'fort', QOUT(:qcount)

   end subroutine Q_Average

   subroutine cone_cut(NQ, radius, Q, D, D4, NT, NS, angleIN, xiIN, IRCutIN, IRRadiusIN, QOut, Dout, D4Out, qcount)
    !! Performs a cone cut along the BCD axis of the data
      ! Arguments
      integer, intent(in) :: NQ
    !! The total number of Q and D, D4 values
      integer, intent(in) :: NT, NS
    !! lattice dimensions
      integer, intent(in) :: radius
    !! keeps points within this radius
      real(kind=WP), dimension(NQ, 4), intent(in) :: q
      !integer, dimension(NQ, 4), intent(in) :: q
    !! The Q values to be cut ove
      real(kind=WP), dimension(NQ), intent(in) :: D, D4
    !! The D, D4 values to be cut over
      integer, intent(out) :: Qcount
    !! the final number of values retained
      real(kind=WP), dimension(NQ, 4), intent(out) :: QOUT
    !! the output cut values of q. Note only the first Qcount are used
      real(kind=WP), dimension(NQ), intent(out) :: DOut, D4out
    !! the output cut values of D, D4. Note only the first Qcount are used
      ! Optional variables
      real(kind=WP), optional, intent(in) :: angleIN, xiIN, IRCutIN, IRRadiusIN
      real(kind=WP) :: angle, xi, IRCut, IRRadius
      ! working variables
      ! The cone_mask is identical in the second index
      logical, dimension(NQ, 4) :: cone_mask
      real(kind=WP), parameter, dimension(4) :: BCD_norm = (/0.5_wp, 0.5_wp, 0.5_wp, 0.5_wp/)
      real(kind=WP), dimension(4) :: qhat
      real(kind=WP) :: r, theta, q_norm
      integer :: qq
      ! set the optional vars
      if (present(angleIN)) then
         angle = angleIN
      else
         angle = pi / 2.0_wp
      end if
      if (present(xiIN)) then
         xi = xiIN
      else
         xi = 1.0_wp
      end if
      if (present(IRRadiusIN)) then
         IRRadius = IRRadiusIN
         !write(*,*) 'IRRadius', IRRadiusIN
      else
         IRRadius = radius
      end if
      if (present(IRCutIN)) then
         IRCut = IRCutIN
      else
         IRCut = 0.0_wp
      end if
      !write(*,*) IRRadius, IRCut, angle
      cone_mask = .false.
      do qq = 1, NQ
         !call get_qhat(q(qq,:), (/NT,NS,NS,NS/), qhat)
         !write(*,*) qq, q(qq, :), qhat
         qhat = q(qq, :)
         q_norm = sqrt(sum(real(qhat, kind=WP)**2.0_wp))

         !if (qq == -1) then
         !   write(*,*) 'coordbase', q(qq,:)
         !   write(*,*) 'qhat', qhat
         !   write(*,*) 'q_norm', q_norm
         !end if
         ! This doesn't really make sense
         ! Except i guess in small angle approx which we are in
         ! correct for anisotropy
         qhat(1) = qhat(1) * xi
         !if (qq == -1) then
         !   write(*,*) 'qhat *xi', qhat
         !end if
         if (all(qhat == qhat(1))) then
            ! handle on-diagonal coordinates
            r = 0
            theta = 0
         else
            theta = acos(dot_product(BCD_norm, qhat) / q_norm)
            r = q_norm * sin(theta)
         end if
         ! Do the checks to construct a mask
         if (q_norm <= IRCut) then
            cone_mask(qq, :) = r <= IRRadius
         else
            cone_mask(qq, :) = r <= IRRadius .and. theta < angle
         end if
         !if (qq == -1) then
         !   write(*,*) 'r', r
         !   write(*,*) 'theta', theta
         !   write(*,*) 'q_norm here', q_norm
         !   write(*,*) 'dot', dot_product(BCD_norm, qhat)
         !   write(*,*) 'BCD', BCD_norm
         !end if
      end do
      ! and now do the cut using the mask
      ! pack selets the elements which are true
      qcount = count(cone_mask(:, 1))
      qout = 0.0_wp
      dout = 0.0_wp
      d4out = 0.0_wp
      dout(:qcount) = pack(d, cone_mask(:, 1))
      qout(:qcount, 1) = pack(q(:, 1), cone_mask(:, 1))
      qout(:qcount, 2) = pack(q(:, 2), cone_mask(:, 1))
      qout(:qcount, 3) = pack(q(:, 3), cone_mask(:, 1))
      qout(:qcount, 4) = pack(q(:, 4), cone_mask(:, 1))
      d4out(:qcount) = pack(d4, cone_mask(:, 1))

   end subroutine cone_cut

end module FLUE_gpManip
