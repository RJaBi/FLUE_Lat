module FLUE_gluonProp
  use FLUE_constants, only: WP, WC, PI
  use FLUE_SU3MatrixOps, only: TracelessConjgSubtract, colourDecomp
  use FLUE_mom, only: get_qhat
  implicit none
  ! private

  ! public :: calc_mom_space_scalarD
    
  !public :: qhat, qhat_1, qhat_a, qhat_a4

contains
 
  pure subroutine calc_mom_space_scalarD(U, D, mu_start, xi)
    ! Fourier transformed U already...
    ! First half in each of space-time is pos freq
    ! second half is neg freq
    complex(kind=WC), dimension(:,:,:,:,:,:,:), intent(in)  :: U ! NT, NS, NS, NS, ND, NC, NC
    integer,                                    intent(in)  :: mu_start
    real(kind=WP),                              intent(in)  :: xi
    complex(kind=WC), dimension(:,:,:,:),       intent(out) :: D
    !
    integer, dimension(7) :: U_shape
    complex(kind=WC), dimension(3,3) :: U_pos, U_neg, A_pos, A_neg
    integer, dimension(4) :: coord
    real(kind=WP), dimension(4) :: qhat
    complex(kind=WC), dimension(8) :: com_pos, com_neg
    real(kind=WC) :: comSum
    real(kind=WP) :: prefac
    ! counters
    integer :: tt, qx, qy, qz, mu, aa

    ! Initialise these to zero
    ! Probably not actually needed...
    D = (0.0_WP, 0.0_WP)

    
    ! Get the size of the dimensions from U
    U_shape = shape(U)
    do concurrent(tt=0: U_shape(1)/2, qx=0: U_shape(2)/2, qy=0: U_shape(3)/2, qz=0: U_shape(4)/2, mu=mu_start: U_shape(5))
       ! Get the gaugefield
       coord = [tt, qx, qy, qz]
       call get_qhat(real(coord, kind=WP), U_shape, qhat, a=1.0_WP)
       U_pos = U(tt, qx, qy, qz, mu, :, :)
       U_neg = U(U_shape(1) - tt, U_shape(2) - qx, U_shape(3) - qy, U_shape(4) - qz, mu, :, :)
       ! TrSub = A - B^dagger - Tr(A-B^dagger) / 3.0
       call TraceLessConjgSubtract(A_pos, U_pos, U_neg)
       call TraceLessConjgSubtract(A_neg, U_neg, U_pos)
       ! (i/2) * exp(-pi Q/ N) * A_
       A_pos = (0.0_WP, 0.5_WP) * exp(-(0.0_WP, 1.0_WP) * 0.5_WP * qhat(mu) ) * A_pos
       A_neg = (0.0_WP, 0.5_WP) * exp( (0.0_WP, 1.0_WP) * 0.5_WP * qhat(mu) ) * A_neg
       ! Perform colour decomposition
       call colourDecomp(com_pos, A_pos)
       call colourDecomp(com_neg, A_neg)
       ! Update running two point-sum
       comSum = (0.0_WP, 0.0_WP)
       ! The reduction can not be done on array elements
       !do concurrent(aa=1:8) reduce(+:comSum)
       do aa=1, 8
          if (mu == 0) then  ! i.e. mu = t
             comSum = comSum + xi*xi * com_pos(aa) * com_neg(aa)
          else
             comSum = comSum + com_pos(aa) * com_neg(aa)     
          end if
       end do
       D(tt, qx, qy, qx) = comSum
    end do
    ! Normallise
    do concurrent(tt=0: U_shape(1)/2, qx=0: U_shape(2)/2, qy=0: U_shape(3)/2, qz=0: U_shape(4)/2)
       ! Determine pre-factor
       if (mu_start == 1 .and. qx + qy + qx < 1) then
          ! i.e. 1/8 * vol
          ! (nc*nc -1) * nd *NT*NS*NS*NS
          prefac = 2.0_WP / (8.0_WP * &
               real(U_shape(1), kind=WP) * &
               real(U_shape(2), kind=WP) * &
               real(U_shape(3), kind=WP) * &
               real(U_shape(4), kind=WP) * &
               real(U_shape(5), kind=WP))
       else if (mu_start == 0 .and. tt + qx + qy + qz + qx < 1) then
          ! i.e. 1/8 * vol
          ! (nc*nc -1) * nd *NT*NS*NS*NS
          prefac = 2.0_WP / (8.0_WP * &
               real(U_shape(1), kind=WP) * &
               real(U_shape(2), kind=WP) * &
               real(U_shape(3), kind=WP) * &
               real(U_shape(4), kind=WP) * &
               real(U_shape(5), kind=WP))
       else
          ! i.e. 1/8 * 3 * NT*NS*NS*NS
          ! (nc*nc -1) * (nd-1) *NT*NS*NS*NS
          prefac = 2.0_WP / (8.0_WP * &
               real(U_shape(1), kind=WP) * &
               real(U_shape(2), kind=WP) * &
               real(U_shape(3), kind=WP) * &
               real(U_shape(4), kind=WP) * 3.0_WP)
       end if
       ! apply
       D(tt, qx, qy, qx) = D(tt, qx, qy, qx) * prefac
    end do

  end subroutine calc_mom_space_scalarD


  
end module FLUE_gluonProp
