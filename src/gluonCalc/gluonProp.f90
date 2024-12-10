! Where we actually calculate the gluon propagator
module FLUE_gluonProp
   use FLUE_constants, only: WP, WC, PI
   use FLUE_SU3MatrixOps, only: TracelessConjgSubtract, colourDecomp
   use FLUE_mom, only: get_qhat
   implicit none(external)
   ! private

   ! public :: calc_mom_space_scalarD

   !public :: qhat, qhat_1, qhat_a, qhat_a4

   public

contains

   subroutine scalarGluonProp(NS, NT, U, xi, D)
      ! U is in position space
      integer, intent(in) :: NS, NT
      complex(kind=WC), dimension(NT, NS, NS, NS, 4, 3, 3), intent(in) :: U
      real(kind=WP), intent(in) :: xi
      real(kind=WP), dimension(NS, NS, NS, NS), intent(out) :: D

      complex(kind=WC), dimension(NT, NS, NS, NS, 4, 3, 3) :: UMom  !, UMom

      complex(kind=WC), dimension(3, 3) :: T1, T2, Umu, trSub

      integer :: tt, xx, yy, zz, mu, aa
      integer :: tt1, xx1, yy1, zz1

      integer, dimension(4) :: coord
      real(kind=WP), dimension(4) :: qhat
      complex(kind=WC) :: prefactor

      complex(kind=WC), dimension(8) :: com_pos, com_neg
      complex(kind=WC) :: comSum

      ! hep-lat/9803105

      ! U_\mu(qhat) = \sum_{x} exp(-iqhat*x)U_\mu(x)
      ! for qhat = 2pin_\mu / aN_\mu
      ! for n_\mu = 0, 1, ..., N_\mu - 1
      ! So we will calculate qhat_\mu and -qhat_\mu at same time
      ! So loop over all space-time
      ! and A_mu(qhat) = (exp(-iqhat_mu*a/2)/2ig0) * ([U_\mu(qhat) - U_\mu(-qhat)&^\dagger] - 1/3Trace)
      write (*, *) "BLAH"
      UMom = (0.0_WP, 0.0_WP)
      UMom = (0.0_WP, 0.0_WP)
      do mu = 1, 4
         !do tt=1, NT
         !   do xx = 1, NS
         !      do yy = 1, NS
         !         do zz = 1, NS
         !D3IR$ PARALLEL ALWAYS ASSERT
       !!!do concurrent(tt=1:NT, xx=1:NS, yy=1:NS, zz=1:NS) default(none) local_init(tt,xx,yy,zz) local(coord, qhat, T1, T2, Umu, trSub, prefactor) shared(U, NS, NT)
         !$omp parallel do collapse(4) private(tt, xx, yy, zz, coord, qhat, T1, T2, Umu, trSub, prefactor) shared(U, NS, NT)
         do tt = 1, NT
            do xx = 1, NS
               do yy = 1, NS
                  do zz = 1, NS
                     coord = (/tt, xx, yy, zz/)
                     call get_qhat(real(coord, kind=WP), (/NT, NS, NS, NS/), qhat, a=1.0_WP)
                     ! write(*,*) 'coord', coord, 'qhat', qhat
                     Umu = U(tt, xx, yy, zz, mu, :, :)
                     T1 = (0.0_WP, 0.0_WP)
                     T2 = (0.0_WP, 0.0_WP)
                     ! Here we do the sum to make U_mu(qhat) and U_mu(-qhat)
                     do tt1 = 1, NT; do xx1 = 1, NS; do yy1 = 1, NS; do zz1 = 1, NS
                                 T1 = T1 + EXP(-(0.0_WP, 1.0_WP) * real(coord(mu), kind=WP) * qhat(mu)) * Umu
                                 T2 = T2 + EXP((0.0_WP, 1.0_WP) * real(coord(mu), kind=WP) * qhat(mu)) * Umu
                              end do; end do; end do; end do
                     ! and then calculate A_\mu
                     call TraceLessConjgSubtract(trSub, T1, T2)
                     UMom(tt, xx, yy, zz, mu, :, :) = trSub
                     call TraceLessConjgSubtract(trSub, T2, T1)
                     UMom(tt, xx, yy, zz, mu, :, :) = trSub
                     ! Pos mom
                     prefactor = (0.0_WP, 0.5_WP) * EXP(-(0.0_WP, 1.0_WP) * qhat(mu) * 0.5_WP)
                     UMom(tt, xx, yy, zz, mu, :, :) = prefactor * UMom(tt, xx, yy, zz, mu, :, :)
                     ! Neg mom
                     prefactor = (0.0_WP, 0.5_WP) * EXP((0.0_WP, 1.0_WP) * qhat(mu) * 0.5_WP)
                     UMom(tt, xx, yy, zz, mu, :, :) = prefactor * UMom(tt, xx, yy, zz, mu, :, :)
                  end do
               end do
            end do
         end do
         !$omp end parallel do
      end do

      ! Now that we have calculate A_mu(qhat), A_\mu(-qhat)
      do mu = 1, 4
         do tt = 1, NT
            do xx = 1, NS
               do yy = 1, NS
                  do zz = 1, NS
                     call colourDecomp(com_pos, UMom(tt, xx, yy, zz, mu, :, :))
                     call colourDecomp(com_neg, UMom(tt, xx, yy, zz, mu, :, :))
                     comSum = (0.0_WP, 0.0_WP)
                     do aa = 1, 8
                        if (mu == 0) then  ! i.e. mu = t
                           comSum = comSum + xi * xi * com_pos(aa) * com_neg(aa)
                        else
                           comSum = comSum + com_pos(aa) * com_neg(aa)
                        end if
                     end do
                     D(tt, xx, yy, zz) = D(tt, xx, yy, zz) + comSum
                  end do
               end do
            end do
         end do
      end do
      D = (1.0_WP / 3.0_WP) * (1.0_WP / 8.0_WP) * real(NT * NS * NS * NS, kind=WP) * D

   end subroutine scalarGluonProp

   subroutine calc_mom_space_scalarD(NS, NT, U, D, mu_start, xi)
      ! Fourier transformed U already...
      ! First half in each of space-time is pos freq
      ! second half is neg freq
      !complex(kind=WC), dimension(:,:,:,:,:,:,:), intent(in)  :: U ! NT, NS, NS, NS, ND, NC, NC
      integer, intent(in) :: NS, NT
      complex(kind=WC), dimension(0:NT, 0:NS, 0:NS, 0:NS, 0:3, 0:2, 0:2), intent(in) :: U  ! NT, NS, NS, NS, ND, NC, NC
      integer, intent(in) :: mu_start
      real(kind=WP), intent(in) :: xi
      !complex(kind=WC), dimension(:,:,:,:),       intent(out) :: D
      complex(kind=WC), dimension(0:NT/2, 0:NS/2, 0:NS/2, 0:NS/2), intent(out) :: D
      !
      integer, dimension(7) :: U_shape
      complex(kind=WC), dimension(3, 3) :: U_pos, U_neg, A_pos, A_neg
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
      write (*, *) 'U0', U(0, 0, 0, 0, 0, 0, 0)
      write (*, *) 'U0', U(1, 0, 0, 0, 0, 0, 0)
      !write(*,*) "HERE I AM\n\n\n\n\n"  , 'lowU', lbound(U), '\nlowD', lbound(D), '\nupU', ubound(U), '\nupD', ubound(D), '\nShapeU', shape(U)&
      !   , size(U), size(U(1, 1, 1, 1, 1, :, :))
      ! Get the size of the dimensions from U
      U_shape = SHAPE(U)
      !do concurrent(tt=0: U_shape(1)/2, qx=0: U_shape(2)/2, qy=0: U_shape(3)/2, qz=0: U_shape(4)/2, mu=mu_start: U_shape(5))
      ! write(*,*) '3,1'
      ! write(*,*) U(0,0,0,1,1, 3,1)
      ! write(*,*) '3,2'
      ! write(*,*) U(0,0,0,1,1, 3,2)
      ! write(*,*) '3,3'
      ! write(*,*) U(0,0,0,1,1, 3,3)
      ! write(*,*) '1,1'
      ! write(*,*) U(0,0,0,1,1, 1,1)
      ! write(*,*) '1,2'
      ! write(*,*) U(0,0,0,1,1, 1,2)
      ! write(*,*) '1,3'
      ! write(*,*) U(0,0,0,1,1, 1,3)
      ! write(*,*) '2,1'
      ! write(*,*) U(0,0,0,1,1, 2,1)
      ! write(*,*) '2,2'
      ! write(*,*) U(0,0,0,1,1, 2,2)
      ! write(*,*) '2,3'
      ! write(*,*) U(0,0,0,1,1, 2,3)
      D = (0.0_WP, 0.0_WP)
      do concurrent(tt=0:U_shape(1) / 2, qx=0:U_shape(2) / 2, qy=0:U_shape(3) / 2, qz=0:U_shape(4) / 2, mu=mu_start:U_shape(5) - 1)
         ! Get the gaugefield
         coord = [tt, qx, qy, qz]
         write (*, *) 'coord my', coord, mu
         call get_qhat(real(coord, kind=WP), U_shape, qhat, a=1.0_WP)
         ! write(*,*) 'qhat', qhat
         U_pos = U(tt, qx, qy, qz, mu, :, :)
         write (*, *) 'got U_pos'
         !write(*,*) U_shape(1) - tt - 1, U_shape(2) - qx - 1, U_shape(3) - qy - 1, U_shape(4) - qz - 1, mu
         !U_neg = U(U_shape(1) - tt - 1, U_shape(2) - qx - 1, U_shape(3) - qy - 1, U_shape(4) - qz - 1, mu, :, :)
         U_neg = U(U_shape(1) - tt, U_shape(2) - qx, U_shape(3) - qy, U_shape(4) - qz, mu, :, :)
         ! write(*,*) 'got U_neg'
         ! TrSub = A - B^dagger - Tr(A-B^dagger) / 3.0
         call TraceLessConjgSubtract(A_pos, U_pos, U_neg)
         ! write(*,*) 'done tr1'
         call TraceLessConjgSubtract(A_neg, U_neg, U_pos)
         ! write(*,*) 'done tr2'
         ! (i/2) * exp(-pi Q/ N) * A_
         A_pos = (0.0_WP, 0.5_WP) * EXP(-(0.0_WP, 1.0_WP) * 0.5_WP * qhat(mu)) * A_pos
         A_neg = (0.0_WP, 0.5_WP) * EXP((0.0_WP, 1.0_WP) * 0.5_WP * qhat(mu)) * A_neg
         ! Perform colour decomposition
         call colourDecomp(com_pos, A_pos)
         ! write(*,*) 'col1'
         call colourDecomp(com_neg, A_neg)
         ! write(*,*) 'col2'
         ! Update running two point-sum
         comSum = (0.0_WP, 0.0_WP)
         ! The reduction can not be done on array elements
         !do concurrent(aa=1:8) reduce(+:comSum)
         do aa = 1, 8
            if (mu == 0) then  ! i.e. mu = t
               comSum = comSum + xi * xi * com_pos(aa) * com_neg(aa)
            else
               comSum = comSum + com_pos(aa) * com_neg(aa)
            end if
         end do
         D(tt, qx, qy, qz) = D(tt, qx, qy, qz) + comSum

         if (tt == 0 .AND. qx == 0 .AND. qy == 0 .AND. qz == 0 .AND. mu == 0) then
            write (*, *) 'U', U_pos, 'N', U_neg
            write (*, *) 'A', A_pos, 'N', A_neg
            write (*, *) 'com', com_pos, 'N', com_neg
            write (*, *) 'sum', comSum
            !stop
         end if
      end do
      ! Normallise
      do concurrent(tt=0:U_shape(1) / 2, qx=0:U_shape(2) / 2, qy=0:U_shape(3) / 2, qz=0:U_shape(4) / 2)
         ! Determine pre-factor
         if (mu_start == 1 .AND. qx + qy + qx < 1) then
            ! i.e. 1/8 * vol
            ! (nc*nc -1) * nd *NT*NS*NS*NS
            prefac = 2.0_WP / (8.0_WP * &
                               real(U_shape(1), kind=WP) * &
                               real(U_shape(2), kind=WP) * &
                               real(U_shape(3), kind=WP) * &
                               real(U_shape(4), kind=WP) * &
                               real(U_shape(5), kind=WP))
         else if (mu_start == 0 .AND. tt + qx + qy + qz + qx < 1) then
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
         D(tt, qx, qy, qz) = D(tt, qx, qy, qz) * prefac
         if (tt + qx + qy + qz < 5) then
            write (*, *) tt, qx, qy, qz, D(tt, qx, qy, qz), prefac
         end if
      end do

   end subroutine calc_mom_space_scalarD

end module FLUE_gluonProp
