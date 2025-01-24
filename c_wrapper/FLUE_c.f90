module FLUE_c
   use ISO_C_BINDING, only: C_DOUBLE, C_INT, C_DOUBLE_COMPLEX
   use FLUE, only: Q_Average, get_qhat, cone_cut, MultiplyMatMat, calc_mom_space_scalarD, scalarGluonProp
   implicit none(external)
   public

contains

   subroutine scalarGluonProp_c(NS, NT, U, xi, D)
      ! U is in position space
      integer(kind=C_INT), intent(in) :: NS, NT
      complex(kind=C_DOUBLE_COMPLEX), dimension(NT, NS, NS, NS, 4, 3, 3), intent(in) :: U
      real(kind=C_DOUBLE), intent(in) :: xi
      real(kind=C_DOUBLE), dimension(NS, NS, NS, NS), intent(out) :: D
      call scalarGluonProp(NS, NT, U, xi, D)
   end subroutine scalarGluonProp_c

   subroutine calc_mom_space_scalarD_c(NS, NT, U, D, mu_start, xi)
      integer(kind=C_INT), intent(in) :: NS
      integer(kind=C_INT), intent(in) :: NT
      !complex(kind=c_double_complex), dimension(:,:,:,:,:,:,:), intent(in)  :: U
      complex(kind=C_DOUBLE_COMPLEX), dimension(0:NT, 0:NS, 0:NS, 0:NS, 0:3, 0:2, 0:2), intent(in) :: U
      integer(kind=C_INT), intent(in) :: mu_start
      real(kind=C_DOUBLE), intent(in) :: xi
      !complex(kind=c_double_complex), dimension(:,:,:,:),       intent(out) :: D
      complex(kind=C_DOUBLE_COMPLEX), dimension(0:NT/2, 0:NS/2, 0:NS/2, 0:NS/2), intent(out) :: D

      call calc_mom_space_scalarD(NS, NT, U, D, mu_start, xi)
   end subroutine calc_mom_space_scalarD_c

   subroutine MultiplyMatMat_c(MM, left, right)
      complex(kind=C_DOUBLE_COMPLEX), dimension(3, 3), intent(in) :: left, right
      complex(kind=C_DOUBLE_COMPLEX), dimension(3, 3), intent(out) :: MM

      call MultiplyMatMat(MM, left, right)
   end subroutine MultiplyMatMat_c

   subroutine Q_Average_c(NQ, q, D, q_slices, DOUT, QOUT, Qcount)
      integer(C_INT), intent(in) :: NQ
      real(C_DOUBLE), intent(in), dimension(NQ) :: q, D
      integer(C_INT), intent(in) :: q_slices
      integer(C_INT), intent(out) :: Qcount
      real(C_DOUBLE), intent(out), dimension(1:NQ + 1) :: Dout, Qout

      call Q_Average(NQ, q, D, q_slices, DOUT, QOUT, Qcount)

   end subroutine Q_Average_c

   subroutine get_qhat_c(coord, shape, qhat, a)
      real(C_DOUBLE), dimension(4), intent(in) :: coord
      integer(C_INT), dimension(4), intent(in) :: shape
      real(C_DOUBLE), dimension(4), intent(out) :: qhat
      real(C_DOUBLE), optional, intent(in) :: a

      call get_qhat(coord, shape, qhat, a)
   end subroutine get_qhat_c

   subroutine cone_cut_c(NQ, radius, Q, D, D4, NT, NS, angleIN, xiIN, IRCutIN, IRRadiusIN, QOUT, Dout, D4out, qcount)
      integer(C_INT), intent(in) :: NQ
      integer(C_INT), intent(in) :: NT, NS
      integer(C_INT), intent(in) :: radius
      real(kind=C_DOUBLE), dimension(NQ, 4), intent(in) :: q
      real(kind=C_DOUBLE), dimension(NQ), intent(in) :: D, D4
      integer(C_INT), intent(out) :: Qcount
      real(kind=C_DOUBLE), dimension(NQ, 4), intent(out) :: QOUT
      real(kind=C_DOUBLE), dimension(NQ), intent(out) :: DOut, D4out
      !real(kind=C_DOUBLE), optional, intent(in) :: angleIN, xiIN, IRCutIN, IRRadiusIN
      real(kind=C_DOUBLE), intent(in) :: angleIN, xiIN, IRCutIN, IRRadiusIN

      call cone_cut(NQ, radius, Q, D, D4, NT, NS, angleIN, xiIN, IRCutIN, IRRadiusIN, Qout, Dout, D4out, qcount)
   end subroutine cone_cut_c

!circ  subroutine calculate_area_c(radius, area)
!circ    !use iso_c_binding, only: c_double
!circ    real(c_double), intent(in) :: radius
!circ    real(c_double), intent(out) :: area
!circ
!circ    call calculate_area(radius, area)
!circ  end subroutine calculate_area_c
!circ
!circ  subroutine calculate_perimeter_c(radius, perimeter)
!circ    !use iso_c_binding, only: c_double
!circ    real(c_double), intent(in) :: radius
!circ    real(c_double), intent(out) :: perimeter
!circ
!circ    call calculate_perimeter(radius, perimeter)
!circ  end subroutine calculate_perimeter_c
end module FLUE_c
