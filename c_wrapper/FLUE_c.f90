module FLUE_c
  use iso_c_binding, only: c_double, c_int, c_double_complex
  use FLUE, only:  Q_Average, get_qhat, cone_cut, MultiplyMatMat, calc_mom_space_scalarD
  !calculate_area, calculate_perimeter,
contains


  subroutine calc_mom_space_scalarD_c(U, D, mu_start, xi)
    complex(kind=c_double_complex), dimension(:,:,:,:,:,:,:), intent(in)  :: U
    integer(kind=c_int),                                      intent(in)  :: mu_start
    real(kind=c_double),                                      intent(in)  :: xi
    complex(kind=c_double_complex), dimension(:,:,:,:),       intent(out) :: D

    call calc_mom_space_scalarD(U, D, mu_start, xi)
  end subroutine calc_mom_space_scalarD_c
    
  subroutine MultiplyMatMat_c(MM, left, right)
    complex(kind=c_double_complex), dimension(3, 3), intent(in) :: left, right
    complex(kind=c_double_complex), dimension(3, 3), intent(out) :: MM

    call MultiplyMatMat(MM, left, right)
  end subroutine MultiplyMatMat_c
    
  
  subroutine Q_Average_c(NQ, q, D, q_slices, DOUT, QOUT, Qcount)
    integer(c_int), intent(in) :: NQ
    real(c_double), intent(in), dimension(NQ) :: q, D
    integer(c_int), intent(in) :: q_slices
    integer(c_int), intent(out) :: Qcount
    real(c_double), intent(out), dimension(1:NQ+1) :: Dout, Qout

    call Q_Average(NQ, q, D, q_slices, DOUT, QOUT, Qcount)

  end subroutine Q_Average_c

  subroutine get_qhat_c(coord, shape, qhat, a)
    real(c_double), dimension(4), intent(in)  :: coord
    integer(c_int), dimension(4), intent(in)  :: shape
    real(c_double), dimension(4), intent(out) :: qhat
    real(c_double), optional,     intent(in)  :: a

    call get_qhat(coord, shape, qhat, a)
  end subroutine get_qhat_c


  subroutine cone_cut_c(NQ, radius, Q, D, D4, NT, NS, angleIN, xiIN, IRCutIN, IRRadiusIN, QOUT, Dout, D4out, qcount)
    integer(c_int), intent(in) :: NQ
    integer(c_int), intent(in) :: NT, NS
    integer(c_int), intent(in) :: radius
    real(kind=C_DOUBLE), dimension(NQ, 4), intent(in) :: q
    real(kind=C_DOUBLE), dimension(NQ), intent(in) :: D, D4
    integer(c_int), intent(out) :: Qcount
    real(kind=C_DOUBLE), dimension(NQ, 4), intent(out) :: QOUT
    real(kind=C_DOUBLE), dimension(NQ), intent(out) :: DOut, D4out
    real(kind=C_DOUBLE), optional, intent(in) :: angleIN, xiIN, IRCutIN, IRRadiusIN

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
 
