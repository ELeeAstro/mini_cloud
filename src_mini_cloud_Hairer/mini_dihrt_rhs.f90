module mini_dihrt_rhs
  use mini_dihrt_precision
  use mini_dihrt_class
  implicit none

  public :: form_RHS

contains

  subroutine form_RHS(n_dust, k, k3, VMR, f)
    implicit none

    integer, intent(in) :: n_dust
    real(dp), dimension(4), intent(inout) :: k
    real(dp), dimension(n_dust), intent(inout) :: k3, VMR

    real(dp), dimension(4+n_dust+n_dust), intent(out) :: f

    real(dp) :: J_net, chi_net

    J_net = sum(d_sp(:)%Js) + sum(d_sp(:)%sevap)
    chi_net = sum(d_sp(:)%chis)

    f(1) = J_net
    f(2) = sum(a_seed*(d_sp(:)%Js + d_sp(:)%sevap)) + chi_net*k(1)
    f(3) = sum(a_seed**2*(d_sp(:)%Js + d_sp(:)%sevap)) + 2.0_dp*chi_net*k(2)
    f(4) = sum(a_seed**3*(d_sp(:)%Js + d_sp(:)%sevap)) + 3.0_dp*chi_net*k(3)

    f(5:5+n_dust-1) = a_seed**3*(d_sp(:)%Js + d_sp(:)%sevap) + 3.0_dp*d_sp(:)%chis*k(3)

    f(5+n_dust:5+n_dust+n_dust-1) = - d_sp(:)%Nl*(d_sp(:)%Js + d_sp(:)%sevap) - 4.0_dp*pi/d_sp(:)%dV * k(3) * d_sp(:)%chis
    ! f(5+n_dust:5+n_dust+n_dust-1) = -fourpi3 * nd_atm/d_sp(:)%dV * (f(5:5+n_dust-1)/nd_atm - &
    ! & d_sp(:)%a0**3*(d_sp(:)%Js + d_sp(:)%sevap))

    f(5+n_dust:5+n_dust+n_dust-1) = f(5+n_dust:5+n_dust+n_dust-1)/nd_atm

    !print*, f(:)
    !stop

  end subroutine form_RHS


end module mini_dihrt_rhs
