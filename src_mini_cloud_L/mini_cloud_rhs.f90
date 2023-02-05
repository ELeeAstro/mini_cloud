module mini_cloud_rhs
  use mini_cloud_precision
  use mini_cloud_class
  implicit none

  public :: form_RHS

contains

  subroutine form_RHS(n_dust, LL, f)
    implicit none

    integer, intent(in) :: n_dust
    real(dp), dimension(4), intent(inout) :: LL

    real(dp), dimension(4+n_dust+n_dust), intent(out) :: f

    real(dp) :: J_net, chi_net

    J_net = sum(d_sp(:)%Js) + sum(d_sp(:)%sevap)
    chi_net = sum(d_sp(:)%chis)

    f(1) = J_net/rho
    f(2) = V_seed**(1.0_dp/3.0_dp)*J_net/rho + (1.0_dp/3.0_dp) * chi_net*LL(1)
    f(3) = V_seed**(2.0_dp/3.0_dp)*J_net/rho + (2.0_dp/3.0_dp) * chi_net*LL(2)
    f(4) = V_seed*J_net/rho + chi_net*LL(3)

    f(5:5+n_dust-1) = V_seed*(d_sp(:)%Js + d_sp(:)%sevap)/rho + 3.0_dp*d_sp(:)%chis*LL(3)

    f(5+n_dust:5+n_dust+n_dust-1) = -d_sp(:)%Nl*(d_sp(:)%Js + d_sp(:)%sevap) - &
    & rho * LL(3) * d_sp(:)%chis/d_sp(:)%dV

    f(5+n_dust:5+n_dust+n_dust-1) = f(5+n_dust:5+n_dust+n_dust-1)/nd_atm

  end subroutine form_RHS

end module mini_cloud_rhs
