module mini_cloud_vf
  use mini_cloud_precision
  use mini_cloud_class
  implicit none

  !! Diameter, LJ potential and molecular weight for background gases
  real(dp), parameter :: d_He = 2.511e-8_dp, LJ_He = 10.22_dp * kb, molg_He = 4.002602_dp
  real(dp), parameter :: d_CH4 = 3.758e-8_dp, LJ_CH4 = 148.6_dp * kb, molg_CH4 = 16.0425_dp
  real(dp), parameter :: d_CO = 3.690e-8_dp, LJ_CO = 91.7_dp * kb, molg_CO = 28.0101_dp
  real(dp), parameter :: d_CO2 = 3.941e-8_dp, LJ_CO2 = 195.2_dp * kb, molg_CO2 = 44.0095_dp
  real(dp), parameter :: d_H2 = 2.827e-8_dp, LJ_H2 = 59.7_dp * kb, molg_H2 = 2.01588_dp
  real(dp), parameter :: d_H2O = 2.641e-8_dp, LJ_H2O = 809.1_dp * kb, molg_H2O = 18.01528_dp
  real(dp), parameter :: d_NH3 = 2.900e-8_dp, LJ_NH3 = 558.3_dp * kb, molg_NH3 = 17.03052_dp
  real(dp), parameter :: d_N2 = 3.798e-8_dp, LJ_N2 = 71.4_dp * kb, molg_N2 = 14.0067_dp
  real(dp), parameter :: d_H = 2.5e-8_dp, LJ_H =  30.0_dp * kb, molg_H = 1.00794 * amu

  ! Manual mixing parameters for background gas
  logical, parameter :: gmix = .True.
  integer, parameter :: ngmix = 2

  real(dp), dimension(ngmix) :: mixr_g = (/0.85_dp, 0.15_dp/)
  real(dp), dimension(ngmix) :: d_g = (/d_H2, d_He/)
  real(dp), dimension(ngmix) :: LJ_g = (/LJ_H2, LJ_He/)
  real(dp), dimension(ngmix) :: molg_g = (/molg_H2, molg_He/)
  real(dp), dimension(ngmix) :: nu_g

contains

  subroutine calc_vf(n_dust, T_in, P_in, mu_in, grav, k, k3, vf)
    implicit none

    integer, intent(in) :: n_dust
    real(dp), intent(in) :: T_in, P_in, mu_in, grav
    real(dp), dimension(4), intent(in) :: k
    real(dp), dimension(n_dust), intent(in) :: k3

    real(dp), intent(out) :: vf

    integer :: g, i, j, n
    real(dp), dimension(n_dust) :: b_mix
    real(dp) :: T, mu, P_cgs, a_av, rho,  rho_mix
    real(dp) :: nu_mix, nu_sum, phi_ij, phi_ij_top, phi_ij_bot
    real(dp) :: l_scale, beta, Kn

    T = T_in
    mu = mu_in
    P_cgs = P_in * 10.0_dp
    rho = (P_cgs * mu * amu)/(kb * T)

    a_av = max(k(2)/k(1),a_seed)
    a_av = min(a_av, 100.0_dp)

    ! Find the mixing ratio of the grain composition
    b_mix(:) = max(k3(:)/k(4),1e-30_dp)
    b_mix(:) = min(k3(:)/k(4),1.0_dp)
    
    ! Find weighted bulk density of mixed grain composition using volume fraction
    rho_mix = 0.0_dp
    do n = 1, n_dust
        rho_mix = rho_mix + b_mix(n) * d_sp(n)%bulk_den
    end do

    ! Find the dynamical viscosity of each gas
    do g = 1, ngmix
      ! Dynamical viscosity formula - Rosner (2000/2012) using Ackerman & Marley (2001) constants
      nu_g(g) = (5.0_dp/16.0_dp) * (sqrt(pi*(molg_g(g)*amu)*kb*T)/(pi*d_g(g)**2)) &
        & * (((kb*T)/LJ_g(g))**(0.16_dp)/1.22_dp)
    end do

    !! Now we loop twice over to find the mixture value using the Wilke (1950) mixing rule
    nu_mix = 0.0_dp
    do i = 1, ngmix
      nu_sum = 0.0_dp
      do j = 1, ngmix
        phi_ij_top = (1.0_dp + sqrt(nu_g(i)/nu_g(j)) * (molg_g(j)/molg_g(i))**(0.25_dp))**2
        phi_ij_bot = (4.0_dp/sqrt(2.0_dp)) * sqrt(1.0_dp + (molg_g(i)/molg_g(j)))
        phi_ij = phi_ij_top  / phi_ij_bot
        nu_sum = nu_sum + mixr_g(j) * phi_ij
      end do
      nu_mix = nu_mix + (mixr_g(i) * nu_g(i)) / nu_sum
    end do

    !! For consistency, we now use the dynamical viscosity to find the mean free path of the layer
    l_scale = (nu_mix/P_cgs) * sqrt((pi * kb * T) / (2.0_dp * mu * amu))

    !! Knudsen number and Cunningham slip factor at mean grain size
    Kn = l_scale / a_av
    beta = 1.0_dp + Kn * (1.256_dp + 0.4_dp * exp(-1.1_dp/Kn))

    ! Final v_f is negative (downward)
    vf = -(2.0_dp * beta * a_av**2 * grav * (rho_mix - rho)) / (9.0_dp * nu_mix)

  end subroutine calc_vf

end module mini_cloud_vf