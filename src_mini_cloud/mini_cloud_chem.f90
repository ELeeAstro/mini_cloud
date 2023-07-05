module mini_cloud_chem
  use mini_cloud_precision
  use mini_cloud_class
  implicit none

  public :: calc_chem

contains

  subroutine calc_chem(n_dust, k, k3, VMR)
    implicit none

    integer, intent(in) :: n_dust
    real(dp), dimension(4), intent(in) :: k
    real(dp), dimension(n_dust), intent(in) :: k3, VMR

    integer :: n
    real(dp) :: S, top_term, therm_term, bmix, k_fac, a_av
    real(dp) :: stab

    a_av = max(k(4)/k(3),a_seed)
    a_av = min(a_av, 100.0_dp)

    do n = 1, n_dust

      !! Kelvin Factor
      k_fac = exp((2.0_dp * d_sp(n)%sig * d_sp(n)%dV)/(a_av * kb * T))

      !! Stability factor
      stab = k_fac/d_sp(n)%sat
      
      if (stab <= 1.0_dp) then
        S = 1.0_dp - stab
      else
        bmix = max(k3(n)/k(4),1e-30_dp)
        bmix = min(bmix,1.0_dp)
        S = (1.0_dp - stab) * bmix
      end if

      top_term = d_sp(n)%dV * d_sp(n)%v_stoi * d_sp(n)%alpha * VMR(n) * P_cgs
      therm_term = sqrt(twopi * d_sp(n)%mass * kb * T)

      d_sp(n)%chis = top_term/therm_term * S

    end do

  end subroutine calc_chem

end module mini_cloud_chem
