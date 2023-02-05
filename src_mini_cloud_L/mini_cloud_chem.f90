module mini_cloud_chem
  use mini_cloud_precision
  use mini_cloud_class
  implicit none

  public :: calc_chem

contains

  subroutine calc_chem(n_dust, LL, L3s, nd)
    implicit none

    integer, intent(in) :: n_dust
    real(dp), dimension(4), intent(in) :: LL
    real(dp), dimension(n_dust), intent(in) :: L3s, nd

    integer :: n
    real(dp) :: S, top_term, therm_term, bmix, k_fac, a_av
    real(dp) :: stab

    a_av = max(LL(2)/LL(1)* vol2rad,a_seed)
    a_av = min(a_av, 1.0_dp)

    do n = 1, n_dust

      k_fac = exp((2.0_dp * d_sp(n)%sig * d_sp(n)%dV)/(a_av * kb * T))

      stab = k_fac/d_sp(n)%sat
      
      if (stab <= 1.0_dp) then
        S = 1.0_dp - stab
      else
        bmix = max(L3s(n)/sum(L3s(:)),1e-30_dp)
        bmix = min(bmix,1.0_dp)
        S = (1.0_dp - stab) * bmix
      end if

      top_term = Afak * d_sp(n)%dV * d_sp(n)%alpha * nd(n) * d_sp(n)%v_stoi
      therm_term = sqrt((kb * T)/(twopi * d_sp(n)%mass))

      d_sp(n)%chis = top_term*therm_term * S

    end do

  end subroutine calc_chem

end module mini_cloud_chem
