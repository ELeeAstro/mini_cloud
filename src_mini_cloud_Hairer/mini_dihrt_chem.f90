module mini_dihrt_chem
  use mini_dihrt_precision
  use mini_dihrt_class
  implicit none


  public :: calc_chem

contains

  subroutine calc_chem(n_dust, k, k3, VMR)
    implicit none

    integer, intent(in) :: n_dust
    real(dp), dimension(4), intent(in) :: k
    real(dp), dimension(n_dust), intent(in) :: k3, VMR

    integer :: n
    real(dp) :: v_stoi, alpha, S, top_term, therm_term, bmix, p_par, k_fac, a_av
    real(dp), dimension(n_dust) :: nd

    nd(:) = nd_atm * VMR(:)

    v_stoi = 1.0_dp
    alpha = 1.0_dp

    a_av = max(k(2)/k(1),a_seed)
    a_av = min(a_av, 1.0_dp)

    do n = 1, n_dust

      k_fac = exp((2.0_dp * d_sp(n)%sig * d_sp(n)%dV)/(a_av * kb * T))

      S = (1.0_dp - k_fac/d_sp(n)%sat)
      
      if (d_sp(n)%sat < 1.0_dp) then
        bmix = max(abs(k3(n)/sum(k3(:))),1e-99_dp)
        bmix = min(bmix,1.0_dp)
        S = S * bmix
      end if

      top_term = d_sp(n)%dV * v_stoi * alpha * nd(n)
      therm_term = sqrt((kb * T)/(twopi * d_sp(n)%mass))

      d_sp(n)%chis = top_term*therm_term * S

    end do

  end subroutine calc_chem


end module mini_dihrt_chem
