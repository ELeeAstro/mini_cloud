module mini_cloud_nucleation
  use mini_cloud_precision
  use mini_cloud_class
  implicit none

  public :: calc_nucleation, calc_seed_evap

contains

  subroutine calc_nucleation(n_dust, nd)
    implicit none

    integer, intent(in) :: n_dust
    real(dp), dimension(n_dust), intent(in) :: nd
    
    integer :: n

    do n = 1, n_dust
      
      if ((d_sp(n)%sat <= 1.0_dp) .or. (d_sp(n)%inuc == 0)) then

        call calc_sig(n, T, d_sp(n)%sig)
        d_sp(n)%Js = 0.0_dp

      else

        select case(d_sp(n)%name)
        case('SiO')
          call nuc_SiO(T, d_sp(n)%sat, nd(n), d_sp(n)%Js)
        case default
          call calc_sig(n, T, d_sp(n)%sig)
          call nuc_MCNT(T, d_sp(n)%sat, nd(n), d_sp(n)%Nf, d_sp(n)%sig, &
          & d_sp(n)%a0, d_sp(n)%alpha, d_sp(n)%mol_wght, d_sp(n)%Js)
        end select

        if (d_sp(n)%Js < 1e-10_dp) then
          d_sp(n)%Js = 0.0_dp
        end if

        !d_sp(n)%Js = 0.0_dp
        
      end if

      ! print*, n, T, d_sp(n)%sat, nd(n), d_sp(n)%Nf, d_sp(n)%sig, &
      ! & d_sp(n)%a0, d_sp(n)%alpha, d_sp(n)%mol_wght, d_sp(n)%Js
      ! stop

    end do

  end subroutine calc_nucleation

  ! Modified classical nucleation theory,
  ! ref: Gail & Sedlmayr (1984), Helling & Fomins (2014), Lee et al. (2015), Lee et al. (2018)
  subroutine nuc_MCNT(TT, SS, nsp, Nf, sig, a0, alpha, gmol, J_s)
    implicit none

    real(dp), intent(in) :: TT, SS, nsp
    real(dp), intent(in) :: Nf, sig, a0, alpha, gmol

    real(dp), intent(out) :: J_s

    real(dp) :: ln_ss, theta_inf, N_inf, N_star, N_star_1, dg_rt, Zel, tau_gr
    real(dp) :: f0, kbT

    ! Efficency Variables
    ln_ss = log(SS) ! Natural log of saturation ratio
    f0 = 4.0_dp * pi * a0**2 ! Monomer Area
    kbT = kb * TT         ! kb * T

    !! Find Nstar, theta_inf -> Dg/RT (Eq. 11 Lee et al. (2015a))
    theta_inf = (f0 * sig)/(kbT)  !Theta_infty eq.8 (Lee et al. (2015a)
    N_inf = (((twothird) * theta_inf) / ln_ss)**3
    !Gail et al. (2014) ! note, typo in Lee et al. (2015a)
    N_star = 1.0_dp + (N_inf / 8.0_dp) &
      & * (1.0_dp + sqrt(1.0_dp + 2.0_dp*(Nf/N_inf)**third) &
      & - 2.0_dp*(Nf/N_inf)**third)**3
    N_star = max(1.000000001_dp, N_star) ! Ensure no div 0
    N_star_1 = N_star - 1.0_dp

    !! delta G (N_star) / RT (Eq. 11 Lee et al. (2015a))
    dg_rt = theta_inf * (N_star_1 / (N_star_1**third + Nf**third))

    !! Calculate Zeldovich factor at N_star (Eq. 9 Lee et al. (2015a))
    Zel = sqrt((theta_inf / (9.0_dp * pi * (N_star_1)**(4.0_dp/3.0_dp))) &
      & * ((1.0_dp + 2.0_dp*(Nf/N_star_1)**third)/(1.0_dp + (Nf/N_star_1)**third)**3))

    !! Calculate !!inverse!! of tau_gr
    tau_gr = (f0 * N_star**(twothird)) * alpha * sqrt(kbT &
      & / (2.0_dp * pi * gmol * amu)) * nsp

    !! Finally calculate J_star ! Note underfloat limiter here
    J_s = nsp * tau_gr * Zel * exp(max(-300.0_dp, N_star_1*ln_ss - dg_rt))
    
  end subroutine nuc_MCNT

  subroutine calc_seed_evap(n_dust, LL)
    implicit none

    integer, intent(in) :: n_dust
    real(dp), dimension(4), intent(inout) :: LL

    integer :: n
    real(dp) :: tau_evap, a_av

    a_av = max(LL(2)/LL(1)*vol2rad,a_seed)
    a_av = min(a_av, 1.0_dp)

    do n = 1, n_dust

      if ((d_sp(n)%chis >= 0.0_dp) .or. (d_sp(n)%inuc == 0)) then

        d_sp(n)%sevap = 0.0_dp

      else

        tau_evap = a_av/abs(d_sp(n)%chis)
        d_sp(n)%sevap = -((LL(1)*rho)/tau_evap)

        if (d_sp(n)%sevap > -1e-10_dp) then
          d_sp(n)%sevap = 0.0_dp
        end if

        !d_sp(n)%sevap = 0.0_dp

      end if

    end do

  end subroutine calc_seed_evap

  ! SiO nucleation
  ! ref: Gail et al. (2013)
  pure subroutine nuc_SiO(TT, SS, nsp, J_s)
    implicit none

    real(dp), intent(in) :: TT, SS, nsp
    real(dp), intent(out) :: J_s

    J_s = nsp**2 * exp(1.33_dp - 4.40e12_dp / (TT**3 * log(SS)**2))

  end subroutine nuc_SiO

  pure subroutine calc_sig(n, T, sig)
    implicit none

    integer, intent(in) :: n
    real(dp), intent(in) :: T

    real(dp), intent(out) :: sig

    real(dp) :: TC

    TC = T - 273.15_dp

    select case (d_sp(n)%name)
    case('TiO2')
      ! Sindel et al. (2022)
      sig = 589.79_dp - 0.0708_dp * T
    case('Fe')
      sig = 1862.0_dp - 0.39_dp * (TC - 1530.0_dp)
      ! Pradhan et al. (2009)
      !sig = 2858.0_dp - 0.51_dp * T
    case('Al2O3')
      ! Pradhan et al. (2009)
      sig = 1024.0_dp - 0.177_dp * T
    case('Cr')
      sig = 1642.0_dp - 0.20_dp * (TC - 1860.0_dp)
    case('SiO2')
      ! Pradhan et al. (2009)
      sig = 243.2_dp - 0.013_dp * T

      ! Pradhan et al. (2009):
      !Si : 732 - 0.086*(T - 1685.0)
      !MgO : 1170 - 0.636*T
      !CaO : 791 - 0.0935*T
    case('KCl')
      sig = 160.4_dp - 0.070_dp*TC
    case('NaCl')
      sig = 171.5_dp - 0.0719_dp*TC
    case('H2O')
      sig = 141.0_dp - 0.15_dp*TC
    case default
      sig = d_sp(n)%sig
    end select

  end subroutine calc_sig

end module mini_cloud_nucleation
