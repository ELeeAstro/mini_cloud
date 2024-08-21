module mini_cloud_saturation
  use mini_cloud_precision
  use mini_cloud_class
  implicit none

  public :: calc_saturation

contains

  subroutine calc_saturation(n_dust, VMR)
    implicit none

    integer, intent(in) :: n_dust
    real(dp), dimension(n_dust), intent(in) :: VMR

    integer :: n
    real(dp) :: TC, p_par, p_vap

    do n = 1, n_dust
    !! Calculate vapour pressure in dyne cm-2
    select case(d_sp(n)%name)
    case('C')
      p_vap = exp(3.27860e1_dp - 8.65139e4_dp/(T + 4.80395e-1_dp))
    !case('TiC')
    !case('SiC')
    !case('CaTiO3')
    case('TiO2')
      ! NIST 5 param fit
      p_vap = exp(-7.70443e4_dp/T + 4.03144e1_dp - 2.59140e-3_dp*T &
      & + 6.02422e-7_dp*T**2 - 6.86899e-11_dp*T**3)
    case('VO')
      ! NIST 5 param fit
      p_vap = exp(-6.74603e4_dp/T + 3.82717e1_dp - 2.78551e-3_dp*T &
      & + 5.72078e-7_dp*T**2 - 7.41840e-11_dp*T**3)
    case('Al2O3')
      ! Kozasa et al. (1989)
      p_vap = exp(-73503.0_dp/T + 22.01_dp) * atm
    case('Fe')
      ! Elspeth note: Changed to Ackerman & Marley et al. (2001) expression
      if (T > 1800.0_dp) then
        p_vap = exp(9.86_dp - 37120.0_dp/T) * bar
      else
        p_vap = exp(15.71_dp - 47664.0_dp/T) * bar
      end if
    case('FeO')
      ! NIST 5 param fit
      p_vap = exp(-6.30018e4_dp/T + 3.66364e1_dp - 2.42990e-3_dp*T &
      & + 3.18636e-7_dp*T**2 - 0.0_dp*T**3)
    case('Mg2SiO4')
      ! Kozasa et al. (1989)
      p_vap = exp(-62279_dp/T + 20.944_dp) * atm
    case('MgSiO3')
      ! A&M (2001)
      p_vap = exp(-58663.0_dp/T + 25.37_dp) * bar
    case('SiO2')
      ! NIST 5 param fit
      p_vap = exp(-7.28086e4_dp/T + 3.65312e1_dp - 2.56109e-4_dp*T &
      & - 5.24980e-7_dp*T**2 + 1.53343E-10_dp*T**3)
    case('SiO')
      ! Gail et al. (2013)
      p_vap = exp(-49520.0_dp/T + 32.52_dp)
    case('Cr')
      ! NIST 5 param fit
      p_vap = exp(-4.78455e4_dp/T + 3.22423e1_dp - 5.28710e-4_dp*T &
      & - 6.17347e-8_dp*T**2 + 2.88469e-12_dp*T**3)
    case('MnS')
      ! Morley et al. (2012)
      p_vap = 10.0_dp**(11.532_dp - 23810.0_dp/T) * bar
    case('Na2S')
      ! Morley et al. (2012)
      p_vap =  10.0_dp**(8.550_dp - 13889.0_dp/T) * bar
    case('ZnS')
      ! Morley et al. (2012)
      p_vap = 10.0_dp**(12.812_dp - 15873.0_dp/T) * bar
    case('KCl')
      ! NIST 5 param fit
      p_vap = exp(-2.69250e+4_dp/T + 3.39574e1_dp - 2.04903e-3_dp*T &
      & - 2.83957e-7_dp*T**2 + 1.82974e-10_dp*T**3)
    case('NaCl')
      ! NIST 5 param fit
      p_vap = exp(-2.79146e+4_dp/T + 3.46023e1_dp - 3.11287e-3_dp*T &
      & + 5.30965e-07_dp*T**2 - 2.59584e-12_dp*T**3)
    case('NH4Cl')
      p_vap = 10.0**(7.0220 - 4302.0/T) * bar
    case('H2O')
      ! Ackerman & Marley et al. (2001) H2O liquid & ice vapour pressure expressions
      TC = T - 273.15_dp
      if (T > 1048.0_dp) then
        p_vap = 6.0e8_dp
      else if (T < 273.16_dp) then
        p_vap = 6111.5_dp * exp((23.036_dp * TC - TC**2/333.7_dp)/(TC + 279.82_dp))
      else
        p_vap = 6112.1_dp * exp((18.729_dp * TC - TC**2/227.3_dp)/(TC + 257.87_dp))
      end if
    case('NH3')
      ! Ackerman & Marley et al. (2001) NH3 ice vapour pressure expression from Weast (1971)
      p_vap = exp(10.53_dp - 2161.0_dp/T - 86596.0_dp/T**2)  * bar
    case('CH4')
      !--- Prydz, R.; Goodwin, R.D., J. Chem. Thermodyn., 1972, 4,1 ---
      p_vap = 10.0_dp**(3.9895_dp - 443.028_dp/(T-0.49_dp)) * bar
    case('NH4SH')
      !--- E.Lee's fit to Walker & Lumsden (1897) ---
      p_vap = 10.0_dp**(7.8974_dp - 2409.4_dp/T) * bar
    case('H2S')
      !--- Stull (1947) ---
      if (T < 30.0_dp) then ! Limiter for very cold T
        p_vap = 10.0_dp**(4.43681_dp - 829.439_dp/(30.0_dp-25.412_dp)) * bar
      end if
      if (T < 212.8_dp) then
        p_vap = 10.0_dp**(4.43681_dp - 829.439_dp/(T-25.412_dp)) * bar
      else
        p_vap = 10.0_dp**(4.52887_dp - 958.587_dp/(T-0.539_dp)) * bar
      end if
    case('S2')
      !--- Zahnle et al. (2016) ---
      if (T < 413.0_dp) then
        p_vap = exp(27.0_dp - 18500.0_dp/T) * bar
      else
        p_vap = exp(16.1_dp - 14000.0_dp/T) * bar
      end if
    case('S8')
      !--- Zahnle et al. (2016) ---
      if (T < 413.0_dp) then
        p_vap = exp(20.0_dp - 11800.0_dp/T) * bar
      else
        p_vap = exp(9.6_dp - 7510.0_dp/T) * bar
      end if
    case default
      print*, 'Saturation: dust species not found: ', d_sp(n)%name, 'Stopping!'
      stop
    end select

    p_par = VMR(n) * P_cgs
    d_sp(n)%sat = p_par / p_vap

    ! Limiters
    d_sp(n)%sat = max(d_sp(n)%sat,1.0e-99_dp)
    d_sp(n)%sat = min(d_sp(n)%sat,1.0e30_dp)

    !print*, n, p_par, p_vap, d_sp(n)%sat

    end do

  end subroutine calc_saturation

end module mini_cloud_saturation
