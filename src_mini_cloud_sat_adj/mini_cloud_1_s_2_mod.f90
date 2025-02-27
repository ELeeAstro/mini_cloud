module mini_cloud_1_s_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  integer, parameter :: dp = REAL64

  !! Input values required from GCM
  real(dp) :: rho_c
  real(dp) :: tau_chem
  real(dp) :: rm 
  real(dp) :: sigma 
  real(dp) :: mol_w_sp
  real(dp) :: Kzz_deep
  real(dp) :: q_v_deep

  real(dp) :: pl, p_vap, q_s, tau_deep, cp, L, mw, Tl
  real(dp) :: delt

  real(dp), parameter :: bar = 1.0e6_dp ! bar to dyne
  real(dp), parameter :: atm = 1.01325e6_dp ! atm to dyne
  real(dp), parameter :: pa = 10.0_dp ! pa to dyne
  real(dp), parameter :: mmHg = 1333.22387415_dp  ! mmHg to dyne

  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp) ! value of pi
  real(dp), parameter :: twopi = pi * 2.0_dp

  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: amu = 1.66053906660e-24_dp ! g mol-1 (note, has to be cgs g mol-1 !!!)
  real(dp), parameter :: R_gas = 8.31446261815324e7_dp


  public :: mini_cloud_1_s
  private :: p_vap_sp, l_heat_sp, dqdt

contains 

  subroutine mini_cloud_1_s(deep_flag, latent_flag, T_in, P_in, grav_in, mu_in, cp_in, &
    & t_end, sp, q_v, q_c, dTl)
    implicit none

    logical, intent(in) :: deep_flag, latent_flag
    character(len=20), intent(in) :: sp
    real(dp), intent(in) :: t_end
    real(dp), intent(in) :: T_in, P_in, grav_in, mu_in, cp_in
    
    real(dp), intent(inout) :: q_v, q_c, dTl

    real(dp) :: mu, grav
    real(dp) :: Hp, eps

    real(dp) :: dt, t_now, dt_max

    integer :: n_it, n, accept, ierr
    integer, parameter :: nq = 3

    real(dp), dimension(nq) :: y, yc, y_in, y_new, y_em
    real(dp), dimension(nq) :: k1, k2, k3

    real(dp), dimension(nq) :: tol, err
    real(dp), parameter ::  pow = 0.2_dp, safe = 0.9_dp
    real(dp), parameter ::  atol = 1e-30_dp, rtol = 1e-3_dp

    !! Convert inputs to cgs units
    Tl = T_in
    pl = P_in * 10.0_dp ! Pa to dyne
    mu = mu_in 
    grav = grav_in * 100.0_dp ! m s-2 to cm s-2
    cp = cp_in * 1e4_dp ! J kg-1 K-1 to erg g-1 K-1

    !! Calculate deep replenishment rate of vapour
    Hp = (kb * Tl) / (mu * amu * grav)
    tau_deep = Hp**2/Kzz_deep

    !! Conversion factor between MMR and VMR (VMR -> MMR) for cloud species
    eps = mol_w_sp/mu

    !! Vapour pressure and saturation vapour pressure
    p_vap = p_vap_sp(sp, Tl)

    !! Latent heat release from sublimation/vapourisation
    if (latent_flag .eqv. .True.) then
      L = l_heat_sp(sp, Tl)
    else
      L = 0.0_dp
    end if

    !! Equilibrium (sat = 1) vapour pressure fraction
    q_s = min(1.0_dp, p_vap/pl)

    !! Calculate timescale with latent heat lag term
    delt = tau_chem * (1.0_dp + (L**2*q_s*mol_w_sp)/(cp*R_gas*Tl**2))

    !! initial conditions - convert to VMR from MMR
    y(1) = q_v/eps
    y(2) = q_c/eps

    ! Initial atmospheric temperature
    y(3) = Tl

    !! Copy input array to work array
    yc(:) = y(:)

    !! Start Runge-Kutta timestepping
    dt_max = t_end
    dt = dt_max/5.0_dp

    !! Begin timestepping routine
    t_now = 0.0_dp
    n_it = 0

    do while ((t_now < t_end) .and. (n_it < 100000))

      !! If next time step overshoots - last time step is equal tend
      if ((t_now + dt) > t_end) then
        dt = t_end - t_now
      end if

      !! Euler method
      !call dqdt(nq,yc,k1,deep_flag)
      !y_new(:) = yc(:) +  dt * k1(:)

      !! Bogacki–Shampine (order 3) Runge-Kutta method
      call dqdt(nq,yc,k1,deep_flag)
      y_in(:) = yc(:) + 0.5_dp * dt * k1(:)
      call dqdt(nq,y_in,k2,deep_flag)
      y_in(:) = yc(:) + 0.75_dp * dt * k2(:)
      call dqdt(nq,y_in,k3,deep_flag)
      y_new(:) = yc(:) + dt * (2.0_dp/9.0_dp * k1(:) + 1.0_dp/3.0_dp * k2(:) + 4.0_dp/9.0_dp * k3(:))

      !print*, y_new(:), y_em(:)

      yc(:) = y_new(:)

      t_now = t_now + dt
      n_it = n_it + 1

    end do

    !! Convert work array to final answer array
    y(:) = yc(:)

    !! Final results, convert back to MMR
    q_v = y(1) * eps
    q_c = y(2) * eps

    ! Change in atmospheric temperature from latent heat
    dTl = y(3) - Tl 

  end subroutine mini_cloud_1_s

  !! Calculate tracer fluxes
  subroutine dqdt(n,y,f,deep_flag)
    implicit none

    logical , intent(in) :: deep_flag
    integer, intent(in) :: n

    real(dp), dimension(n), intent(inout) :: y

    real(dp), dimension(n), intent(out) :: f

    real(dp) :: sat

    !y(1) = q_v, y(2) = q_c, y(3) = Tl

    y(1) = max(y(1),1e-30_dp)
    y(2) = max(y(2),1e-30_dp)
    y(3) = max(y(3),1e-30_dp)

    !! Calculate dqdt for vapour and condensate
    sat = max((y(1) * pl)/p_vap, 1e-99_dp)

    !! Calculate dqdt given the supersaturation ratio
    if (sat < 1.0_dp) then

      ! Evaporate from q_c
      f(1) = min(q_s - y(1), y(2))/delt

    else if (sat > 1.0_dp) then

      ! Condense q_v toward the saturation ratio
      f(1) = -(y(1) - q_s)/delt

    else
      f(1) = 0.0_dp
    end if

    f(2) = -f(1)

    ! Add replenishment to lower boundary at the tau_deep rate
    if (deep_flag .eqv. .True.) then
      f(1) = f(1) - (y(1) - q_v_deep)/tau_deep
    end if

    !! Change in atmospheric temperature due to latent heat
    f(3) = L/cp * f(2)

  end subroutine dqdt

  !! Vapour pressure for each species
  function p_vap_sp(sp, T) result(p_vap)
    implicit none

    character(len=20), intent(in) :: sp
    real(dp), intent(in) :: T

    real(dp) :: p_vap

    real(dp) :: TC

    ! Return vapour pressure in dyne
    select case(trim(sp))
    case('C')
      ! Kimura et al. (2023)
      p_vap = 10.0_dp**(-41523.0_dp/T + 10.609_dp) * atm
    case('TiC')
      ! Kimura et al. (2023)
      p_vap = 10.0_dp**(-33600.0_dp/T + 7.652_dp) * atm
    case('SiC')
      ! Elspeth 5 polynomial JANAF-NIST fit
      p_vap =  exp(-9.51431385e4_dp/T + 3.72019157e1_dp + 1.09809718e-3_dp*T &
        & -5.63629542e-7_dp*T**2 + 6.97886017e-11_dp*T**3)
    case('CaTiO3')
      ! Kozasa et al. (1987)
      p_vap = exp(-79568.2_dp/T + 42.0204_dp) * atm  
    case('Al2O3')
      ! Kozasa et al. (1989)
      p_vap = exp(-73503.0_dp/T + 22.01_dp) * atm      
    case('TiO2')
      ! GGChem 5 polynomial NIST fit
      p_vap = exp(-7.70443e4_dp/T +  4.03144e1_dp - 2.59140e-3_dp*T &
        &  + 6.02422e-7_dp*T**2 - 6.86899e-11_dp*T**3)
    case('VO')
      ! GGChem 5 polynomial NIST fit
      p_vap =  exp(-6.74603e4_dp/T + 3.82717e1_dp - 2.78551e-3_dp*T &
        & + 5.72078e-7_dp*T**2 - 7.41840e-11_dp*T**3)         
    case('Fe')
      ! Elspeth note: Changed to Ackerman & Marley et al. (2001) expression
      if (T > 1800.0_dp) then
        p_vap = exp(9.86_dp - 37120.0_dp/T) * bar
      else
        p_vap = exp(15.71_dp - 47664.0_dp/T) * bar
      end if
    case('FeS')
      ! GGChem 5 polynomial NIST fit
      p_vap = exp(-5.69922e4_dp/T + 3.86753e1_dp - 4.68301e-3_dp*T &
        & + 1.03559e-6_dp*T**2 - 8.42872e-11_dp*T**3)
    case('FeO')
      ! GGChem 5 polynomial NIST fit
      p_vap = exp(-6.30018e4_dp/T + 3.66364e1_dp - 2.42990e-3_dp*T &
        & + 3.18636e-7_dp*T**2)     
    case('MgS')
      ! Elspeth 5 polynomial JANAF-NIST fit
      p_vap = exp(-5.92010440e4_dp/T +  3.58148138e1_dp - 1.93822353e-3_dp*T &
        &  + 9.77971406e-7_dp*T**2 - 11.74813601e-10_dp*T**3)   
    case('Mg2SiO4')
      ! Kozasa et al. (1989)
      p_vap = exp(-62279_dp/T + 20.944_dp) * atm
    case('MgSiO3','MgSiO3_amorph')
      ! Ackerman & Marley (2001)
      p_vap = exp(-58663.0_dp/T + 25.37_dp) * bar
    case('MgO')
      ! GGChem 5 polynomial NIST fit
      p_vap = exp(-7.91838e4_dp/T + 3.57312e1_dp + 1.45021e-4_dp*T &
        & - 8.47194e-8_dp*T**2 + 4.49221e-12_dp*T**3)
    case('SiO2','SiO2_amorph')
      ! NIST 5 param fit
      p_vap = exp(-7.28086e4_dp/T + 3.65312e1_dp - 2.56109e-4_dp*T &
      & - 5.24980e-7_dp*T**2 + 1.53343E-10_dp*T**3) * bar
    case('SiO')
      ! Gail et al. (2013)
      p_vap = exp(-49520.0_dp/T + 32.52_dp)
    case('Cr')
      ! GGChem 5 polynomial NIST fit
      p_vap = exp(-4.78455e4_dp/T + 3.22423e1_dp - 5.28710e-4_dp*T & 
        &  - 6.17347e-8_dp*T**2 + 2.88469e-12_dp*T**3)
    case('MnS')
      ! Morley et al. (2012)
      p_vap = 10.0_dp**(11.532_dp - 23810.0_dp/T) * bar
    case('Na2S')
      ! Morley et al. (2012)
      p_vap =  10.0_dp**(8.550_dp - 13889.0_dp/T) * bar
    case('ZnS')
      ! Elspeth 5 polynomial Barin data fit
      p_vap = exp(-4.75507888e4_dp/T + 3.66993865e1_dp - 2.49490016e-3_dp*T &
        &  + 7.29116854e-7_dp*T**2 - 1.12734453e-10_dp*T**3)
    case('KCl')
      ! GGChem 5 polynomial NIST fit
      p_vap = exp(-2.69250e4_dp/T + 3.39574e1_dp - 2.04903e-3_dp*T &
          & -2.83957e-7_dp*T**2 + 1.82974e-10_dp*T**3)
    case('NaCl')
      ! GGChem 5 polynomial NIST fit
      p_vap = exp(-2.79146e4_dp/T + 3.46023e1_dp - 3.11287e3_dp*T & 
        & + 5.30965e-7_dp*T**2 -2.59584e-12_dp*T**3)
    case('NH4Cl')
      p_vap = 10.0_dp**(7.0220_dp - 4302.0_dp/T) * bar
    case('H2O')
      ! Ackerman & Marley (2001) H2O liquid & ice vapour pressure expressions
      TC = T - 273.15_dp
      if (T > 1048.0_dp) then
        p_vap = 6.0e8_dp
      else if (T < 273.16_dp) then
        p_vap = 6111.5_dp * exp((23.036_dp * TC - TC**2/333.7_dp)/(TC + 279.82_dp))
      else
        p_vap = 6112.1_dp * exp((18.729_dp * TC - TC**2/227.3_dp)/(TC + 257.87_dp)) 
      end if
    case('NH3')
      ! Ackerman & Marley (2001) NH3 ice vapour pressure expression from Weast (1971)
      p_vap = exp(10.53_dp - 2161.0_dp/T - 86596.0_dp/T**2)  * bar
    case('CH4')
      !--- Prydz, R.; Goodwin, R.D., J. Chem. Thermodyn., 1972, 4,1 ---
      if (T < 0.5_dp) then ! Limiter for very very very cold T
        p_vap = 10.0_dp**(3.9895_dp - 443.028_dp/(0.5_dp-0.49_dp)) * bar       
      else
        p_vap = 10.0_dp**(3.9895_dp - 443.028_dp/(T-0.49_dp)) * bar
      end if
    case('NH4SH')
      !--- E.Lee's fit to Walker & Lumsden (1897) ---
      p_vap = 10.0_dp**(7.8974_dp - 2409.4_dp/T) * bar
    case('H2S')
      !--- Stull (1947) ---
      if (T < 30.0_dp) then ! Limiter for very cold T
        p_vap = 10.0_dp**(4.43681_dp - 829.439_dp/(30.0_dp-25.412_dp)) * bar
      else if (T < 212.8_dp) then
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
    case('CO')
      ! Yaws
      p_vap = 10.0_dp**(51.8145e0_dp - 7.8824e2_dp/T - 2.2734e1_dp*log10(T) &
        & + 5.1225e-2_dp*T + 4.6603e-11_dp*T**2) * mmHg
    case('CO2')
      ! Yaws
      p_vap = 10.0_dp**(35.0187e0_dp - 1.5119e3_dp/T - 1.1335e1_dp*log10(T) &
          & + 9.3383e-3_dp*T + 7.7626e-10_dp*T**2) * mmHg
    case('H2SO4')
      ! GGChem 5 polynomial NIST fit
      p_vap = exp(-1.01294e4_dp/T + 3.55465e1_dp - 8.34848e-3_dp*T)           
    case default
      print*, 'Saturation: dust species not found: ', trim(sp), 'Stopping!'
      stop
    end select

  end function p_vap_sp

  !! latent heat for each species
  function l_heat_sp(sp, T) result(L_heat)
    implicit none

    character(len=20), intent(in) :: sp
    real(dp), intent(in) :: T

    real(dp) :: L_heat


    ! Return latent heat in erg g-1
    ! Get value from vapour pressure expression (or special function)
    select case(trim(sp))
    case('C')
      L_heat = 41523.0_dp * log(10.0_dp) * R_gas / mol_w_sp
    case('TiC')
      L_heat = 33600.0_dp * log(10.0_dp) * R_gas / mol_w_sp
    case('SiC')
      L_heat = 9.51431385e4_dp * R_gas / mol_w_sp
    case('CaTiO3')
      L_heat = 79568.2_dp * R_gas / mol_w_sp
    case('Al2O3')
      L_heat = 73503.0_dp * R_gas / mol_w_sp
    case('TiO2')
      L_heat = 7.70443e4_dp * R_gas / mol_w_sp
    case('VO')
      L_heat = 6.74603e4_dp * R_gas / mol_w_sp
    case('Fe')
      L_heat = 37120.0_dp * R_gas / mol_w_sp
    case('FeS')
      L_heat = 5.69922e4_dp * R_gas / mol_w_sp
    case('FeO')
      L_heat = 6.30018e4_dp * R_gas / mol_w_sp
    case('MgS')
      L_heat = 5.92010440e4_dp * R_gas / mol_w_sp
    case('Mg2SiO4')
      L_heat = 62279.0_dp * R_gas / mol_w_sp
    case('MgSiO3','MgSiO3_amorph')
      L_heat = 58663.0_dp * R_gas / mol_w_sp
    case('MgO')
      L_heat = 7.91838e4_dp * R_gas / mol_w_sp
    case('SiO2','SiO2_amorph')
      L_heat = 7.28086e4_dp * R_gas / mol_w_sp
    case('SiO')
      L_heat = 49520.0_dp * R_gas / mol_w_sp
    case('Cr')
      L_heat = 4.78455e4_dp * R_gas / mol_w_sp
    case('MnS')
      L_heat = 23810.0_dp * log(10.0_dp) * R_gas / mol_w_sp
    case('Na2S')
      L_heat = 13889.0_dp * log(10.0_dp) * R_gas / mol_w_sp
    case('ZnS')
      L_heat = 4.75507888e4_dp * R_gas / mol_w_sp
    case('KCl')
      L_heat = 2.69250e4_dp * R_gas / mol_w_sp
    case('NaCl')
      L_heat = 2.79146e4_dp * R_gas / mol_w_sp
    case('NH4Cl')
      L_heat = 4302.0_dp * log(10.0_dp) * R_gas / mol_w_sp
    case('H2O')
      L_heat = 2257.0e7_dp
    case('NH3')
      L_heat = 1371.0e7_dp
    case('CH4')
      L_heat = 480.6e7_dp
    case('NH4SH')
      L_heat = 2409.4_dp * log(10.0_dp) * R_gas / mol_w_sp
    case('H2S')
      L_heat = 958.587_dp * log(10.0_dp) * R_gas / mol_w_sp
    case('S2')
      L_heat = 14000.0_dp * R_gas / mol_w_sp
    case('S8')
      L_heat = 7510.0_dp * R_gas / mol_w_sp
    case('CO')
      L_heat = 7.8824e2_dp * log(10.0_dp) * R_gas / mol_w_sp
    case('CO2')
      L_heat = 1.5119e3_dp * log(10.0_dp) * R_gas / mol_w_sp
    case('H2SO4')
      L_heat = 1.01294e4_dp * R_gas / mol_w_sp
    case default
      print*, 'Latent heat: dust species not found: ', trim(sp), 'Stopping!'
      stop
    end select

  end function l_heat_sp

  !! Dummy subroutine required for dopri5
  subroutine solout(nr,xold,x,y,n,con,icomp,nd,rpar,ipar,irtrn)
    implicit none
    integer, intent(in) :: nr, n, nd, irtrn, ipar(*)
    real(dp), intent(in) ::  Y(*),CON(*),ICOMP(*),RPAR(*), xold, x
  end subroutine solout

end module mini_cloud_1_s_mod

