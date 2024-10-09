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

  real(dp) :: pl, p_vap, q_s, tau_deep

  real(dp), parameter :: bar = 1.0e6_dp ! bar to dyne
  real(dp), parameter :: atm = 1.01325e6_dp ! atm to dyne
  real(dp), parameter :: pa = 10.0_dp ! pa to dyne
  real(dp), parameter :: mmHg = 1333.22387415_dp  ! mmHg to dyne

  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp) ! value of pi
  real(dp), parameter :: twopi = pi * 2.0_dp

  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: amu = 1.66053906660e-24_dp ! g mol-1 (note, has to be cgs g mol-1 !!!)
  real(dp), parameter :: R_gas = 8.31446261815324e7_dp


  public :: mini_cloud_1_s, dqdt
  private :: p_vap_sp

contains 

  subroutine mini_cloud_1_s(deep_flag, T_in, P_in, grav_in, mu_in, t_end, sp, q_v, q_c)
    implicit none

    logical, intent(in) :: deep_flag
    character(len=20), intent(in) :: sp
    real(dp), intent(in) :: t_end
    real(dp), intent(in) :: T_in, P_in, grav_in, mu_in
    
    real(dp), intent(inout) :: q_v, q_c

    real(dp) :: Tl, mu, grav
    real(dp) :: Hp, eps

    real(dp) :: t_now

    integer :: itol, iout, idid
    integer :: ipar
    real(dp) :: rtol, atol
    real(dp) :: rpar

    integer, parameter :: n = 2
    integer, parameter :: lwork = 37 !8*2
    integer, parameter :: liwork = 21
    real(dp), dimension(lwork) :: work
    integer, dimension(liwork) :: iwork

    real(dp), dimension(n) :: y

    !! Convert inputs to cgs units
    Tl = T_in
    pl = P_in * 10.0_dp
    mu = mu_in 
    grav = grav_in * 100.0_dp

    !! Calculate deep replenishment rate of vapour
    Hp = (kb * Tl) / (mu * amu * grav)
    tau_deep = Hp**2/Kzz_deep

    !! Conversion factor between MMR and VMR (VMR -> MMR) for cloud species
    eps = mol_w_sp/mu

    !! Vapour pressure and saturation vapour pressure
    p_vap = p_vap_sp(sp, Tl)

    !! Equilibrium (sat = 1) vapour pressure fraction
    q_s = min(1.0_dp, p_vap/pl)

    !! Start time integration
    t_now = 0.0_dp

    !! Parameters for dopri5 Runge-Kutta
    itol = 0
    rtol = 1e-3_dp
    atol = 1e-30_dp

    iout = 0

    !! Set ipar and rpar to send through variables into dopri5 to calculate dqdt
    ipar = 0
    if (deep_flag .eqv. .True.) then
      ipar = 1
    end if
    rpar = 0.0_dp


    work(:) = 0.0_dp
    iwork(:) = 0

    !! initial conditions - convert to VMR from MMR
    y(1) = q_v/eps
    y(2) = q_c/eps

    call dopri5(n,dqdt,t_now,y,t_end, &
      &         rtol,atol,itol, &
      &         solout,iout, &
      &         work,lwork,iwork,liwork,rpar,ipar,idid)

    !! Final results, convert back to MMR
    q_v = y(1) * eps
    q_c = y(2) * eps

    if (idid /= 1) then
      print*, 'error in dopri5: ', idid
    end if

  end subroutine mini_cloud_1_s

  !! Calculate tracer fluxes
  subroutine dqdt(n,x,y,f,rpar,ipar)
    implicit none

    integer, intent(in) :: n
    integer, intent(in) :: ipar

    real(dp), intent(in) :: x 
    real(dp), dimension(n), intent(inout) :: y
    real(dp), intent(in) :: rpar

    real(dp), dimension(n), intent(out) :: f

    real(dp) :: sat

    !y(1) = q_v, y(2) = q_c

    y(1) = max(y(1),1e-30_dp)
    y(2) = max(y(2),1e-30_dp)

    !! Calculate dqdt for vapour and condensate
    sat = max((y(1) * pl)/p_vap, 1e-99_dp)

    !! Calculate dqdt given the supersaturation ratio
    if (sat < 0.99_dp) then

      ! Evaporate from q_c
      f(1) = min(q_s - y(1), y(2))/tau_chem

    else if (sat > 1.01_dp) then

      ! Condense q_v toward the saturation ratio
      f(1) = -(y(1) - q_s)/tau_chem

    else
      f(1) = 0.0_dp
    end if

    f(2) = -f(1)

    ! Add replenishment to lower boundary at the tau_deep rate
    if (ipar == 1) then
      f(1) = f(1) - (y(1) - q_v_deep)/tau_deep
    end if

  end subroutine dqdt

  !! Dummy subroutine required for dopri5
  subroutine solout(nr,xold,x,y,n,con,icomp,nd,rpar,ipar,irtrn)
    implicit none
    integer, intent(in) :: nr, n, nd, irtrn, ipar(*)
    real(dp), intent(in) ::  Y(*),CON(*),ICOMP(*),RPAR(*), xold, x
  end subroutine solout

  !! Vapour pressure for each species
  real(dp) function p_vap_sp(sp, T)
    implicit none

    character(len=20), intent(in) :: sp
    real(dp), intent(in) :: T

    real(dp) :: TC

    ! Return vapour pressure in dyne
    select case(sp)
    case('C')
      ! Kimura et al. (2023)
      p_vap_sp = 10.0_dp**(-41523.0_dp/T + 10.609_dp) * atm
    case('TiC')
      ! Kimura et al. (2023)
      p_vap_sp = 10.0_dp**(-33600.0_dp/T + 7.652_dp) * atm
    case('SiC')
      ! Elspeth 5 polynomial JANAF-NIST fit
      p_vap_sp =  exp(-9.51431385e4_dp/T + 3.72019157e1_dp + 1.09809718e-3_dp*T &
        & -5.63629542e-7_dp*T**2 + 6.97886017e-11_dp*T**3)
    case('CaTiO3')
      ! Kozasa et al. (1987)
      p_vap_sp = exp(-79568.2_dp/T + 42.0204_dp) * atm  
    case('Al2O3')
      ! Kozasa et al. (1989)
      p_vap_sp = exp(-73503.0_dp/T + 22.01_dp) * atm      
    case('TiO2')
      ! GGChem 5 polynomial NIST fit
      p_vap_sp = exp(-7.70443e4_dp/T +  4.03144e1_dp - 2.59140e-3_dp*T &
        &  + 6.02422e-7_dp*T**2 - 6.86899e-11_dp*T**3)
    case('VO')
      ! GGChem 5 polynomial NIST fit
      p_vap_sp =  exp(-6.74603e4_dp/T + 3.82717e1_dp - 2.78551e-3_dp*T &
        & + 5.72078e-7_dp*T**2 - 7.41840e-11_dp*T**3)         
    case('Fe')
      ! Elspeth note: Changed to Ackerman & Marley et al. (2001) expression
      if (T > 1800.0_dp) then
        p_vap_sp = exp(9.86_dp - 37120.0_dp/T) * bar
      else
        p_vap_sp = exp(15.71_dp - 47664.0_dp/T) * bar
      end if
    case('FeS')
      ! GGChem 5 polynomial NIST fit
      p_vap_sp = exp(-5.69922e4_dp/T + 3.86753e1_dp - 4.68301e-3_dp*T &
        & + 1.03559e-6_dp*T**2 - 8.42872e-11_dp*T**3)
    case('FeO')
      ! GGChem 5 polynomial NIST fit
      p_vap_sp = exp(-6.30018e4_dp/T + 3.66364e1_dp - 2.42990e-3_dp*T &
        & + 3.18636e-7_dp*T**2)     
    case('MgS')
      ! Elspeth 5 polynomial JANAF-NIST fit
      p_vap_sp = exp(-5.92010440e4_dp/T +  3.58148138e1_dp - 1.93822353e-3_dp*T &
        &  + 9.77971406e-7_dp*T**2 - 11.74813601e-10_dp*T**3)   
    case('Mg2SiO4')
      ! Kozasa et al. (1989)
      p_vap_sp = exp(-62279_dp/T + 20.944_dp) * atm
    case('MgSiO3','MgSiO3_amorph')
      ! Ackerman & Marley (2001)
      p_vap_sp = exp(-58663.0_dp/T + 25.37_dp) * bar
    case('SiO2','SiO2_amorph')
      ! NIST 5 param fit
      p_vap_sp = exp(-7.28086e4_dp/T + 3.65312e1_dp - 2.56109e-4_dp*T &
      & - 5.24980e-7_dp*T**2 + 1.53343E-10_dp*T**3) * bar
    case('SiO')
      ! Gail et al. (2013)
      p_vap_sp = exp(-49520.0_dp/T + 32.52_dp)
    case('Cr')
      ! GGChem 5 polynomial NIST fit
      p_vap_sp = exp(-4.78455e4_dp/T + 3.22423e1_dp - 5.28710e-4_dp*T & 
        &  - 6.17347e-8_dp*T**2 + 2.88469e-12_dp*T**3)
    case('MnS')
      ! Morley et al. (2012)
      p_vap_sp = 10.0_dp**(11.532_dp - 23810.0_dp/T) * bar
    case('Na2S')
      ! Morley et al. (2012)
      p_vap_sp =  10.0_dp**(8.550_dp - 13889.0_dp/T) * bar
    case('ZnS')
      ! Elspeth 5 polynomial Barin data fit
      p_vap_sp = exp(-4.75507888e4_dp/T + 3.66993865e1_dp - 2.49490016e-3_dp*T &
        &  + 7.29116854e-7_dp*T**2 - 1.12734453e-10_dp*T**3)
    case('KCl')
      ! GGChem 5 polynomial NIST fit
      p_vap_sp = exp(-2.69250e4_dp/T + 3.39574e1_dp - 2.04903e-3_dp*T &
          & -2.83957e-7_dp*T**2 + 1.82974e-10_dp*T**3)
    case('NaCl')
      ! GGChem 5 polynomial NIST fit
      p_vap_sp = exp(-2.79146e4_dp/T + 3.46023e1_dp - 3.11287e3_dp*T & 
        & + 5.30965e-7_dp*T**2 -2.59584e-12_dp*T**3)
    case('NH4Cl')
      p_vap_sp = 10.0_dp**(7.0220_dp - 4302.0_dp/T) * bar
    case('H2O')
      ! Ackerman & Marley (2001) H2O liquid & ice vapour pressure expressions
      TC = T - 273.15_dp
      if (T > 1048.0_dp) then
        p_vap_sp = 6.0e8_dp
      else if (T < 273.16_dp) then
        p_vap_sp = 6111.5_dp * exp((23.036_dp * TC - TC**2/333.7_dp)/(TC + 279.82_dp))
      else
        p_vap_sp = 6112.1_dp * exp((18.729_dp * TC - TC**2/227.3_dp)/(TC + 257.87_dp)) 
      end if
    case('NH3')
      ! Ackerman & Marley (2001) NH3 ice vapour pressure expression from Weast (1971)
      p_vap_sp = exp(10.53_dp - 2161.0_dp/T - 86596.0_dp/T**2)  * bar
    case('CH4')
      !--- Prydz, R.; Goodwin, R.D., J. Chem. Thermodyn., 1972, 4,1 ---
      if (T < 0.5_dp) then ! Limiter for very very very cold T
        p_vap_sp = 10.0_dp**(3.9895_dp - 443.028_dp/(0.5_dp-0.49_dp)) * bar       
      else
        p_vap_sp = 10.0_dp**(3.9895_dp - 443.028_dp/(T-0.49_dp)) * bar
      end if
    case('NH4SH')
      !--- E.Lee's fit to Walker & Lumsden (1897) ---
      p_vap_sp = 10.0_dp**(7.8974_dp - 2409.4_dp/T) * bar
    case('H2S')
      !--- Stull (1947) ---
      if (T < 30.0_dp) then ! Limiter for very cold T
        p_vap_sp = 10.0_dp**(4.43681_dp - 829.439_dp/(30.0_dp-25.412_dp)) * bar
      else if (T < 212.8_dp) then
        p_vap_sp = 10.0_dp**(4.43681_dp - 829.439_dp/(T-25.412_dp)) * bar
      else
        p_vap_sp = 10.0_dp**(4.52887_dp - 958.587_dp/(T-0.539_dp)) * bar
      end if
    case('S2')
      !--- Zahnle et al. (2016) ---
      if (T < 413.0_dp) then
        p_vap_sp = exp(27.0_dp - 18500.0_dp/T) * bar
      else
        p_vap_sp = exp(16.1_dp - 14000.0_dp/T) * bar
      end if
    case('S8')
      !--- Zahnle et al. (2016) ---
      if (T < 413.0_dp) then
        p_vap_sp = exp(20.0_dp - 11800.0_dp/T) * bar
      else
        p_vap_sp = exp(9.6_dp - 7510.0_dp/T) * bar
      end if
    case('CO')
      ! Yaws
      p_vap_sp = 10.0_dp**(51.8145e0_dp - 7.8824e2_dp/T - 2.2734e1_dp*log10(T) &
        & + 5.1225e-2_dp*T + 4.6603e-11_dp*T**2) * mmHg
    case('CO2')
      ! Yaws
      p_vap_sp = 10.0_dp**(35.0187e0_dp - 1.5119e3_dp/T - 1.1335e1_dp*log10(T) &
          & + 9.3383e-3_dp*T + 7.7626e-10_dp*T**2) * mmHg
    case('H2SO4')
      ! GGChem 5 polynomial NIST fit
      p_vap_sp = exp(-1.01294e4_dp/T + 3.55465e1_dp - 8.34848e-3_dp*T)           
    case default
      print*, 'Saturation: dust species not found: ', trim(sp), 'Stopping!'
      stop
    end select

  end function p_vap_sp

end module mini_cloud_1_s_mod

