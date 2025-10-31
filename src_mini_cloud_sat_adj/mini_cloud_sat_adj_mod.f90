module mini_cloud_sat_adj_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: bar = 1.0e6_dp ! bar to dyne
  real(dp), parameter :: atm = 1.01325e6_dp ! atm to dyne
  real(dp), parameter :: pa = 10.0_dp ! pa to dyne
  real(dp), parameter :: mmHg = 1333.22387415_dp  ! mmHg to dyne

  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp) ! value of pi
  real(dp), parameter :: twopi = pi * 2.0_dp

  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: amu = 1.66053906660e-24_dp ! g mol-1 (note, has to be cgs g mol-1 !!!)
  real(dp), parameter :: R = 8.31446261815324e7_dp


  public :: mini_cloud_sat_adj, dqdt
  private :: p_vap_sp! , l_heat_sp

contains 

  subroutine mini_cloud_sat_adj(i, t_step, mol_w_sp, sp, rho_d, & 
    & T_in, P_in, rho_in, met, tau_cond, q_v, q_c)
    implicit none

    integer, intent(in) :: i
    character(len=20), intent(in) :: sp
    real(dp), intent(in) :: t_step
    real(dp), intent(in) :: mol_w_sp, rho_d, T_in, P_in, rho_in, met, tau_cond
    
    real(dp), intent(inout) :: q_v, q_c

    real(dp) :: t_now

    integer :: itol, iout, idid
    real(dp) :: rtol, atol

    integer, parameter :: n = 2
    integer, parameter :: lwork = 8*n + 21 !8*N+5*NRDENS+21
    integer, parameter :: liwork = 21
    real(dp), dimension(lwork) :: work
    integer, dimension(liwork) :: iwork

    !! Work variables
    real(dp) :: Tl, pl, rho, Rd_v, p_vap, q_s
    real(dp), dimension(n) :: y
    integer :: ipar
    real(dp), dimension(2) :: rpar

    !! Convert inputs to cgs units
    Tl = T_in
    pl = P_in * 10.0_dp ! Pa to dyne
    rho = rho_in

    !! Vapour pressure and saturation vapour pressure
    p_vap = p_vap_sp(sp, Tl, pl, met)

    !! Equilibrium (sat = 1) vapour mass density rho_s [g cm^-3]
    Rd_v = R/mol_w_sp
    q_s = (p_vap/(Rd_v * Tl))/rho

    !! Start time integration
    t_now = 0.0_dp

    !! Parameters for dopri5 Runge-Kutta
    itol = 0
    rtol = 1e-3_dp
    atol = 1e-30_dp

    iout = 0

    !! Set ipar and rpar to send through variables into dopri5 to calculate dqdt
    ipar = 0
    rpar(1) = q_s
    rpar(2) = tau_cond

    work(:) = 0.0_dp
    iwork(:) = 0

    !! initial conditions - convert to VMR from MMR
    y(1) = max(q_v,1e-30_dp)
    y(2) = max(q_c,1e-30_dp)

    call dopri5(n,dqdt,t_now,y,t_step, &
      &         rtol,atol,itol, &
      &         solout,iout, &
      &         work,lwork,iwork,liwork,rpar,ipar,idid)

    !! Final results, convert back to MMR
    q_v = max(y(1),1e-30_dp)
    q_c = max(y(2),1e-30_dp)

    if (idid /= 1) then
      print*, 'error in dopri5: ', idid
    end if

  end subroutine mini_cloud_sat_adj

  !! Calculate tracer fluxes
  subroutine dqdt(n,x,y,f,rpar,ipar)
    implicit none

    integer, intent(in) :: n
    integer, intent(in) :: ipar

    real(dp), intent(in) :: x 
    real(dp), dimension(n), intent(inout) :: y
    real(dp), dimension(2), intent(in) :: rpar

    real(dp), dimension(n), intent(out) :: f

    real(dp) :: q_v, q_c, q_s, tau_cond

    !! y(1) = q_v
    !! y(2) = q_c
    !! rpar(1) = q_s
    !! rpar(2) = tau_cond

    q_v = max(y(1),1e-30_dp)
    q_c = max(y(2),1e-30_dp)
    q_s = rpar(1)
    tau_cond = rpar(2)

    if (q_v < q_s) then
      !! Vapour is undersaturated - adjust vapour by evaporating from q_c
      f(1) = min(q_s - q_v, q_c)/tau_cond
    else if (q_v > q_s) then
      !! Vapour is supersaturated - adjust vapour by condensing from q_v
      f(1) = -(q_v - q_s)/tau_cond
    else 
      !! Do nothing as q_v = q_s
      f(1) = 0.0_dp
    end if

    !! RHS of condensate is negative vapour RHS
    f(2) = -f(1)

  end subroutine dqdt

  !! Vapour pressure for each species
  function p_vap_sp(sp, T, p, met) result(p_vap)
    implicit none

    character(len=20), intent(in) :: sp
    real(dp), intent(in) :: T, p, met

    real(dp) :: TC, f
    !real(dp) :: A, B, C
    real(dp) :: p_vap

    ! Return vapour pressure in dyne
    select case(trim(sp))
    case('C')
      ! Gail & Sedlmayr (2013) - I think...
      p_vap = exp(3.27860e1_dp - 8.65139e4_dp/(T + 4.80395e-1_dp))
    case('TiC')
      ! Kimura et al. (2023)
      p_vap = 10.0_dp**(-33600.0_dp/T + 7.652_dp) * atm
    case('SiC')
      ! Elspeth 5 polynomial JANAF-NIST fit
      p_vap =  exp(-9.51431385e4_dp/T + 3.72019157e1_dp + 1.09809718e-3_dp*T &
        & -5.63629542e-7_dp*T**2 + 6.97886017e-11_dp*T**3)
    case('CaTiO3')
      ! Wakeford et al. (2017) -  taken from VIRGA
      p_vap = 10.0_dp**(-72160.0_dp/T + 30.24_dp - log10(p/1e6_dp) - 2.0_dp*met) * bar
      ! Kozasa et al. (1987)
      !p_vap = exp(-79568.2_dp/T + 42.0204_dp) * atm
    case('TiO2')
      ! GGChem 5 polynomial NIST fit
      p_vap = exp(-7.70443e4_dp/T +  4.03144e1_dp - 2.59140e-3_dp*T &
        &  + 6.02422e-7_dp*T**2 - 6.86899e-11_dp*T**3)
    case('VO')
      ! NIST 5 param fit
      p_vap = exp(-6.74603e4_dp/T + 3.82717e1_dp - 2.78551e-3_dp*T &
        & + 5.72078e-7_dp*T**2 - 7.41840e-11_dp*T**3)
    case('Al2O3')
      ! Wakeford et al. (2017) - taken from CARMA
      p_vap = 10.0_dp**(17.7_dp - 45892.6_dp/T - 1.66_dp*met) * bar
      ! Kozasa et al. (1987)
      !p_vap = exp(-73503.0_dp/T + 22.005_dp) * atm
    case('Fe')
      ! Visscher et al. (2010) - taken from CARMA
      p_vap = 10.0_dp**(7.23_dp - 20995.0_dp/T) * bar
      ! Elspeth note: Changed to Ackerman & Marley et al. (2001) expression
      ! if (T > 1800.0_dp) then
      !   p_vap = exp(9.86_dp - 37120.0_dp/T) * bar
      ! else
      !   p_vap = exp(15.71_dp - 47664.0_dp/T) * bar
      ! end if
    case('FeS')
      ! GGChem 5 polynomial NIST fit
      p_vap = exp(-5.69922e4_dp/T + 3.86753e1_dp - 4.68301e-3_dp*T &
        & + 1.03559e-6_dp*T**2 - 8.42872e-11_dp*T**3)
    case('FeO')
      ! GGChem 5 polynomial NIST fit
      p_vap = exp(-6.30018e4_dp/T + 3.66364e1_dp - 2.42990e-3_dp*T &
        & + 3.18636e-7_dp*T**2)
    case('Mg2SiO4')
      ! Visscher et al. (2010)/Visscher notes - taken from CARMA
      p_vap = 10.0_dp**(14.88_dp - 32488.0_dp/T - 1.4_dp*met - 0.2_dp*log10(p/1e6_dp)) * bar
      ! Kozasa et al. (1989) - Seems to be too high
      !p_vap = p*10.0**(8.25_dp -  27250.0_dp/T - log10(p/1e6_dp) + 3.58_dp)
    case('MgSiO3','MgSiO3_amorph')
      ! Visscher - taken from VIRGA
      p_vap = 10.0_dp**(13.43_dp - 28665.0_dp/T - met) * bar
      ! Ackerman & Marley (2001)
      !p_vap = exp(-58663.0_dp/T + 25.37_dp) * bar
    case('MgO')
      ! GGChem 5 polynomial NIST fit
      p_vap = exp(-7.91838e4_dp/T + 3.57312e1_dp + 1.45021e-4_dp*T &
        &  - 8.47194e-8*T**2 + 4.49221e-12_dp*T**3)
    case('SiO2','SiO2_amorph')
      ! GGChem 5 polynomial NIST fit
      p_vap = exp(-7.28086e4_dp/T + 3.65312e1_dp - 2.56109e-4_dp*T &
        & - 5.24980e-7_dp*T**2 + 1.53343E-10_dp*T**3) 
    case('SiO')
      ! Gail et al. (2013)
      p_vap = exp(-49520.0_dp/T + 32.52_dp)
    case('Cr')
      ! GGChem 5 polynomial NIST fit
      p_vap = exp(-4.78455e+4_dp/T + 3.22423e1_dp - 5.28710e-4_dp*T & 
        &  - 6.17347e-8_dp*T**2 + 2.88469e-12_dp*T**3)
     case('MnS')
      ! Morley et al. (2012)
      p_vap = 10.0_dp**(11.532_dp - 23810.0_dp/T - met) * bar
    case('Na2S')
      ! Morley et al. (2012)
      p_vap =  10.0_dp**(8.550_dp - 13889.0_dp/T - 0.5_dp*met) * bar
    case('ZnS')
      ! Elspeth 5 polynomial Barin data fit
      p_vap = exp(-4.75507888e4_dp/T + 3.66993865e1_dp - 2.49490016e-3_dp*T &
        &  + 7.29116854e-7_dp*T**2 - 1.12734453e-10_dp*T**3)
      ! Morley et al. (2012)
      !p_vap = 10.0_dp**(12.812_dp - 15873.0_dp/T - met) * bar        
    case('KCl')
      ! GGChem 5 polynomial NIST fit
      p_vap = exp(-2.69250e4_dp/T + 3.39574e+1_dp - 2.04903e-3_dp*T &
        & -2.83957e-7_dp*T**2 + 1.82974e-10_dp*T**3)
      ! Morley et al. (2012)
      !p_vap =  10.0_dp**(7.611_dp - 11382.0_dp/T) * bar
    case('NaCl')
      ! GGChem 5 polynomial NIST fit
      p_vap = exp(-2.79146e4_dp/T + 3.46023e1_dp - 3.11287e3_dp*T & 
        & + 5.30965e-7_dp*T**2 -2.59584e-12_dp*T**3)
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
    case('NH4Cl')
      ! Unknown - I think I fit this?
      p_vap = 10.0_dp**(7.0220_dp - 4302.0_dp/T) * bar
    case('H2O')
      TC = T - 273.15_dp
      ! Huang (2018) - A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice
      if (TC < 0.0_dp) then
        f = 0.99882_dp * exp(0.00000008_dp * p/pa)
        p_vap = exp(43.494_dp - (6545.8_dp/(TC + 278.0_dp)))/(TC + 868.0_dp)**2.0_dp * pa * f
      else
        f = 1.00071_dp * exp(0.000000045_dp * p/pa)
        p_vap = exp(34.494_dp - (4924.99_dp/(TC + 237.1_dp)))/(TC + 105.0_dp)**1.57_dp * pa * f
      end if
      ! Ackerman & Marley (2001) H2O liquid & ice vapour pressure expressions
      !if (T > 1048.0_dp) then
      !  p_vap = 6.0e8_dp
      !else if (T < 273.16_dp) then
      !  p_vap = 6111.5_dp * exp((23.036_dp * TC - TC**2/333.7_dp)/(TC + 279.82_dp))
      !else
      !  p_vap = 6112.1_dp * exp((18.729_dp * TC - TC**2/227.3_dp)/(TC + 257.87_dp)) 
      !end if
    case('NH3')
      ! Blakley et al. (2024) - experimental to low T and pressure
      p_vap = exp(-5.55_dp - 3605.0_dp/T + 4.82792_dp*log(T) - 0.024895_dp*T + 2.1669e-5_dp*T**2 - 2.3575e-8_dp *T**3) * bar
      ! Ackerman & Marley (2001) NH3 ice vapour pressure expression fit from Weast (1971) data
      !p_vap = exp(10.53_dp - 2161.0_dp/T - 86596.0_dp/T**2)  * bar
    case('CH4')
      ! Frey & Schmitt (2009)
      p_vap = exp(1.051e1_dp - 1.110e3_dp/T - 4.341e3_dp/T**2 + 1.035e5_dp/T**3 - 7.910e5_dp/T**4) * bar
      ! Lodders & Fegley (1998) - directly taken from VIRGA
      ! if (T < 90.68_dp) then
      !   C = -16.043_dp/8.3143_dp * (2.213_dp - 2.650_dp)
      !   B = -16.043_dp/8.3143_dp * (611.10_dp + (2.213_dp - 2.650_dp) * 90.68_dp)
      !   A = 0.11719_dp * 90.68_dp**(-C) * exp(-B/90.68_dp)
      ! else
      !   C = -16.043_dp/8.3143_dp * (2.213_dp - 3.370_dp)
      !   B = -16.043_dp/8.3143_dp * (552.36_dp + (2.213_dp - 3.370_dp) * 90.68_dp)
      !   A = 0.11719_dp * 90.68_dp**(-C) * exp(-B/90.68_dp)
      ! end if
      ! p_vap = A * T**C * exp(B/T) * bar

    case('NH4SH')
      !--- E.Lee's fit to Walker & Lumsden (1897) ---
      p_vap = 10.0_dp**(7.8974_dp - 2409.4_dp/T) * bar
    case('H2S')
      ! Frey & Schmitt (2009)
      p_vap = exp(12.98_dp - 2.707e3_dp/T) * bar
    case('H2SO4')
      ! GGChem 5 polynomial NIST fit
      p_vap = exp(-1.01294e4_dp/T + 3.55465e1_dp - 8.34848e-3_dp*T)      
    case('CO')
      ! Frey & Schmitt (2009)
      if (T < 61.55_dp) then
        p_vap = exp(1.043e1_dp - 7.213e2_dp/T - 1.074e4_dp/T**2 + 2.341e5_dp/T**3 - 2.392e6_dp/T**4 + 9.478e6_dp/T**5) * bar
      else
        p_vap = exp(1.025e1_dp - 7.482e2_dp/T - 5.843e3_dp/T**2 + 3.939e4_dp/T**3) * bar
      end if
      ! Yaws
      !p_vap = 10.0_dp**(51.8145e0_dp - 7.8824e2_dp/T - 2.2734e1_dp*log10(T) &
      !  & + 5.1225e-2_dp*T + 4.6603e-11_dp*T**2) * mmHg
    case('CO2')
      ! Frey & Schmitt (2009)
      if (T < 194.7_dp) then
        p_vap = exp(1.476e1_dp - 2.571e3_dp/T - 7.781e4_dp/T**2 + 4.325e6_dp/T**3 - 1.207e8_dp/T**4 + 1.350e9_dp/T**5) * bar
      else
        p_vap = exp(1.861e1_dp - 4.154e3_dp/T + 1.041e5_dp/T**2) * bar
      end if      
      ! Yaws
      !p_vap = 10.0_dp**(35.0187e0_dp - 1.5119e3_dp/T - 1.1335e1_dp*log10(T) &
      !  & + 9.3383e-3_dp*T + 7.7626e-10_dp*T**2) * mmHg
    case('O2')
      ! Blakley et al. (2024) - experimental to low T and pressure (beta O2)
      p_vap = exp(15.29_dp - 1166.2_dp/T - 0.75587_dp*log(T) + 0.14188_dp*T - 1.8665e-3_dp*T**2 + 7.582e-6_dp *T**3) * bar
    case default
      print*, 'Vapour pressure species not found: ', trim(sp)
      print*, 'STOP'
      stop
    end select

  end function p_vap_sp

  !! latent heat for each species
  !function l_heat_sp(sp, T) result(L_heat)
  !   implicit none

  !   character(len=20), intent(in) :: sp
  !   real(dp), intent(in) :: T

  !   real(dp) :: L_heat


  !   ! Return latent heat in erg g-1
  !   ! Get value from vapour pressure expression (or special function)
  !   select case(trim(sp))
  !   case('C')
  !     L_heat = 41523.0_dp * log(10.0_dp) * R_gas / mol_w_sp
  !   case('TiC')
  !     L_heat = 33600.0_dp * log(10.0_dp) * R_gas / mol_w_sp
  !   case('SiC')
  !     L_heat = 9.51431385e4_dp * R_gas / mol_w_sp
  !   case('CaTiO3')
  !     L_heat = 79568.2_dp * R_gas / mol_w_sp
  !   case('Al2O3')
  !     L_heat = 73503.0_dp * R_gas / mol_w_sp
  !   case('TiO2')
  !     L_heat = 7.70443e4_dp * R_gas / mol_w_sp
  !   case('VO')
  !     L_heat = 6.74603e4_dp * R_gas / mol_w_sp
  !   case('Fe')
  !     L_heat = 37120.0_dp * R_gas / mol_w_sp
  !   case('FeS')
  !     L_heat = 5.69922e4_dp * R_gas / mol_w_sp
  !   case('FeO')
  !     L_heat = 6.30018e4_dp * R_gas / mol_w_sp
  !   case('MgS')
  !     L_heat = 5.92010440e4_dp * R_gas / mol_w_sp
  !   case('Mg2SiO4')
  !     L_heat = 62279.0_dp * R_gas / mol_w_sp
  !   case('MgSiO3','MgSiO3_amorph')
  !     L_heat = 58663.0_dp * R_gas / mol_w_sp
  !   case('MgO')
  !     L_heat = 7.91838e4_dp * R_gas / mol_w_sp
  !   case('SiO2','SiO2_amorph')
  !     L_heat = 7.28086e4_dp * R_gas / mol_w_sp
  !   case('SiO')
  !     L_heat = 49520.0_dp * R_gas / mol_w_sp
  !   case('Cr')
  !     L_heat = 4.78455e4_dp * R_gas / mol_w_sp
  !   case('MnS')
  !     L_heat = 23810.0_dp * log(10.0_dp) * R_gas / mol_w_sp
  !   case('Na2S')
  !     L_heat = 13889.0_dp * log(10.0_dp) * R_gas / mol_w_sp
  !   case('ZnS')
  !     L_heat = 4.75507888e4_dp * R_gas / mol_w_sp
  !   case('KCl')
  !     L_heat = 2.69250e4_dp * R_gas / mol_w_sp
  !   case('NaCl')
  !     L_heat = 2.79146e4_dp * R_gas / mol_w_sp
  !   case('NH4Cl')
  !     L_heat = 4302.0_dp * log(10.0_dp) * R_gas / mol_w_sp
  !   case('H2O')
  !     L_heat = 2257.0e7_dp
  !   case('NH3')
  !     L_heat = 1371.0e7_dp
  !   case('CH4')
  !     L_heat = 480.6e7_dp
  !   case('NH4SH')
  !     L_heat = 2409.4_dp * log(10.0_dp) * R_gas / mol_w_sp
  !   case('H2S')
  !     L_heat = 958.587_dp * log(10.0_dp) * R_gas / mol_w_sp
  !   case('S2')
  !     L_heat = 14000.0_dp * R_gas / mol_w_sp
  !   case('S8')
  !     L_heat = 7510.0_dp * R_gas / mol_w_sp
  !   case('CO')
  !     L_heat = 7.8824e2_dp * log(10.0_dp) * R_gas / mol_w_sp
  !   case('CO2')
  !     L_heat = 1.5119e3_dp * log(10.0_dp) * R_gas / mol_w_sp
  !   case('H2SO4')
  !     L_heat = 1.01294e4_dp * R_gas / mol_w_sp
  !   case default
  !     print*, 'Latent heat: dust species not found: ', trim(sp), 'Stopping!'
  !     stop
  !   end select

  ! end function l_heat_sp

  !! Dummy subroutine required for dopri5
  subroutine solout(nr,xold,x,y,n,con,icomp,nd,rpar,ipar,irtrn)
    implicit none
    integer, intent(in) :: nr, n, nd, irtrn, ipar(*)
    real(dp), intent(in) ::  Y(*),CON(*),ICOMP(*),RPAR(*), xold, x
  end subroutine solout

end module mini_cloud_sat_adj_mod
