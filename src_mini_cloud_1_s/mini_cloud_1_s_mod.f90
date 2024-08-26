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

  real(dp), parameter :: bar = 1.0e6_dp ! bar to dyne
  real(dp), parameter :: atm = 1.01325e6_dp ! atm to dyne
  real(dp), parameter :: pa = 10.0_dp ! pa to dyne

  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp) ! value of pi
  real(dp), parameter :: twopi = pi * 2.0_dp

  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: amu = 1.66053906660e-24_dp ! g mol-1 (note, has to be cgs g mol-1 !!!)
  real(dp), parameter :: R_gas = 8.31446261815324e7_dp

  !! Diameter, LJ potential and molecular weight for background gases ! Do everything in cgs for vf calculation
  real(dp), parameter :: d_OH = 3.06e-8_dp, LJ_OH = 100.0_dp * kb, molg_OH = 17.00734_dp  ! estimate
  real(dp), parameter :: d_H2 = 2.827e-8_dp, LJ_H2 = 59.7_dp * kb, molg_H2 = 2.01588_dp
  real(dp), parameter :: d_H2O = 2.641e-8_dp, LJ_H2O = 809.1_dp * kb, molg_H2O = 18.01528_dp
  real(dp), parameter :: d_H = 2.5e-8_dp, LJ_H =  30.0_dp * kb, molg_H = 1.00794_dp
  real(dp), parameter :: d_CO = 3.690e-8_dp, LJ_CO = 91.7_dp * kb, molg_CO = 28.0101_dp
  real(dp), parameter :: d_CO2 = 3.941e-8_dp, LJ_CO2 = 195.2_dp * kb, molg_CO2 = 44.0095_dp
  real(dp), parameter :: d_O = 2.66e-8_dp, LJ_O = 70.0_dp * kb, molg_O = 15.99940_dp
  real(dp), parameter :: d_CH4 = 3.758e-8_dp, LJ_CH4 = 148.6_dp * kb, molg_CH4 = 16.0425_dp
  real(dp), parameter :: d_C2H2 = 4.033e-8_dp, LJ_C2H2 = 231.8_dp * kb, molg_C2H2 = 26.0373_dp
  real(dp), parameter :: d_NH3 = 2.900e-8_dp, LJ_NH3 = 558.3_dp * kb, molg_NH3 = 17.03052_dp
  real(dp), parameter :: d_N2 = 3.798e-8_dp, LJ_N2 = 71.4_dp * kb, molg_N2 = 14.0067_dp
  real(dp), parameter :: d_HCN = 3.630e-8_dp, LJ_HCN = 569.1_dp * kb, molg_HCN = 27.0253_dp
  real(dp), parameter :: d_He = 2.511e-8_dp, LJ_He = 10.22_dp * kb, molg_He = 4.002602_dp

  !! Constuct required arrays for calculating gas mixtures
  real(dp), allocatable, dimension(:) :: d_g, LJ_g, molg_g

  public :: mini_cloud_1_s, dqdt
  private :: p_vap_sp, eta_construct

contains 

  subroutine mini_cloud_1_s(deep_flag, T_in, P_in, grav_in, mu_in, bg_VMR_in, t_end, sp, sp_bg, q_v, q_c, v_f)
    implicit none

    logical, intent(in) :: deep_flag
    character(len=20), intent(in) :: sp
    character(len=20), dimension(:), intent(in) :: sp_bg
    real(dp), intent(in) :: t_end
    real(dp), intent(in) :: T_in, P_in, grav_in, mu_in
    real(dp), dimension(:), intent(in)::  bg_VMR_in
    
    real(dp), intent(inout) :: q_v, q_c, v_f

    integer :: k, n_gas
    real(dp) :: Tl, pl, mu, grav
    real(dp) :: p_vap, q_s, tau_deep, Hp, rV, eps
    real(dp), allocatable, dimension(:) :: VMR_g

    real(dp) :: top, bot, eta, mfp, beta, rho, eta_g, Kn

    real(dp) :: t_now

    integer :: itol, iout, idid
    integer :: ipar
    real(dp) :: rtol, atol
    real(dp), dimension(5) :: rpar

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

    n_gas = size(bg_VMR_in)
    allocate(VMR_g(n_gas))
    VMR_g(:) = bg_VMR_in(:)

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

    rpar(1) = pl
    rpar(2) = p_vap
    rpar(3) = q_s
    rpar(4) = q_v_deep
    rpar(5) = tau_deep


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

    !! Calculate settling velocity !!

    !! Mass density of layer
    rho = (pl*mu*amu)/(kb * Tl) ! Mass density [g cm-3]

    !! Calculate dynamical viscosity for this layer - do square root mixing law from Rosner 2012
    top = 0.0_dp
    bot = 0.0_dp
    call eta_construct(sp_bg)
    do k = 1, n_gas
      eta_g = (5.0_dp/16.0_dp) * (sqrt(pi*(molg_g(k)*amu)*kb*Tl)/(pi*d_g(k)**2)) &
        & * ((((kb*Tl)/LJ_g(k))**(0.16_dp))/1.22_dp)
      top = top + sqrt(molg_g(k))*VMR_g(k)*eta_g
      bot = bot + sqrt(molg_g(k))*VMR_g(k)
    end do

    !! Mixture dynamical viscosity
    eta = top/bot

    !! Calculate mean free path for this layer
    mfp = (2.0_dp*eta/rho) * sqrt((pi * mu)/(8.0_dp*R_gas*Tl))

    !! Volume (or mass) weighted mean radius of particle assuming log-normal distribution
    rV = rm * exp(7.0_dp/2.0_dp * log(sigma)**2)

    !! Knudsen number
    Kn = mfp/rV

    !! Cunningham slip factor
    beta = 1.0_dp + Kn*(1.257_dp + 0.4_dp * exp(-1.1_dp/Kn))

    !! Settling velocity
    v_f = (2.0_dp * beta * grav * rV**2 * rho_c)/(9.0_dp * eta) & 
      & * (1.0_dp + ((0.45_dp*grav*rV**3*rho*rho_c)/(54.0_dp*eta**2))**(0.4_dp))**(-1.25_dp)

    !! Convert vf to mks
    v_f = v_f / 100.0_dp 

    deallocate(d_g, LJ_g, molg_g)

  end subroutine mini_cloud_1_s

  !! Calculate tracer fluxes
  subroutine dqdt(n,x,y,f,rpar,ipar)
    implicit none

    integer, intent(in) :: n
    integer, intent(in) :: ipar

    real(dp), intent(in) :: x 
    real(dp), dimension(n), intent(inout) :: y
    real(dp), dimension(*), intent(in) :: rpar

    real(dp), dimension(n), intent(out) :: f

    real(dp) :: sat

    !y(1) = q_v, y(2) = q_c
    !rpar(1) = pl(k), rpar(2) = p_vap, rpar(3) = q_s, rpar(4) = q_v_deep, rpar(5) = tau_deep

    y(1) = max(y(1),1e-30_dp)
    y(2) = max(y(2),1e-30_dp)

    !! Calculate dqdt for vapour and condensate
    sat = (y(1) * rpar(1))/rpar(2)

    !! Calculate dqdt given the supersaturation ratio
    if (sat < 0.99_dp) then

      ! Evaporate from q_c
      f(1) = min(rpar(3) - y(1), y(2))/tau_chem

    else if (sat > 1.01_dp) then

      ! Condense q_v toward the saturation ratio
      f(1) = -(y(1) - rpar(3))/tau_chem

    else
      f(1) = 0.0_dp
    end if

    f(2) = -f(1)

    ! Add replenishment to lower boundary at the tau_deep rate
    if (ipar == 1) then
      f(1) = f(1) - (y(1) - rpar(4))/rpar(5)
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
      p_vap_sp = exp(3.27860e1_dp - 8.65139e4_dp/(T + 4.80395e-1_dp))
    !case('TiC')
    !case('SiC')
    !case('CaTiO3')
    case('TiO2')
      ! Woitke & Helling (2004)
      p_vap_sp = exp(35.8027_dp - 74734.7_dp/T)
    case('Al2O3')
      ! Kozasa et al. (1989)
      p_vap_sp = exp(-73503.0_dp/T + 22.01_dp) * atm
    case('Fe')
      ! Elspeth note: Changed to Ackerman & Marley et al. (2001) expression
      if (T > 1800.0_dp) then
        p_vap_sp = exp(9.86_dp - 37120.0_dp/T) * bar
      else
        p_vap_sp = exp(15.71_dp - 47664.0_dp/T) * bar
      end if
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
      ! Morley et al. (2012)
      p_vap_sp = 10.0_dp**(7.490_dp - 20592_dp/T) * bar
     case('MnS')
      ! Morley et al. (2012)
      p_vap_sp = 10.0_dp**(11.532_dp - 23810.0_dp/T) * bar
    case('Na2S')
      ! Morley et al. (2012)
      p_vap_sp =  10.0_dp**(8.550_dp - 13889.0_dp/T) * bar
    case('ZnS')
      ! Morley et al. (2012)
      p_vap_sp = 10.0_dp**(12.812_dp - 15873.0_dp/T) * bar
    case('KCl')
      ! Morley et al. (2012)
      p_vap_sp = 10.0_dp**(7.611_dp - 11382_dp/T) * bar
    case('NaCl')
      ! Stull (1947)
      if (T < 83.0_dp) then ! Limiter for very cold T
        p_vap_sp = 10.0_dp**(5.07184_dp - 8388.497_dp / (83.0_dp - 82.638_dp)) * bar
      else
        p_vap_sp = 10.0_dp**(5.07184_dp - 8388.497_dp / (T - 82.638_dp)) * bar
      end if
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
    case default
      print*, 'Saturation: dust species not found: ', trim(sp), 'Stopping!'
      stop
    end select

  end function p_vap_sp

  !! eta for background gas
  subroutine eta_construct(sp_bg)
    implicit none

    character(len=20), dimension(:), intent(in) :: sp_bg
    
    integer :: n_bg, i

    n_bg = size(sp_bg)

    !! Construct the background gas arrays for eta (dynamical viscosity) calculation
    allocate(d_g(n_bg), LJ_g(n_bg), molg_g(n_bg))

    do i = 1, n_bg
      select case(sp_bg(i))

      case('OH')
        d_g(i) = d_OH
        LJ_g(i) = LJ_OH
        molg_g(i) = molg_OH
      case('H2')
        d_g(i) = d_H2
        LJ_g(i) = LJ_H2
        molg_g(i) = molg_H2
      case('H2O')
        d_g(i) = d_H2O
        LJ_g(i) = LJ_H2O
        molg_g(i) = molg_H2O
      case('H')
        d_g(i) = d_H
        LJ_g(i) = LJ_H
        molg_g(i) = molg_H
      case('CO')
        d_g(i) = d_CO
        LJ_g(i) = LJ_CO
        molg_g(i) = molg_CO
      case('CO2')
        d_g(i) = d_CO2
        LJ_g(i) = LJ_CO2
        molg_g(i) = molg_CO2
      case('O')
        d_g(i) = d_O
        LJ_g(i) = LJ_O
        molg_g(i) = molg_O
      case('CH4')
        d_g(i) = d_CH4
        LJ_g(i) = LJ_CH4
        molg_g(i) = molg_CH4
      case('C2H2')
        d_g(i) = d_C2H2
        LJ_g(i) = LJ_C2H2
        molg_g(i) = molg_C2H2
      case('NH3')
        d_g(i) = d_NH3
        LJ_g(i) = LJ_NH3
        molg_g(i) = molg_NH3
      case('N2')
        d_g(i) = d_N2
        LJ_g(i) = LJ_N2
        molg_g(i) = molg_N2 
      case('HCN')
        d_g(i) = d_HCN
        LJ_g(i) = LJ_HCN
        molg_g(i) = molg_HCN
      case('He')
        d_g(i) = d_He
        LJ_g(i) = LJ_He
        molg_g(i) = molg_He
      case default
        print*, 'Background gas species data not found: ', trim(sp_bg(i)), 'STOP'
        stop
      end select

    end do
    
  end subroutine eta_construct

end module mini_cloud_1_s_mod

