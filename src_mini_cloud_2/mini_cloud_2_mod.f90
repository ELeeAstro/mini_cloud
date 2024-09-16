module mini_cloud_2_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  ! Fortran 2008 intrinsic precisions - reccomended if possible
  !integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  !integer, parameter :: qp = REAL128

  !! Constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: R_gas = 8.31446261815324e7_dp
  real(dp), parameter :: amu = 1.66053906892e-24_dp
  real(dp), parameter :: Avo = 6.02214076e23_dp
  real(dp), parameter :: eV = 1.60217663e-12_dp
  real(dp), parameter :: third = 1.0_dp/3.0_dp
  real(dp), parameter :: twothird = 2.0_dp/3.0_dp

  !! Conversions to dyne
  real(dp), parameter :: bar = 1.0e6_dp ! bar to dyne
  real(dp), parameter :: atm = 1.01325e6_dp ! atm to dyne
  real(dp), parameter :: pa = 10.0_dp ! pa to dyne

  !! Global variables
  real(dp) :: T, mu, nd_atm, rho, p, grav, Rd
  real(dp) :: p_vap, rho_s, vth, sig, D

  !! Cloud global constants - some passed into from main
  real(dp) :: rho_d, mol_w_sp
  real(dp) :: r0, V0, m0, d0 ! Unit mass and volume
  real(dp), parameter :: r_seed = 1e-7_dp
  real(dp) :: V_seed, m_seed

  real(dp) :: mfp, eta

  !! Diameter, LJ potential and molecular weight for background gases
  real(dp), parameter :: d_OH = 3.06e-8_dp, LJ_OH = 100.0_dp * kb, molg_OH = 17.00734_dp  ! estimate
  real(dp), parameter :: d_H2 = 2.827e-8_dp, LJ_H2 = 59.7_dp * kb, molg_H2 = 2.01588_dp
  real(dp), parameter :: d_H2O = 2.641e-8_dp, LJ_H2O = 809.1_dp * kb, molg_H2O = 18.01528_dp
  real(dp), parameter :: d_H = 2.5e-8_dp, LJ_H = 30.0_dp * kb, molg_H = 1.00794_dp
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

  public :: mini_cloud_2, RHS_mom, jac_dum
  private :: calc_coal, calc_coag, calc_cond, calc_hom_nuc, calc_het_nuc, calc_seed_evap, &
    & p_vap_sp, surface_tension, eta_construct

  contains

  subroutine mini_cloud_2(T_in, P_in, grav_in, mu_in, bg_VMR_in, t_end, sp, sp_bg, q_v, q_0, q_1)
    implicit none

    ! Input variables
    character(len=20), intent(in) :: sp
    character(len=20), dimension(:), intent(in) :: sp_bg
    real(dp), intent(in) :: T_in, P_in, mu_in, grav_in, t_end
    real(dp), dimension(:), intent(in) :: bg_VMR_in

    ! Input/Output tracer values
    real(dp), intent(inout) :: q_v, q_0, q_1

    integer :: ncall, n

    ! Time controls
    real(dp) :: t_now

    ! DLSODE variables
    integer :: n_eq
    real(dp), allocatable, dimension(:) :: y
    real(dp), allocatable, dimension(:) :: rwork
    integer, allocatable, dimension(:) :: iwork
    integer :: itol, itask, istate, iopt, mf
    integer :: rworkdim, iworkdim
    real(dp) :: rtol, atol

    !! Work variables
    integer :: n_gas
    real(dp), allocatable, dimension(:) :: VMR_g
    real(dp) :: eta_g, bot, top

    !! Alter input values to mini-cloud units
    !! (note, some are obvious not not changed in case specific models need different conversion factors)
    n_eq = 3

    !! Find the number density of the atmosphere
    T = T_in             ! Convert temperature to K
    p = P_in * 10.0_dp   ! Convert pascal to dyne cm-2

    !! Allocate bg gas VMR
    n_gas = size(bg_VMR_in)
    allocate(VMR_g(n_gas))
    VMR_g(:) = bg_VMR_in(:)

    !! Change mu_in to mu
    mu = mu_in ! Convert mean molecular weight to mu [g mol-1]

    !! Change gravity to cgs [cm s-2]
    grav = grav_in * 100.0_dp

    !! Number density [cm-3] of layer
    nd_atm = p/(kb*T)  

    !! Mass density of layer
    rho = (p*mu*amu)/(kb * T) ! Mass density [g cm-3]

    !! Specific gas constant of layer [erg g-1 K-1]
    Rd = R_gas/mu

    !! Calculate dynamical viscosity for this layer - do square root mixing law from Rosner 2012
    call eta_construct(sp_bg)
    top = 0.0_dp
    bot = 0.0_dp
    do n = 1, n_gas
      eta_g = (5.0_dp/16.0_dp) * (sqrt(pi*(molg_g(n)*amu)*kb*T)/(pi*d_g(n)**2)) &
        & * ((((kb*T)/LJ_g(n))**(0.16_dp))/1.22_dp)
      top = top + sqrt(molg_g(n))*VMR_g(n)*eta_g
      bot = bot + sqrt(molg_g(n))*VMR_g(n)
    end do

    !! Mixture dynamical viscosity
    eta = top/bot

    !! Calculate mean free path for this layer
    mfp = (2.0_dp*eta/rho) * sqrt((pi * mu)/(8.0_dp*R_gas*T))

    !! Basic cloud properties
    m0 = mol_w_sp * amu
    V0 = m0 / rho_d
    r0 = ((3.0_dp*V0)/(4.0_dp*pi))**(third)
    V_seed = 4.0_dp/3.0_dp * pi * r_seed**3
    m_seed = V_seed * rho_d
    d0 = 2.0_dp * r0

    !! Saturation vapour pressure
    p_vap = p_vap_sp(sp, T)

    !! Saturation vapour density
    rho_s = p_vap/(Rd*T)

    !! Thermal velocity
    vth = sqrt((kb*T)/(2.0_dp*pi*m0))

    !! Gaseous diffusion constant
    D = 5.0_dp/(16.0_dp*Avo*d0**2*rho) * sqrt((R_gas*T*mu)/(2.0_dp*pi) * (mol_w_sp + mu)/mol_w_sp)

    !! Surface tension of material
    sig = surface_tension(sp, T)

    ! -----------------------------------------
    ! ***  parameters for the DLSODE solver  ***
    ! -----------------------------------------

    itask = 1
    istate = 1
    iopt = 1

    ! Problem is stiff (usual)
    ! mf = 21 - full jacobian matrix with jacobian save
    ! mf = 22 - internal calculated jacobian
    mf = 22
    rworkdim = 22 + 9*n_eq + n_eq**2
    iworkdim = 20 + n_eq
    allocate(rwork(rworkdim), iwork(iworkdim))

    itol = 1
    rtol = 1.0e-3_dp           ! Relative tolerances for each scalar
    atol = 1.0e-30_dp               ! Absolute tolerance for each scalar (floor value)

    rwork(:) = 0.0_dp
    iwork(:) = 0

    rwork(1) = 0.0_dp               ! Critical T value (don't integrate past time here)
    rwork(5) = 0.0_dp              ! Initial starting timestep (start low, will adapt in DVODE)
    rwork(6) = 0.0_dp       ! Maximum timestep

    iwork(5) = 0               ! Max order required
    iwork(6) = 100000               ! Max number of internal steps
    iwork(7) = 1                ! Number of error messages


    allocate(y(n_eq))

    !! Give tracer values to y
    y(1) = q_0
    y(2) = q_1
    y(3) = q_v 

    !! Limit y values
    y(:) = max(y(:),1e-30_dp)

    t_now = 0.0_dp

    ! Set the printing flag
    ! 0 = no printing, 1 = printing
    call xsetf(1)

    ncall = 0

    do while ((t_now < t_end) .and. (ncall < 100))

      call DLSODE (RHS_mom, n_eq, y, t_now, t_end, itol, rtol, atol, itask, &
        & istate, iopt, rwork, rworkdim, iwork, iworkdim, jac_dum, mf)

      ncall = ncall + 1

      if (mod(ncall,10) == 0) then
        istate = 1
      else  if (istate == -1) then
        istate = 2
      else if (istate < -1) then
        print*, 'dlsode: ', istate
        exit
      end if

    end do


    !! Limit y values
    y(:) = max(y(:),1e-30_dp)

    !! Give y values to tracers
    q_0 = y(1)
    q_1 = y(2)
    q_v = y(3)

    deallocate(y, rwork, iwork, d_g, LJ_g, molg_g)

  end subroutine mini_cloud_2

  subroutine RHS_mom(n_eq, time, y, f)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), intent(inout) :: time
    real(dp), dimension(n_eq), intent(inout) :: y
    real(dp), dimension(n_eq), intent(inout) :: f

    real(dp) :: f_coal, f_coag, f_nuc_hom, f_cond, f_nuc_het, f_seed_evap
    real(dp) :: m_c, r_c, Kn, beta, sat
    real(dp) :: p_v, n_v

    !! In this routine, you calculate the instanenous new fluxes (f) for each moment
    !! The current values of each moment (y) are typically kept constant
    !! Basically, you solve for the RHS of the ODE for each moment

    !! Limit y values
    y(:) = max(y(:),1e-30_dp)

    !! Convert y to real physical numbers to calculate f
    y(1) = y(1)*nd_atm ! Convert to real number density
    y(2) = y(2)*rho   ! Convert to real mass density
    y(3) = y(3)*rho   ! Convert to real mass density

    !! Find the true vapour VMR
    p_v = y(3) * Rd * T     !! Pressure of vapour
    n_v = p_v/(kb*T)        !! Number density of vapour

    !! Mean mass of particle
    m_c = max(y(2)/y(1), m_seed)

    !! Mass weighted mean radius of particle
    r_c = max(((3.0_dp*m_c)/(4.0_dp*pi*rho_d))**(third), r_seed)

    !! Knudsen number
    Kn = mfp/r_c

    !! Cunningham slip factor
    beta = 1.0_dp + Kn*(1.257_dp + 0.4_dp * exp(-1.1_dp/Kn))

    !! Find supersaturation ratio
    sat = p_v/p_vap

    !! Calculate condensation rate
    call calc_cond(n_eq, y, r_c, Kn, n_v, sat, f_cond)

    !! Calculate homogenous nucleation rate
    call calc_hom_nuc(n_eq, y, sat, n_v, f_nuc_hom)

    !! Calculate hetrogenous nucleation rate
    !call calc_het_nuc(n_eq, y, sat, r_c, f_nuc_het)
    f_nuc_het = 0.0_dp

    !! Calculate seed particle evaporation rate
    call calc_seed_evap(n_eq, y, m_c, f_cond, f_seed_evap)

    !! Calculate the coagulation rate
    call calc_coag(n_eq, y, m_c, r_c, beta, f_coag)

    !! Calculate the coalesence rate
    call calc_coal(n_eq, y, r_c, Kn, beta, f_coal)

    !! Calculate final net flux rate for each moment and vapour
    f(1) = (f_nuc_hom + f_nuc_het + f_seed_evap) + f_coag + f_coal
    f(2) = m_seed*(f_nuc_hom + f_nuc_het + f_seed_evap) + f_cond*y(1)
    f(3) = -f(2)

    !! Convert f to ratios
    f(1) = f(1)/nd_atm
    f(2) = f(2)/rho
    f(3) = f(3)/rho
      
    !! Convert y back to ratios
    y(1) = y(1)/nd_atm
    y(2) = y(2)/rho 
    y(3) = y(3)/rho

  end subroutine RHS_mom

  !! Condensation and evaporation
  subroutine calc_cond(n_eq, y, r_c, Kn, n_v, sat, dmdt)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), dimension(n_eq), intent(in) :: y 
    real(dp), intent(in) :: r_c, Kn, n_v, sat

    real(dp), intent(out) :: dmdt

    if (Kn >= 1.0_dp) then
      !! Kinetic regime [g s-1]
      dmdt = 4.0_dp * pi * r_c**2 * vth * m0 * n_v * (1.0_dp - 1.0_dp/sat)
    else
      !! Diffusive limited regime [g s-1]
      dmdt = 4.0_dp * pi * r_c * D * m0 * n_v * (1.0_dp - 1.0_dp/sat)
    end if

  end subroutine calc_cond

  !! Classical nucleation theory (CNT)
  subroutine calc_hom_nuc(n_eq, y, sat, n_v, J_hom)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), dimension(n_eq), intent(in) :: y 
    real(dp), intent(in) :: sat, n_v

    real(dp), intent(out) :: J_hom

    real(dp) :: ln_ss, theta_inf, N_inf, N_star, N_star_1, dg_rt, Zel, tau_gr
    real(dp) :: f0, kbT

    real(dp), parameter :: alpha = 1.0_dp
    real(dp), parameter :: Nf = 10.0_dp

    if (sat > 1.0_dp) then

      ! Efficency Variables
      ln_ss = log(sat) ! Natural log of saturation ratio
      f0 = 4.0_dp * pi * r0**2 ! Monomer Area
      kbT = kb * T         ! kb * T

      !! Find Nstar, theta_inf -> Dg/RT (Eq. 11 Lee et al. (2015a))
      theta_inf = (f0 * sig)/(kbT)  !Theta_infty eq.8 (Lee et al. (2015a)
      N_inf = (((twothird) * theta_inf) / ln_ss)**3

      !! Gail et al. (2014) ! note, typo in Lee et al. (2015a)
      N_star = 1.0_dp + (N_inf / 8.0_dp) &
        & * (1.0_dp + sqrt(1.0_dp + 2.0_dp*(Nf/N_inf)**third) &
        & - 2.0_dp*(Nf/N_inf)**third)**3
      N_star = max(1.00001_dp, N_star) ! Ensure no div 0
      N_star_1 = N_star - 1.0_dp

      !! delta G (N_star) / RT (Eq. 11 Lee et al. (2015a))
      dg_rt = theta_inf * (N_star_1 / (N_star_1**third + Nf**third))

      !! Calculate Zeldovich factor at N_star (Eq. 9 Lee et al. (2015a))
      Zel = sqrt((theta_inf / (9.0_dp * pi * (N_star_1)**(4.0_dp/3.0_dp))) &
        & * ((1.0_dp + 2.0_dp*(Nf/N_star_1)**third)/(1.0_dp + (Nf/N_star_1)**third)**3))

      !! Calculate !!inverse!! of tau_gr
      tau_gr = (f0 * N_star**(twothird)) * alpha * sqrt(kbT &
        & / (2.0_dp * pi * mol_w_sp * amu)) * n_v

      !! Finally calculate J_star [cm-3 s-1] ! Note underfloat limiter here
      J_hom = n_v * tau_gr * Zel * exp(max(-300.0_dp, N_star_1*ln_ss - dg_rt))
    else 
      !! Unsaturated, zero nucleation
      J_hom = 0.0_dp
    end if
    
  end subroutine calc_hom_nuc

  !! Seed particle evaporation
  subroutine calc_seed_evap(n_eq, y, m_c, f_cond, J_evap)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), dimension(n_eq), intent(in) :: y
    real(dp), intent(in) :: m_c, f_cond

    real(dp), intent(out) :: J_evap

    real(dp) :: tau_evap

    if ((f_cond >= 0.0_dp) .or. (y(1)/nd_atm <= 1e-29_dp)) then

      !! If growing or too little number density then evaporation can't take place
      J_evap = 0.0_dp

    else 

      !! Check if average mass is around 0.1% the seed particle mass
      !! This means the core is (probably) exposed to the air and can evaporate freely
      if (m_c <= (1.001_dp * m_seed)) then
        tau_evap = 1.0_dp !max(m_c/abs(f_cond),1.0_dp)
        !! Seed particle evaporation rate [cm-3 s-1]
        J_evap = -y(1)/tau_evap
      else
        !! There is still some mantle to evaporate from
        J_evap = 0.0_dp
      end if

    end if

  end subroutine calc_seed_evap

  !! Particle-particle Brownian coagulation
  subroutine calc_coag(n_eq, y, m_c, r_c, beta, f_coag)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), dimension(n_eq), intent(in) :: y 
    real(dp), intent(in) :: m_c, r_c, beta

    real(dp), intent(out) :: f_coag

    real(dp) :: phi, del_r, D_r, V_r, lam_r

    !! Particle diffusion rate
    D_r = (kb*T*beta)/(6.0_dp*pi*eta*r_c)

    !! Thermal velocity limit rate
    V_r = sqrt((8.0_dp*kb*T)/(pi*m_c))

    !! Ratio fraction
    lam_r = (8.0_dp*D_r)/(pi*V_r)

    !! Interpolation function
    del_r = ((2.0_dp*r_c + lam_r)**3 - (4.0_dp*r_c**2 + lam_r**2)**(1.5_dp))/(6.0_dp*r_c*lam_r) &
      & - 2.0_dp*r_c 

    !! Correction factor
    phi = 2.0_dp*r_c/(2.0_dp*r_c + sqrt(2.0_dp)*del_r) + (4.0_dp*D_r)/(r_c*sqrt(2.0_dp)*V_r) 

    !! Coagulation flux (Zeroth moment) [cm-3 s-1]
    f_coag = -(4.0_dp * kb * T * beta)/(3.0_dp * eta * phi)  * y(1)**2

  end subroutine calc_coag

  !! Particle-particle gravitational coalesence
  subroutine calc_coal(n_eq, y, r_c, Kn, beta, f_coal)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), dimension(n_eq), intent(in) :: y 
    real(dp), intent(in) :: r_c, Kn, beta

    real(dp), intent(out) :: f_coal

    real(dp) :: d_vf, vf, Stk, E
    real(dp), parameter :: eps = 0.5_dp

    !! Settling velocity
    vf = (2.0_dp * beta * grav * r_c**2 * rho_d)/(9.0_dp * eta) & 
      & * (1.0_dp &
      & + ((0.45_dp*grav*r_c**3*rho*rho_d)/(54.0_dp*eta**2))**(0.4_dp))**(-1.25_dp)

    !! Estimate differential velocity
    d_vf = eps * vf

    !! Calculate E
    if (Kn >= 1.0_dp) then
      !! E = 1 when Kn > 1
      E = 1.0_dp
    else
      !! Calculate Stokes number
      Stk = (vf * d_vf)/(grav * r_c)
      E = max(0.0_dp,1.0_dp - 0.42_dp*Stk**(-0.75_dp))
    end if

    !! Coalesence flux (Zeroth moment) [cm-3 s-1]
    f_coal = -2.0_dp*pi*r_c**2*y(1)**2*d_vf*E

  end subroutine calc_coal

  !! Classical hetrogenous nucleation theory
  subroutine calc_het_nuc(n_eq, y, sat, r_c, J_het)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), dimension(n_eq), intent(in) :: y
    real(dp), intent(in) :: sat, r_c

    real(dp), intent(out) :: J_het

    real(dp) :: r_crit, FG, Theta, Zel, Nl
    real(dp) :: f, mu_con, phi, f0, x
    real(dp) :: c_surf, nu, F_des

    if (sat > 1.0_dp) then

      !! Critical particle size
      r_crit = (2.0_dp * mol_w_sp * sig)/(rho_d*R_gas*T*log(sat))

      !! Formation energy of critical cluster - classical approximation
      FG = 4.0_dp/3.0_dp * pi * sig * r_crit**2

      !! Cluster-cluster diffusive interaction rate 
      Theta = p/(sqrt(2.0_dp*m0*kb*T))

      !! Number of vapour units in critical cluster size
      Nl = max((4.0_dp/3.0_dp * pi * r_crit**3)/V0, 1.0_dp)

      !! Zeldovich factor
      Zel = sqrt(FG/(3.0_dp*pi*kb*T*Nl**2))

      !! Cosine contact angle
      mu_con = 0.5_dp

      !! Begin shape factor calculation
      x = r_c/r_crit
      phi = sqrt(1.0_dp - 2.0_dp*mu_con*x + x**2)
      f0 = (x - mu_con)/phi

      !! Shape factor
      f = 0.5_dp * (1.0_dp + ((1.0_dp-mu_con*x)/phi)**3 + x**3*(2.0_dp - 3.0_dp*f0 + f0**3) &
        & + 3.0_dp*mu_con*x**2*(f0 - 1.0_dp))

      !! Desorption energy
      F_des = 0.5_dp * eV !~ 0.5 eV (Gao & Powell 2021)

      !! Oscillation frequency
      nu = 10e11_dp  * sqrt(F_des/(kb*m0)) ! Gao per. comm.

      !! Number density of condensate molecules on the nucleating surface
      c_surf = Theta/nu * exp(F_des/(kb*T))
      
      !! Heterogenous nucleation rate [cm-3 s-1]
      J_het = 4.0_dp * pi**2 * r_c**2 * r_crit**2 * Theta * c_surf * Zel * exp(-(FG*f)/(kb*T)) * y(1)

    else
      !! Unsaturated, zero nucleation
      J_het = 0.0_dp
    end if

  end subroutine calc_het_nuc

  !! Dummy jacobian subroutine required for dlsode
  subroutine jac_dum (NEQ, X, Y, ML, MU, PD, NROWPD)
    integer, intent(in) :: NEQ, ML, MU, NROWPD
    real(dp), intent(in) :: X
    real(dp), dimension(NEQ), intent(in) :: Y
    real(dp), dimension(NROWPD, NEQ), intent(inout) :: PD
  end subroutine jac_dum

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
      print*, 'Vapour pressure species not found: ', trim(sp), 'STOP'
      stop
    end select

  end function p_vap_sp

  !! Surface tension for each species
  real(dp) function surface_tension(sp, T)
    implicit none

    character(len=20), intent(in) :: sp
    real(dp), intent(in) :: T

    real(dp) :: TC

    TC = T - 273.15_dp

    ! Return surface tension in erg cm-2
    select case (sp)
    case('C')
      ! Tabak et al. (1995)
      sig = 1400.0_dp
    case('CaTiO3')
      sig = 915.0_dp
    case('TiC')
      sig = 1240.0_dp
    case('TiO2')
      ! Sindel et al. (2022)
      sig = 589.79_dp - 0.0708_dp * T
    case('Fe')
      ! http://www.kayelaby.npl.co.uk/general_physics/2_2/2_2_5.html
      sig = 1862.0_dp - 0.39_dp * (TC - 1530.0_dp)
      ! Pradhan et al. (2009)
      !sig = 2858.0_dp - 0.51_dp * T
    case('Fe2O3')
      sig = 410.0_dp
    case('FeO')
      ! Janz 1988
      sig = 585.0_dp
    case('Al2O3')
      ! Pradhan et al. (2009)
      sig = 1024.0_dp - 0.177_dp * T
      ! Kozasa et al. (1989)
      !sig = 690.0_dp
    case('MgSiO3')
      ! Janz 1988
      sig = 197.3_dp + 0.098_dp * T
    case('Mg2SiO4')
      ! Kozasa et al. (1989)
      sig = 436.0_dp
    case('SiC')
      ! Nozawa et al. (2003)
      sig = 1800.0_dp
    case('SiO')
      ! Gail and Sedlmayr (1986)
      sig = 500.0_dp
    case('SiO2')
      ! Pradhan et al. (2009)
      !sig = 243.2_dp - 0.013_dp * T
      ! Janz 1988
      sig = 243.2_dp + 0.031_dp * T
    case('Cr')
      ! http://www.kayelaby.npl.co.uk/general_physics/2_2/2_2_5.html
      sig = 1642.0_dp - 0.20_dp * (TC - 1860.0_dp)      
    case('MnS')
      sig = 2326.0_dp
    case('Na2S')
      sig = 1033.0_dp
    case('KCl')
      ! Janz 1988
      sig = 175.57_dp - 0.07321_dp * T
    case('NaCl')
      ! Janz 1988
      sig = 191.16_dp - 0.07188_dp * T
    case('ZnS')
      sig = 860.0_dp
    case('H2O')
      ! Hale and Plummer (1974)
      sig = 141.0_dp - 0.15_dp*TC
    case('NH3')
      ! Weast et al. (1988)
      sig = 23.4_dp
    case('NH4Cl')
      sig = 56.0_dp
    case('NH4SH')
      sig = 50.0_dp
    case('CH4')
      ! USCG (1984)
      sig = 14.0_dp
    case('H2S')
      ! Nehb and Vydra (2006)
      sig = 58.1_dp
    case('S2','S8')
      ! Fanelli 1950
      sig = 60.8_dp
    case default
      print*, 'Species surface tension not found: ', trim(sp), 'STOP'
      stop
    end select

    surface_tension = max(10.0_dp, sig)

      ! Pradhan et al. (2009):
      !Si : 732 - 0.086*(T - 1685.0)
      !MgO : 1170 - 0.636*T
      !CaO : 791 - 0.0935*T

  end function surface_tension

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

end module mini_cloud_2_mod