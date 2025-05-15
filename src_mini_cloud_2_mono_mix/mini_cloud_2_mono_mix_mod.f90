module mini_cloud_2_mono_mix_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  ! Fortran 2008 intrinsic precisions - recommended if possible
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
  real(dp), parameter :: mmHg = 1333.22387415_dp  ! mmHg to dyne

  !! Global atmospheric variables
  real(dp) :: T, mu, nd_atm, rho, p, grav, Rd
  real(dp) :: mfp, eta, nu, cT, met

  real(dp), parameter :: r_seed = 1e-7_dp
  real(dp), parameter :: V_seed = 4.0_dp/3.0_dp * pi * r_seed**3

  !! Cloud properties container
  type cld_sp

    character(len=20) :: sp
    real(dp) :: rho_d, mol_w_sp, Rd_v
    real(dp) :: p_vap, vth, sig, D
    real(dp) :: m_seed, Kn_crit, alp_c
    real(dp) :: r0, V0, m0, d0, r_seed, mol_w_v, m_v

    real(dp) :: v2c

    integer :: inuc
    real(dp) :: sat, Nf, alp_nuc

  end type cld_sp

  integer :: ndust

  !! Cloud properties array
  type(cld_sp), dimension(:), allocatable :: cld

  !$omp threadprivate(T, mu, nd_atm, rho, p, grav, Rd, mfp, eta, nu, cT, met, cld, ndust)


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

  !! Construct required arrays for calculating gas mixtures
  real(dp), allocatable, dimension(:) :: d_g, LJ_g, molg_g, eta_g

  !$omp threadprivate(d_g, LJ_g, molg_g, eta_g)

  public :: mini_cloud_2_mono_mix, RHS_mom, jac_dum
  private :: calc_coal, calc_coag, calc_cond, calc_hom_nuc, calc_seed_evap, &
    & p_vap_sp, surface_tension, eta_construct

  contains

  subroutine mini_cloud_2_mono_mix(ilay, T_in, P_in, grav_in, mu_in, met_in, bg_VMR_in, t_end, sp_in, sp_bg, &
    & n_in, q_v, q_0, q_1)
    implicit none

    ! Input variables
    integer, intent(in) :: n_in, ilay
    character(len=20), dimension(n_in), intent(in) :: sp_in
    character(len=20), dimension(:), intent(in) :: sp_bg
    real(dp), intent(in) :: T_in, P_in, mu_in, grav_in, t_end, met_in
    real(dp), dimension(:), intent(in) :: bg_VMR_in

    ! Input/Output tracer values
    real(dp), intent(inout) :: q_0
    real(dp), dimension(n_in), intent(inout) :: q_v, q_1

    integer :: ncall

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
    integer :: n_bg, j
    real(dp), allocatable, dimension(:) :: VMR_bg

    ndust = n_in

    !! Alter input values to mini-cloud units
    !! (note, some are obvious not not changed in case specific models need different conversion factors)
    n_eq = ndust*2+1

    !! Find the number density of the atmosphere
    T = T_in             ! Convert temperature to K
    p = P_in * 10.0_dp   ! Convert pascal to dyne cm-2

    !! Allocate bg gas VMR
    n_bg = size(bg_VMR_in)
    allocate(VMR_bg(n_bg))
    VMR_bg(:) = bg_VMR_in(:)

    !! Change mu_in to mu
    mu = mu_in ! Convert mean molecular weight to mu [g mol-1]

    !! Approx log10 metallicity of atmosphere
    met = met_in

    !! Change gravity to cgs [cm s-2]
    grav = grav_in * 100.0_dp

    !! Number density [cm-3] of layer
    nd_atm = p/(kb*T)  

    !! Mass density of layer
    rho = (p*mu*amu)/(kb * T) ! Mass density [g cm-3]

    !! Thermal velocity
    cT = sqrt((2.0_dp * kb * T) / (mu * amu))

    !! Specific gas constant of layer [erg g-1 K-1]
    Rd = R_gas/mu

    !! Calculate dynamical viscosity for this layer
    call eta_construct(n_bg, sp_bg, VMR_bg, T, eta)

    !! Mixture kinematic viscosity
    nu = eta/rho

    !! Calculate mean free path for this layer
    mfp = (2.0_dp*eta/rho) * sqrt((pi * mu)/(8.0_dp*R_gas*T))

    !! Get the basic cloud species properties
    call cloud_sp(ndust, sp_in, T, rho, mu, mfp)

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
    rwork(5) = 0.0_dp              ! Initial starting timestep (start low, will adapt in DLSODE)
    rwork(6) = 0.0_dp       ! Maximum timestep
    rwork(7) = 0.0_dp       ! Minimum timestep

    iwork(5) = 0               ! Max order required
    iwork(6) = 1000000               ! Max number of internal steps
    iwork(7) = 1                ! Number of error messages


    allocate(y(n_eq))

    !! Give tracer values to y
    y(1) = q_0
    y(2:2+ndust-1) = q_1(:)
    y(2+ndust:) = q_v(:)

    !! Limit y values
    y(:) = max(y(:),1e-30_dp)

    t_now = 0.0_dp

    ! Set the printing flag
    ! 0 = no printing, 1 = printing
    call xsetf(0)

    ncall = 0

    do while ((t_now < t_end) .and. (ncall < 100))

      call DLSODE (RHS_mom, n_eq, y, t_now, t_end, itol, rtol, atol, itask, &
        & istate, iopt, rwork, rworkdim, iwork, iworkdim, jac_dum, mf)

      ncall = ncall + 1

      !if (mod(ncall,10) == 0) then
        !istate = 1
      !else  if (istate == -1) then
        !istate = 2
      if (istate < 0) then
        print*, 'dlsode: ', istate, ilay, t_now
        exit
      end if

    end do

    !! Limit y values
    y(:) = max(y(:),1e-30_dp)

    !! Give y values to tracers
    q_0 = y(1)
    q_1(:) = y(2:2+ndust-1)
    q_v(:) = y(2+ndust:)

    deallocate(y, rwork, iwork, d_g, LJ_g, molg_g, eta_g, cld)

  end subroutine mini_cloud_2_mono_mix

  subroutine RHS_mom(n_eq, time, y, f)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), intent(inout) :: time
    real(dp), dimension(n_eq), intent(inout) :: y
    real(dp), dimension(n_eq), intent(inout) :: f

    real(dp) :: f_coal, f_coag
    real(dp) :: m_c, r_c, Kn, beta, vf_s, vf_e, vf, Kn_b
    real(dp) :: fx

    integer :: j
    real(dp) :: N_c, rho_c_t, rho_d_m, V_tot
    real(dp), dimension(ndust) :: rho_c, rho_v, V_mix
    real(dp), dimension(ndust) :: p_v, n_v
    real(dp), dimension(ndust) :: f_nuc_hom, f_cond, f_seed_evap

    !! In this routine, you calculate the instantaneous new fluxes (f) for each moment
    !! The current values of each moment (y) are typically kept constant
    !! Basically, you solve for the RHS of the ODE for each moment

    !! Limit y values
    y(:) = max(y(:),1e-30_dp)

    !! Convert y to real physical numbers to calculate f
    N_c = y(1)*nd_atm ! Convert to real number density
    rho_c(:) = y(2:2+ndust-1)*rho   ! Convert to real mass density
    rho_v(:) = y(2+ndust:)*rho   ! Convert to real mass density

    !! Find the true vapour VMR
    p_v(:) = rho_v(:) * cld(:)%Rd_v * T     !! Pressure of vapour
    n_v(:) = p_v(:)/(kb*T)        !! Number density of vapour

    !! Mean mass of particle
    rho_c_t = sum(rho_c(:))
    m_c = max(rho_c_t/N_c, cld(1)%m_seed)

    rho_d_m = 0.0_dp
    do j = 1, ndust
      rho_d_m = rho_d_m + (rho_c(j)/rho_c_t) * cld(j)%rho_d
    end do

    !! Bulk material volume mixing ratio
    V_tot = sum(rho_c(:)/cld(:)%rho_d) ! Total condensed volume
    V_mix(:) = (rho_c(:)/cld(:)%rho_d)/V_tot ! Condensed volume mixing ratio
    !V_mix(:) = rho_c(:)/rho_c_t

    !! Mass weighted mean radius of particle
    r_c = max(((3.0_dp*m_c)/(4.0_dp*pi*rho_d_m))**(third), r_seed)

    !! Knudsen number
    Kn = mfp/r_c

    !! Cunningham slip factor (Kim et al. 2005)
    Kn_b = min(Kn, 100.0_dp)
    beta = 1.0_dp + Kn_b*(1.165_dp + 0.483_dp * exp(-0.997_dp/Kn_b))

    !! Settling velocity (Stokes regime)
    vf_s = (2.0_dp * beta * grav * r_c**2 * (rho_d_m - rho))/(9.0_dp * eta) & 
     & * (1.0_dp &
     & + ((0.45_dp*grav*r_c**3*rho*rho_d_m)/(54.0_dp*eta**2))**(0.4_dp))**(-1.25_dp)

    !! Settling velocity (Epstein regime)
    vf_e = (sqrt(pi)*grav*rho_d_m*r_c)/(2.0_dp*cT*rho)

    !! tanh interpolation function
    fx = 0.5_dp * (1.0_dp - tanh(2.0_dp*log10(Kn)))

    !! Interpolation for settling velocity
    vf = fx*vf_s + (1.0_dp - fx)*vf_e

    !! Find supersaturation ratio
    cld(:)%sat = p_v(:)/cld(:)%p_vap

    !! Calculate condensation rate
    call calc_cond(ndust, r_c, Kn, n_v(:), V_mix(:), f_cond)

    !! Calculate homogenous nucleation rate
    call calc_hom_nuc(ndust, n_v(:), f_nuc_hom)

    !! Calculate seed particle evaporation rate
    call calc_seed_evap(ndust, N_c, m_c, f_cond, f_seed_evap)

    !! Calculate the coagulation rate
    call calc_coag(m_c, r_c, beta, f_coag)

    !! Calculate the coalescence rate
    call calc_coal(r_c, Kn, vf, f_coal)

    !! Calculate final net flux rate for each moment and vapour
    f(1) = (f_nuc_hom(1) + f_seed_evap(1)) + (f_coag + f_coal)*N_c**2
    f(2:2+ndust-1) = cld(:)%m_seed*(f_nuc_hom(:)  + f_seed_evap(:)) + f_cond(:)*N_c
    f(2+ndust:) = -cld(:)%m_seed*(f_nuc_hom(:)  + f_seed_evap(:)) - cld(:)%v2c*f_cond(:)*N_c

    !! Convert f to ratios
    f(1) = f(1)/nd_atm
    f(2:) = f(2:)/rho

    ! Check if condensation from vapour is viable
    do j = 1, ndust
      if (rho_v(j)/rho <= 1e-28_dp) then
        if (f(2+j-1) > 0.0_dp) then
          f(2+j-1) = 0.0_dp
          f(2+ndust+j-1) = 0.0_dp
         !print*, j,'in'
        end if
      end if
    end do

    ! Check if evaporation from condensate is viable
    do j = 1, ndust
      if (rho_c(j)/rho <= 1e-28_dp) then
        if (f(2+ndust+j-1) > 0.0_dp) then
          f(2+j-1) = 0.0_dp
          f(2+ndust+j-1) = 0.0_dp
         !print*, j,'in'
        end if
      end if
    end do

  end subroutine RHS_mom

  !! Condensation and evaporation
  subroutine calc_cond(ndust, r_c, Kn, n_v, V_mix, dmdt)
    implicit none

    integer, intent(in) :: ndust
    real(dp), intent(in) :: r_c, Kn
    real(dp), dimension(ndust), intent(in) :: n_v, V_mix

    real(dp), dimension(ndust), intent(out) :: dmdt

    integer :: j
    real(dp) :: dmdt_low, dmdt_high, Knd, Kn_crit, fx

    do j = 1, ndust

      !! Diffusive limited regime (Kn << 1) [g s-1]
      dmdt_low = 4.0_dp * pi * r_c * cld(j)%D * cld(j)%m0 * n_v(j) &
        & *  (1.0_dp - 1.0_dp/cld(j)%sat)

      !! Free molecular flow regime (Kn >> 1) [g s-1]
      dmdt_high = 4.0_dp * pi * r_c**2 * cld(j)%vth * cld(j)%m0 * n_v(j) * cld(j)%alp_c & 
        & * (1.0_dp - 1.0_dp/cld(j)%sat)

      !! If evaporation, weight rate by current condensed volume ratio (Woitke et al. 2020)
      if (cld(j)%sat < 1.0_dp) then
        dmdt_high = dmdt_high * V_mix(j) 
        dmdt_low = dmdt_low * V_mix(j) 
      end if

      !! Critical Knudsen number
      Kn_crit = Kn * (dmdt_high/dmdt_low)

      !! Kn' (Woitke & Helling 2003)
      Knd = Kn/Kn_crit

      !! tanh interpolation function
      fx = 0.5_dp * (1.0_dp - tanh(2.0_dp*log10(Knd)))

      !! Mass change rate
      dmdt(j) = dmdt_low * fx + dmdt_high * (1.0_dp - fx)
    end do

  end subroutine calc_cond

  !! modified classical nucleation theory (MCNT)
  subroutine calc_hom_nuc(ndust, n_v, J_hom)
    implicit none

    integer, intent(in) :: ndust
    real(dp), dimension(ndust), intent(in) :: n_v

    real(dp), dimension(ndust), intent(out) :: J_hom

    integer :: j
    real(dp) :: ln_ss, theta_inf, N_inf, N_star, N_star_1, dg_rt, Zel, tau_gr
    real(dp) :: f0, kbT

    real(dp) :: alpha, Nf

    do j = 1, ndust

      if (cld(j)%inuc == 0) then
        J_hom(j) = 0.0_dp
        cycle
      end if

      if (cld(j)%sat > 1.0_dp) then

        if (cld(j)%sp == 'SiO') then
          !! Special nucleation rate for SiO (Gail et al. 2016)
          J_hom(j) = n_v(j)**2 * exp(1.33_dp - 4.40e12_dp / (T**3 * log(cld(j)%sat)**2))

        else

          ! Efficency Variables
          ln_ss = log(cld(j)%sat) ! Natural log of saturation ratio
          f0 = 4.0_dp * pi * cld(j)%r0**2 ! Monomer Area
          kbT = kb * T         ! kb * T

          alpha = cld(j)%alp_nuc
          Nf = cld(j)%Nf

          !! Find Nstar, theta_inf -> Dg/RT (Eq. 11 Lee et al. (2015a))
          theta_inf = (f0 * cld(j)%sig)/(kbT)  !Theta_infty eq.8 (Lee et al. (2015a)
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
            & / (2.0_dp * pi * cld(j)%mol_w_v * amu)) * n_v(j)

          !! Finally calculate J_star [cm-3 s-1] ! Note underfloat limiter here
          J_hom(j) = n_v(j) * tau_gr * Zel * exp(max(-300.0_dp, N_star_1*ln_ss - dg_rt))
        end if

      else 
        !! Unsaturated, zero nucleation
        J_hom(j) = 0.0_dp
      end if
    end do
    
  end subroutine calc_hom_nuc

  !! Seed particle evaporation
  subroutine calc_seed_evap(ndust, N_c, m_c, f_cond, J_evap)
    implicit none

    integer, intent(in) :: ndust
    real(dp), dimension(ndust), intent(in) :: f_cond
    real(dp), intent(in) :: N_c, m_c

    real(dp), dimension(ndust), intent(out) :: J_evap

    integer :: j
    real(dp) :: tau_evap

    do j = 1, ndust

      if (cld(j)%inuc == 0) then
        J_evap(j) = 0.0_dp
        cycle
      end if

      if ((cld(j)%sat >= 1.0_dp)) then

        !! If growing or too little number density then evaporation can't take place
        J_evap(j) = 0.0_dp

      else 

        !! Check if average mass is around 0.1% the seed particle mass
        !! This means the core is (probably) exposed to the air and can evaporate freely
        if (m_c <= (1.001_dp * cld(j)%m_seed)) then
          tau_evap = 0.1_dp !m_c/abs(f_cond)
          !! Seed particle evaporation rate [cm-3 s-1]
          J_evap(j) = -N_c/tau_evap
        else
          !! There is still some mantle to evaporate from
          J_evap(j) = 0.0_dp
        end if

      end if
    end do

  end subroutine calc_seed_evap

  !! Particle-particle Brownian coagulation
  subroutine calc_coag(m_c, r_c, beta, f_coag)
    implicit none

    real(dp), intent(in) :: m_c, r_c, beta

    real(dp), intent(out) :: f_coag

    real(dp) :: phi, D_r, V_r
    real(dp) :: Knd

    !! Particle diffusion rate
    D_r = (kb*T*beta)/(6.0_dp*pi*eta*r_c)

    !! Thermal velocity limit rate
    V_r = sqrt((8.0_dp*kb*T)/(pi*m_c))

    !! Moran (2022) method using diffusive Knudsen number
    Knd = (8.0_dp*D_r)/(pi*V_r*r_c)
    phi = 1.0_dp/sqrt(1.0_dp + pi**2/8.0_dp * Knd**2)
    f_coag = (-4.0_dp*kb*T*beta)/(3.0_dp*eta) * phi

  end subroutine calc_coag

  !! Particle-particle gravitational coalesence
  subroutine calc_coal(r_c, Kn, vf, f_coal)
    implicit none

    real(dp), intent(in) :: r_c, Kn, vf

    real(dp), intent(out) :: f_coal

    real(dp) :: d_vf,  Stk, E
    real(dp), parameter :: eps = 0.5_dp

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

    !! Coalesence flux (Zeroth moment) [cm3 s-1]
    f_coal = -2.0_dp*pi*r_c**2*d_vf*E

  end subroutine calc_coal

  !! Vapour pressure for each species
  real(dp) function p_vap_sp(sp, T)
    implicit none

    character(len=20), intent(in) :: sp
    real(dp), intent(in) :: T

    real(dp) :: TC, A, B, C, f

    ! Return vapour pressure in dyne
    select case(sp)
    case('C')
      ! Gail & Sedlmayr (2013) - I think...
      p_vap_sp = exp(3.27860e1_dp - 8.65139e4_dp/(T + 4.80395e-1_dp))
    case('TiC')
      ! Kimura et al. (2023)
      p_vap_sp = 10.0_dp**(-33600.0_dp/T + 7.652_dp) * atm
    case('SiC')
      ! Elspeth 5 polynomial JANAF-NIST fit
      p_vap_sp =  exp(-9.51431385e4_dp/T + 3.72019157e1_dp + 1.09809718e-3_dp*T &
        & -5.63629542e-7_dp*T**2 + 6.97886017e-11_dp*T**3)
    case('CaTiO3')
      ! Wakeford et al. (2017) -  taken from VIRGA
      p_vap_sp = 10.0_dp**(-72160.0_dp/T + 30.24_dp - log10(p/1e6_dp) - 2.0_dp*met) * bar
      ! Kozasa et al. (1987)
      !p_vap_sp = exp(-79568.2_dp/T + 42.0204_dp) * atm
    case('TiO2')
      ! GGChem 5 polynomial NIST fit
      p_vap_sp = exp(-7.70443e4_dp/T +  4.03144e1_dp - 2.59140e-3_dp*T &
        &  + 6.02422e-7_dp*T**2 - 6.86899e-11_dp*T**3)
    case('VO')
      ! NIST 5 param fit
      p_vap_sp = exp(-6.74603e4_dp/T + 3.82717e1_dp - 2.78551e-3_dp*T &
        & + 5.72078e-7_dp*T**2 - 7.41840e-11_dp*T**3)
    case('Al2O3')
      ! Wakeford et al. (2017) - taken from CARMA
      p_vap_sp = 10.0_dp**(17.7_dp - 45892.6_dp/T - 1.66_dp*met) * bar
      ! Kozasa et al. (1987)
      !p_vap_sp = exp(-73503.0_dp/T + 22.005_dp) * atm
    case('Fe')
      ! Visscher et al. (2010) - taken from CARMA
      p_vap_sp = 10.0_dp**(7.23_dp - 20995.0_dp/T) * bar
      ! Elspeth note: Changed to Ackerman & Marley et al. (2001) expression
      ! if (T > 1800.0_dp) then
      !   p_vap_sp = exp(9.86_dp - 37120.0_dp/T) * bar
      ! else
      !   p_vap_sp = exp(15.71_dp - 47664.0_dp/T) * bar
      ! end if
    case('FeS')
      ! GGChem 5 polynomial NIST fit
      p_vap_sp = exp(-5.69922e4_dp/T + 3.86753e1_dp - 4.68301e-3_dp*T &
        & + 1.03559e-6_dp*T**2 - 8.42872e-11_dp*T**3)
    case('FeO')
      ! GGChem 5 polynomial NIST fit
      p_vap_sp = exp(-6.30018e4_dp/T + 3.66364e1_dp - 2.42990e-3_dp*T &
        & + 3.18636e-7_dp*T**2)
    case('Mg2SiO4')
      ! Visscher et al. (2010)/Visscher notes - taken from CARMA
      p_vap_sp = 10.0_dp**(14.88_dp - 32488.0_dp/T - 1.4_dp*met - 0.2_dp*log10(p/1e6_dp)) * bar
      ! Kozasa et al. (1989) - Seems to be too high
      !p_vap_sp = p*10.0**(8.25_dp -  27250.0_dp/T - log10(p/1e6_dp) + 3.58_dp)
    case('MgSiO3','MgSiO3_amorph')
      ! Vissher - taken from VIRGA
      p_vap_sp = 10.0_dp**(13.43_dp - 28665.0_dp/T - met) * bar
      ! Ackerman & Marley (2001)
      !p_vap_sp = exp(-58663.0_dp/T + 25.37_dp) * bar
    case('MgO')
      ! GGChem 5 polynomial NIST fit
      p_vap_sp = exp(-7.91838e4_dp/T + 3.57312e1_dp + 1.45021e-4_dp*T &
        &  - 8.47194e-8*T**2 + 4.49221e-12_dp*T**3)
    case('SiO2','SiO2_amorph')
      ! GGChem 5 polynomial NIST fit
      p_vap_sp = exp(-7.28086e4_dp/T + 3.65312e1_dp - 2.56109e-4_dp*T &
      & - 5.24980e-7_dp*T**2 + 1.53343E-10_dp*T**3) 
    case('SiO')
      ! Gail et al. (2013)
      p_vap_sp = exp(-49520.0_dp/T + 32.52_dp)
    case('Cr')
      ! GGChem 5 polynomial NIST fit
      p_vap_sp = exp(-4.78455e+4_dp/T + 3.22423e1_dp - 5.28710e-4_dp*T & 
        &  - 6.17347e-8_dp*T**2 + 2.88469e-12_dp*T**3)
     case('MnS')
      ! Morley et al. (2012)
      p_vap_sp = 10.0_dp**(11.532_dp - 23810.0_dp/T - met) * bar
    case('Na2S')
      ! Morley et al. (2012)
      p_vap_sp =  10.0_dp**(8.550_dp - 13889.0_dp/T - 0.5_dp*met) * bar
    case('ZnS')
      ! Elspeth 5 polynomial Barin data fit
      p_vap_sp = exp(-4.75507888e4_dp/T + 3.66993865e1_dp - 2.49490016e-3_dp*T &
        &  + 7.29116854e-7_dp*T**2 - 1.12734453e-10_dp*T**3)
      ! Morley et al. (2012)
      !p_vap_sp = 10.0_dp**(12.812_dp - 15873.0_dp/T - met) * bar        
    case('KCl')
      ! GGChem 5 polynomial NIST fit
      p_vap_sp = exp(-2.69250e4_dp/T + 3.39574e+1_dp - 2.04903e-3_dp*T &
        & -2.83957e-7_dp*T**2 + 1.82974e-10_dp*T**3)
      ! Morley et al. (2012)
      !p_vap_sp =  10.0_dp**(7.611_dp - 11382.0_dp/T) * bar
    case('NaCl')
      ! GGChem 5 polynomial NIST fit
      p_vap_sp = exp(-2.79146e4_dp/T + 3.46023e1_dp - 3.11287e3_dp*T & 
        & + 5.30965e-7_dp*T**2 -2.59584e-12_dp*T**3)
    case('NH4Cl')
      p_vap_sp = 10.0_dp**(7.0220_dp - 4302.0_dp/T) * bar
    case('H2O')
      TC = T - 273.15_dp
      ! Huang (2018) - A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice
      if (TC < 0.0_dp) then
        ! f = 0.99882_dp * exp(0.00000008_dp * p/pa)
        p_vap_sp = exp(43.494_dp - (6545.8_dp/(TC + 278.0_dp)))/(TC + 868.0_dp)**2.0_dp * pa
      else
        ! f = 1.00071_dp * exp(0.000000045_dp * p/pa
        p_vap_sp = exp(34.494_dp - (4924.99_dp/(TC + 237.1_dp)))/(TC + 105.0_dp)**1.57_dp * pa
      end if
      ! Ackerman & Marley (2001) H2O liquid & ice vapour pressure expressions
      !if (T > 1048.0_dp) then
      !  p_vap_sp = 6.0e8_dp
      !else if (T < 273.16_dp) then
      !  p_vap_sp = 6111.5_dp * exp((23.036_dp * TC - TC**2/333.7_dp)/(TC + 279.82_dp))
      !else
      !  p_vap_sp = 6112.1_dp * exp((18.729_dp * TC - TC**2/227.3_dp)/(TC + 257.87_dp)) 
      !end if
    case('NH3')
      ! Ackerman & Marley (2001) NH3 ice vapour pressure expression from Weast (1971)
      p_vap_sp = exp(10.53_dp - 2161.0_dp/T - 86596.0_dp/T**2)  * bar
    case('CH4')
      ! Lodders & Fegley (1998) - directly taken from VIRGA
      if (T < 90.68_dp) then
        C = -16.043_dp/8.3143_dp * (2.213_dp - 2.650_dp)
        B = -16.043_dp/8.3143_dp * (611.10_dp + (2.213_dp - 2.650_dp) * 90.68_dp )
        A = 0.11719_dp * 90.68_dp**(-C) * exp(-B/90.68_dp)
      else
        C = -16.043_dp/8.3143_dp * (2.213_dp - 3.370_dp)
        B = -16.043_dp/8.3143_dp * (552.36_dp + (2.213_dp - 3.370_dp) * 90.68_dp)
        A = 0.11719_dp * 90.68_dp**(-C) * exp(-B/90.68_dp)
      end if
      p_vap_sp = A * T**C * exp(B/T) * bar

      !--- Prydz, R.; Goodwin, R.D., J. Chem. Thermodyn., 1972, 4,1 ---
      !if (T < 0.5_dp) then ! Limiter for very very very cold T
      !  p_vap_sp = 10.0_dp**(3.9895_dp - 443.028_dp/(0.5_dp-0.49_dp)) * bar       
      !else
      !  p_vap_sp = 10.0_dp**(3.9895_dp - 443.028_dp/(T-0.49_dp)) * bar
      !end if
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
      print*, 'Vapour pressure species not found: ', trim(sp), 'STOP'
      stop
    end select

  end function p_vap_sp

  !! Surface tension for each species
  real(dp) function surface_tension(sp, T)
    implicit none

    character(len=20), intent(in) :: sp
    real(dp), intent(in) :: T

    real(dp) :: TC, sig

    TC = T - 273.15_dp

    ! Return surface tension in erg cm-2
    select case (sp)
    case('C')
      ! Tabak et al. (1995)
      sig = 1400.0_dp
    case('TiC')
      ! Chigai et al. (1999)
      sig = 1242.0_dp
    case ('SiC')
      ! Nozawa et al. (2003)
      sig = 1800.0_dp
    case('CaTiO3')
      ! Kozasa et al. (1987)
      sig = 494.0_dp      
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
      !sig = 175.57_dp - 0.07321_dp * T
      sig = 160.4_dp - 0.07_dp*T
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
  subroutine eta_construct(n_bg, sp_bg, VMR_bg, T, eta_out)
    implicit none

    integer, intent(in) :: n_bg
    character(len=20), dimension(:), intent(in) :: sp_bg
    real(dp), dimension(n_bg), intent(in) :: VMR_bg
    real(dp), intent(in) :: T

    real(dp), intent(out) :: eta_out
    
    integer :: i, j
    real(dp) :: bot, Eij, part
    real(dp), dimension(n_bg) :: y

    !! Construct the background gas arrays for eta (dynamical viscosity) calculation
    allocate(d_g(n_bg), LJ_g(n_bg), molg_g(n_bg), eta_g(n_bg))

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

    !! Davidson (1993) mixing rule
    
    !! First calculate each species eta
    do i = 1, n_bg
      eta_g(i) = (5.0_dp/16.0_dp) * (sqrt(pi*(molg_g(i)*amu)*kb*T)/(pi*d_g(i)**2)) &
        & * ((((kb*T)/LJ_g(i))**(0.16_dp))/1.22_dp)
    end do

    !! Calculate y values
    bot = 0.0_dp
    do i = 1, n_bg
      bot = bot + VMR_bg(i) * sqrt(molg_g(i))
    end do
    y(:) = (VMR_bg(:) * sqrt(molg_g(:)))/bot

    !! Calculate fluidity following Davidson equation
    eta_out = 0.0_dp
    do i = 1, n_bg
      do j = 1, n_bg
        Eij = ((2.0_dp*sqrt(molg_g(i)*molg_g(j)))/(molg_g(i) + molg_g(j)))**0.375
        part = (y(i)*y(j))/(sqrt(eta_g(i)*eta_g(j))) * Eij
        eta_out = eta_out + part
      end do
    end do

    !! Viscosity is inverse fluidity
    eta_out = 1.0_dp/eta_out

  end subroutine eta_construct

  subroutine cloud_sp(ndust, sp_in, T, rho, mu, mfp)
    implicit none

    integer, intent(in) :: ndust
    character(len=20), dimension(ndust), intent(in) :: sp_in
    real(dp), intent(in) :: T, rho, mu, mfp

    integer :: j

    !! Allocate global cld type array
    allocate(cld(ndust))

    !! Look over number of dust species
    do j = 1, ndust

      !! First get the constant values from a select case

      cld(j)%sp = sp_in(j)

      !! Calculate the basic cloud properties
      select case(trim(cld(j)%sp))

      case('C')

        cld(j)%mol_w_sp = 12.0107_dp
        cld(j)%rho_d = 2.27_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 1
        cld(j)%Nf = 5.0_dp
        cld(j)%alp_nuc = 1.0_dp

        ! Key reaction: C -> C[s]
        cld(j)%mol_w_v = cld(j)%mol_w_sp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('TiC')

        cld(j)%mol_w_sp = 59.8777_dp
        cld(j)%rho_d = 4.93_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: 2 Ti + C2H2 -> 2 TiC[s] + H2
        cld(j)%mol_w_v = 47.8670_dp
        cld(j)%v2c = (2.0_dp*cld(j)%mol_w_v)/(2.0_dp*cld(j)%mol_w_sp)

      case('SiC')

        cld(j)%mol_w_sp = 40.0962_dp
        cld(j)%rho_d = 3.21_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: Si2C -> SiC[s] + Si
        cld(j)%mol_w_v = 68.1817_dp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('CaTiO3')

        cld(j)%mol_w_sp = 135.943_dp
        cld(j)%rho_d = 3.98_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: Ca + TiO + 2 H2O -> CaTiO3[s] + 2 H2
        cld(j)%mol_w_v = 40.0780_dp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('Al2O3')

        cld(j)%mol_w_sp = 101.961_dp
        cld(j)%rho_d = 3.986_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: 2 Al + 3 H2O -> Al2O3[s] + 3 H2
        cld(j)%mol_w_v = 26.98153860_dp
        cld(j)%v2c = 2.0_dp*cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('TiO2')

        cld(j)%mol_w_sp = 79.866_dp
        cld(j)%rho_d = 4.23_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 1
        cld(j)%Nf = 0.0_dp
        cld(j)%alp_nuc = 1.0_dp

        ! Key reaction: TiO2 -> TiO2[s]
        cld(j)%mol_w_v = cld(j)%mol_w_sp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('VO')

        cld(j)%mol_w_sp = 66.94090_dp
        cld(j)%rho_d = 5.76_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0
         
        ! Key reaction: VO -> VO[s]
        cld(j)%mol_w_v = cld(j)%mol_w_sp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('Fe')

        cld(j)%mol_w_sp = 55.845_dp
        cld(j)%rho_d = 7.87_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: Fe -> Fe[s]
        cld(j)%mol_w_v = cld(j)%mol_w_sp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('FeO')

        cld(j)%mol_w_sp = 71.8444_dp
        cld(j)%rho_d = 5.99_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: Fe + H2O -> FeO[s] + H2
        cld(j)%mol_w_v = 55.845_dp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('Fe2O3')

        cld(j)%mol_w_sp = 159.6882_dp
        cld(j)%rho_d = 5.25_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: 2 Fe + 3 H2O -> Fe2O3[s] + 3 H2
        cld(j)%mol_w_v = 55.845_dp
        cld(j)%v2c = 2.0_dp*cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('FeS')

        cld(j)%mol_w_sp = 87.91_dp
        cld(j)%rho_d = 4.84_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: Fe + H2S -> FeS[s] + H2
        cld(j)%mol_w_v = 55.845_dp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('Mg2SiO4')

        cld(j)%mol_w_sp = 140.693_dp
        cld(j)%rho_d = 3.21_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: 2 Mg + SiO + 3 H2O -> Mg2SiO4[s] + 3 H2
        cld(j)%mol_w_v = 24.305_dp
        cld(j)%v2c = 2.0_dp*cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('MgSiO3')

        cld(j)%mol_w_sp = 100.389_dp
        cld(j)%rho_d = 3.19_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: Mg + SiO + 2 H2O -> MgSiO3[s] + H2
        cld(j)%mol_w_v = 24.305_dp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('MgO')

        cld(j)%mol_w_sp = 40.3044_dp
        cld(j)%rho_d = 3.6_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: Mg + H2O -> MgO[s] + H2
        cld(j)%mol_w_v = 24.305_dp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp   

      case('MgS')

        cld(j)%mol_w_sp = 56.37_dp
        cld(j)%rho_d = 2.84_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: Mg + H2S -> MgS[s] + H2
        cld(j)%mol_w_v = 24.305_dp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp   

      case('SiO2')

        cld(j)%mol_w_sp = 60.084_dp
        cld(j)%rho_d = 2.648_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: SiO + H2O -> SiO2[s] + H2O
        cld(j)%mol_w_v = 44.085_dp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('SiO')

        cld(j)%mol_w_sp = 44.085_dp
        cld(j)%rho_d = 2.18_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 1
        cld(j)%Nf = 0.0_dp
        cld(j)%alp_nuc = 1.0_dp

        ! Key reaction: SiO -> SiO[s]
        cld(j)%mol_w_v = cld(j)%mol_w_sp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('SiS')

        cld(j)%mol_w_sp = 60.1505_dp
        cld(j)%rho_d = 2.18_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: SiS -> SiS[s]
        cld(j)%mol_w_v = cld(j)%mol_w_sp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp   

      case('Cr')

        cld(j)%mol_w_sp = 51.996_dp
        cld(j)%rho_d = 7.19_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: Cr -> Cr[s]
        cld(j)%mol_w_v = cld(j)%mol_w_sp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('MnS')

        cld(j)%mol_w_sp = 87.003_dp
        cld(j)%rho_d =  4.08_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: Mn + H2S -> MnS[s] + H2
        cld(j)%mol_w_v = 54.938045_dp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('Na2S')

        cld(j)%mol_w_sp = 78.0445_dp
        cld(j)%rho_d =  1.856_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: 2 Na + H2S -> Na2S[s] + H2
        cld(j)%mol_w_v = 22.98976928_dp
        cld(j)%v2c = 2.0_dp*cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('ZnS')

        cld(j)%mol_w_sp = 97.445_dp
        cld(j)%rho_d = 4.09_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: Zn + H2S -> ZnS[s] + H2
        cld(j)%mol_w_v =  65.38_dp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('KCl')

        cld(j)%mol_w_sp = 74.551_dp
        cld(j)%rho_d = 1.99_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 1
        cld(j)%Nf = 5.0_dp
        cld(j)%alp_nuc = 1.0_dp

        ! Key reaction: KCl -> KCl[s]
        cld(j)%mol_w_v = cld(j)%mol_w_sp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('NaCl')

        cld(j)%mol_w_sp =  58.443_dp
        cld(j)%rho_d = 2.165_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 1
        cld(j)%Nf = 5.0_dp
        cld(j)%alp_nuc = 1.0_dp

        ! Key reaction: NaCl -> NaCl[s]
        cld(j)%mol_w_v = cld(j)%mol_w_sp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('NH4Cl')

        cld(j)%mol_w_sp = 53.4915_dp
        cld(j)%rho_d = 1.52_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: NH3 + HCl -> NH4C[s]
        cld(j)%mol_w_v = 36.4609_dp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('H2O')

        cld(j)%mol_w_sp = 18.015_dp
        cld(j)%rho_d = 0.93_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: H2O -> H2O[l/s]
        cld(j)%mol_w_v = cld(j)%mol_w_sp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('NH3')

        cld(j)%mol_w_sp = 17.031_dp
        cld(j)%rho_d = 0.87_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: NH3 -> NH3[s]
        cld(j)%mol_w_v = cld(j)%mol_w_sp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('CH4')

        cld(j)%mol_w_sp = 16.043_dp
        cld(j)%rho_d = 0.656_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: CH4 -> CH4[s]
        cld(j)%mol_w_v = cld(j)%mol_w_sp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('CO')

        cld(j)%mol_w_sp = 28.0101_dp
        cld(j)%rho_d = 1.14_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: CO -> CO[s]
        cld(j)%mol_w_v = cld(j)%mol_w_sp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('CO2')

        cld(j)%mol_w_sp = 44.0095_dp
        cld(j)%rho_d = 1.98_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: CO2 -> CO2[s]
        cld(j)%mol_w_v = cld(j)%mol_w_sp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp


      case('H2SO4')

        cld(j)%mol_w_sp = 98.0785_dp
        cld(j)%rho_d = 1.8302_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: SO3 + H2O -> H2SO4[l]
        cld(j)%mol_w_v = 80.0632_dp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case ('NH4SH')

        cld(j)%mol_w_sp = 51.1114_dp
        cld(j)%rho_d = 1.17_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: NH3 + H2S -> NH4SH[s]
        cld(j)%mol_w_v = 34.08_dp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case('H2S')

        cld(j)%mol_w_sp = 34.08_dp
        cld(j)%rho_d = 1.12_dp
        cld(j)%alp_c = 1.0_dp
        cld(j)%inuc = 0

        ! Key reaction: H2S -> H2S[s]
        cld(j)%mol_w_v = cld(j)%mol_w_sp
        cld(j)%v2c = cld(j)%mol_w_v/cld(j)%mol_w_sp

      case default
        print*, 'Cloud species not found: ',   trim(cld(j)%sp)
        print*, 'STOP'
        stop     
      end select

      !! Find the constant values for the cloud species
      cld(j)%m0 = cld(j)%mol_w_sp * amu
      cld(j)%V0 = cld(j)%m0 / cld(j)%rho_d
      cld(j)%r0 = ((3.0_dp*cld(j)%V0)/(4.0_dp*pi))**(third)
      if (cld(j)%inuc == 0) then
        cld(j)%m_seed = 0.0_dp
      else
        cld(j)%m_seed = V_seed * cld(j)%rho_d
      end if
      cld(j)%d0 = 2.0_dp * cld(j)%r0

      !! Find the constant values for the cloud species
      !! that are specific for the layer

      !! Saturation vapour pressure
      cld(j)%p_vap = p_vap_sp(cld(j)%sp, T)

      !! Thermal velocity of vapour
      cld(j)%vth = sqrt((kb*T)/(2.0_dp*pi*cld(j)%mol_w_v*amu))

      !! Gaseous diffusion constant of vapour - d0 estimated here!!
      cld(j)%D = 5.0_dp/(16.0_dp*Avo*cld(j)%d0**2*rho) * &
        & sqrt((R_gas*T*mu)/(2.0_dp*pi) * (cld(j)%mol_w_v + mu)/cld(j)%mol_w_v)

      !! Surface tension of material
      cld(j)%sig = surface_tension(cld(j)%sp, T)

      !! Specific gas constant of vapour [erg g-1 K-1]
      cld(j)%Rd_v = R_gas/cld(j)%mol_w_v

    end do


  end subroutine cloud_sp

  !! Dummy jacobian subroutine required for dlsode
  subroutine jac_dum (NEQ, X, Y, ML, MU, PD, NROWPD)
    integer, intent(in) :: NEQ, ML, MU, NROWPD
    real(dp), intent(in) :: X
    real(dp), dimension(NEQ), intent(in) :: Y
    real(dp), dimension(NROWPD, NEQ), intent(inout) :: PD
  end subroutine jac_dum

end module mini_cloud_2_mono_mix_mod