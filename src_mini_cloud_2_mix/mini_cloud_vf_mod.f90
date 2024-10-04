module mini_cloud_vf_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: R_gas = 8.31446261815324e7_dp
  real(dp), parameter :: amu = 1.66053906892e-24_dp
  
  real(dp), parameter :: third = 1.0_dp/3.0_dp
  real(dp), parameter :: twothird = 2.0_dp/3.0_dp

  real(dp), parameter :: r_seed = 1e-7_dp

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

  public :: mini_cloud_vf
  private :: eta_construct

  contains

  subroutine mini_cloud_vf(T_in, P_in, grav_in, mu_in, bg_VMR_in, rho_s, sp_bg, q_0, q_1s, v_f)
    implicit none

    ! Input variables
    character(len=20), dimension(:), intent(in) :: sp_bg
    real(dp), intent(in) :: T_in, P_in, mu_in, grav_in, q_0
    real(dp), dimension(:), intent(in) :: bg_VMR_in, q_1s, rho_s

    real(dp), intent(out) :: v_f

    integer :: n_gas, n, n_dust
    real(dp) :: T, mu, nd_atm, rho, p, grav, mfp, eta, rho_d, rho_t
    real(dp), allocatable, dimension(:) :: VMR_g, V_frac, m_c_s, V_c_s
    real(dp) :: top, bot, eta_g, m_c, r_c, Kn, beta

    n_dust = size(q_1s)

    !! Find the number density of the atmosphere
    T = T_in             ! Convert temperature to K
    p = P_in * 10.0_dp   ! Convert pascal to dyne cm-2

    !! Number density [cm-3] of layer
    nd_atm = p/(kb*T)  

    !! Zero velocity if little amount of clouds
    if (q_0*nd_atm < 1e-10_dp) then
      v_f = 0.0_dp
      return
    end if

    n_gas = size(bg_VMR_in)
    allocate(VMR_g(n_gas))
    VMR_g(:) = bg_VMR_in(:)

    !! Change mu_in to mu
    mu = mu_in ! Convert mean molecular weight to mu [g mol-1]

    !! Change gravity to cgs [cm s-2]
    grav = grav_in * 100.0_dp

    !! Mass density of layer
    rho = (p*mu*amu)/(kb * T) ! Mass density [g cm-3]

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

    !! Calculate vf from final results of interation

    allocate(V_frac(n_dust),m_c_s(n_dust),V_c_s(n_dust))

    !! Mean mass of particle per species
    m_c_s(:) = (q_1s(:)*rho)/(q_0*nd_atm)
    !! Mean volumes of particle per species
    V_c_s(:) = m_c_s(:)/rho_s(:)
    !! Volume fraction of species in grain
    V_frac(:) = max(V_c_s(:)/sum(V_c_s(:)),1e-99_dp)
    !! Volume weighted bulk density
    rho_d = sum(V_frac(:)*rho_s(:))

    !! Mean mass of particle
    rho_t = sum(q_1s(:))
    m_c = (rho_t*rho)/(q_0*nd_atm)

    !! Mass weighted mean radius of particle
    r_c = max(((3.0_dp*m_c)/(4.0_dp*pi*rho_d))**(third), r_seed)

    !! Knudsen number
    Kn = mfp/r_c

    !! Cunningham slip factor
    beta = 1.0_dp + Kn*(1.257_dp + 0.4_dp * exp(-1.1_dp/Kn))

    !! Settling velocity
    v_f = (2.0_dp * beta * grav * r_c**2 * rho_d)/(9.0_dp * eta) & 
      & * (1.0_dp + & 
      & ((0.45_dp*grav*r_c**3*rho*rho_d)/(54.0_dp*eta**2))**(0.4_dp))**(-1.25_dp)

    deallocate(d_g, LJ_g, molg_g, V_frac, m_c_s, V_c_s)

  end subroutine mini_cloud_vf

    !! eta for background gas
  subroutine eta_construct(sp_bg)
    implicit none

    character(len=20), dimension(:), intent(in) :: sp_bg
    
    integer :: n_bg, i

    n_bg = size(sp_bg)

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

end module mini_cloud_vf_mod