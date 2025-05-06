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
  real(dp), parameter :: V_seed = 4.0_dp/3.0_dp * pi * r_seed**3

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

  public :: mini_cloud_vf
  private :: eta_construct

  contains

  subroutine mini_cloud_vf(T_in, P_in, grav_in, mu_in, bg_VMR_in, rho_d, sp_bg, q_0, q_1, q_2, v_f)
    implicit none

    ! Input variables
    character(len=20), dimension(:), intent(in) :: sp_bg
    real(dp), intent(in) :: T_in, P_in, mu_in, grav_in, rho_d, q_0, q_1, q_2
    real(dp), dimension(:), intent(in) :: bg_VMR_in

    real(dp), dimension(3), intent(out) :: v_f

    integer :: n_bg
    real(dp) :: T, mu, nd_atm, rho, p, grav, mfp, eta, cT
    real(dp), allocatable, dimension(:) :: VMR_bg
    real(dp) :: N_c, rho_c, Z_c, sig2, lam, nu
    real(dp) :: m_c, r_c, m_c2, r_c2, r_n, m_seed
    real(dp) :: vf_s, vf_e, fx, Rey, St, Ep, gam_fac, lgnu, lgnu1, lgnu2

    real(dp) :: nu_n, Kn, Kn_m, Kn_m2, Kn_b, Kn_n

    real(dp), parameter :: A = 1.639_dp

    !! Find the number density of the atmosphere
    T = T_in             ! Convert temperature to K
    p = P_in * 10.0_dp   ! Convert pascal to dyne cm-2

    !! Number density [cm-3] of layer
    nd_atm = p/(kb*T)  

    !! Zero velocity if little amount of clouds
    ! if (q_0*nd_atm < 1e-10_dp) then
    !   v_f = 0.0_dp
    !   return
    ! end if

    n_bg = size(bg_VMR_in)
    allocate(VMR_bg(n_bg))
    VMR_bg(:) = bg_VMR_in(:)

    !! Change mu_in to mu
    mu = mu_in ! Convert mean molecular weight to mu [g mol-1]

    !! Change gravity to cgs [cm s-2]
    grav = grav_in * 100.0_dp

    !! Mass density of layer
    rho = (p*mu*amu)/(kb * T) ! Mass density [g cm-3]

    !! Thermal velocity
    cT = sqrt((2.0_dp * kb * T) / (mu * amu))

    !! Calculate dynamical viscosity for this layer - do square root mixing law from Rosner 2012
    call eta_construct(n_bg, sp_bg, VMR_bg, T, eta)

    !! Calculate mean free path for this layer
    mfp = (2.0_dp*eta/rho) * sqrt((pi * mu)/(8.0_dp*R_gas*T))

    m_seed = V_seed * rho_d

    N_c = q_0 * nd_atm
    rho_c = q_1 * rho
    Z_c = q_2 * rho**2

    !! Calculate vf from final results of interaction
    !! Mean mass of particle
    m_c = rho_c/N_c
    m_c2 = Z_c/rho_c

    !! Calculate lambda and nu gamma distribution parameters
    sig2 = max(Z_c/rho_c - (rho_c/N_c)**2,m_seed**2)
    nu = max(m_c**2/sig2,0.01_dp)
    nu = min(nu,100.0_dp)
    lam = m_c/nu

    !! Mass weighted mean radius of particle
    r_c = max(((3.0_dp*m_c)/(4.0_dp*pi*rho_d))**(third), r_seed)

    !! Mass^2 weighted mean radius of particle
    r_c2 = max(((3.0_dp*m_c2)/(4.0_dp*pi*rho_d))**(third), r_seed)

    lgnu  = log_gamma(nu)
    lgnu1 = log_gamma(nu + 1.0_dp)
    lgnu2 = log_gamma(nu + 2.0_dp)

    !! Number weighted mean radius
    r_n = max(r_c * nu**(-1.0_dp/3.0_dp) * exp(log_gamma(nu + 1.0_dp/3.0_dp) - lgnu), r_seed)

    !! Monodisperse Knudsen number
    Kn = mfp/r_c
    Kn_b = min(Kn, 100.0_dp)

    !! Population averaged Knudsen number for n, m and m^2
    nu_n = max(nu,0.3334_dp)
    Kn_n = Kn * nu_n**(1.0_dp/3.0_dp) * &
      & exp(log_gamma(nu_n - 1.0_dp/3.0_dp) - log_gamma(nu_n))
    Kn_m = Kn * nu**(1.0_dp/3.0_dp) * &
      & exp(log_gamma(nu + 2.0_dp/3.0_dp) - lgnu1)
    Kn_m2 = Kn * nu**(1.0_dp/3.0_dp) * & 
      & exp(log_gamma(nu + 5.0_dp/3.0_dp) - lgnu2)

    !! Now find moment dependent settling velocities
    St = (2.0_dp * grav * r_c**2 * (rho_d - rho))/(9.0_dp * eta) 
    Ep = (sqrt(pi)*grav*rho_d*r_c)/(2.0_dp*cT*rho)

    !! Zeroth moment
    !! Settling velocity (Stokes regime)
    Rey = (1.0_dp + ((0.45_dp*grav*r_n**3*rho*rho_d)/(54.0_dp*eta**2))**(0.4_dp))**(-1.25_dp)
    gam_fac = (nu**(-2.0/3.0) * exp(log_gamma(nu + 2.0_dp/3.0_dp) - lgnu) & 
      & + A*Kn_b*nu**(-1.0/3.0) * exp(log_gamma(nu + 1.0_dp/3.0_dp) - lgnu))
    vf_s = St * gam_fac * Rey

    !! Settling velocity (Epstein regime)
    gam_fac =  (nu**(-1.0/3.0) * exp(log_gamma(nu + 1.0_dp/3.0_dp) - lgnu))
    vf_e = Ep * gam_fac

    !! tanh interpolation function
    fx = 0.5_dp * (1.0_dp - tanh(2.0_dp*log10(Kn_n)))

    !! Interpolation for settling velocity
    v_f(1) = fx*vf_s + (1.0_dp - fx)*vf_e

    !! First moment
    !! Settling velocity (Stokes regime)
    Rey = (1.0_dp + ((0.45_dp*grav*r_c**3*rho*rho_d)/(54.0_dp*eta**2))**(0.4_dp))**(-1.25_dp)
    gam_fac = (nu**(-2.0/3.0) * exp(log_gamma(nu + 5.0_dp/3.0_dp) - lgnu1) & 
      & + A*Kn_b*nu**(-1.0/3.0) * exp(log_gamma(nu + 4.0_dp/3.0_dp) - lgnu1))
    vf_s = St * gam_fac * Rey

    !! Settling velocity (Epstein regime)
    gam_fac =  (nu**(-1.0/3.0) * exp(log_gamma(nu + 4.0_dp/3.0_dp) - lgnu1))
    vf_e = Ep * gam_fac

    !! tanh interpolation function
    fx = 0.5_dp * (1.0_dp - tanh(2.0_dp*log10(Kn_m)))

    !! Interpolation for settling velocity
    v_f(2) = fx*vf_s + (1.0_dp - fx)*vf_e

    !! Second moment
    !! Settling velocity (Stokes regime)
    Rey = (1.0_dp + ((0.45_dp*grav*r_c2**3*rho*rho_d)/(54.0_dp*eta**2))**(0.4_dp))**(-1.25_dp)
    gam_fac = (nu**(-2.0/3.0) * exp(log_gamma(nu + 8.0_dp/3.0_dp) - lgnu2) & 
      & + A*Kn_b*nu**(-1.0/3.0) * exp(log_gamma(nu + 7.0_dp/3.0_dp) - lgnu2))
    vf_s = St * gam_fac * Rey

    !! Settling velocity (Epstein regime)
    gam_fac =  (nu**(-1.0/3.0) * exp(log_gamma(nu + 7.0_dp/3.0_dp) - lgnu2))
    vf_e = Ep * gam_fac

    !! tanh interpolation function
    fx = 0.5_dp * (1.0_dp - tanh(2.0_dp*log10(Kn_m2)))

    !! Interpolation for settling velocity
    v_f(3) = fx*vf_s + (1.0_dp - fx)*vf_e

    deallocate(d_g, LJ_g, molg_g, eta_g)

  end subroutine mini_cloud_vf

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

end module mini_cloud_vf_mod