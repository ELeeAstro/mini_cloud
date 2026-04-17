module mini_cloud_vf_sat_adj_mod
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

  !! Construct required arrays for calculating gas mixtures
  real(dp), allocatable, dimension(:) :: d_g, LJ_g, molg_g, eta_g

  public :: mini_cloud_vf_sat_adj
  private :: eta_construct

  contains

  subroutine mini_cloud_vf_sat_adj(T_in, P_in, grav_in, mu_in, bg_VMR_in, rho_d, sp_bg, r_c, sigma, v_f, dist)
    implicit none

    ! Input variables
    character(len=20), dimension(:), intent(in) :: sp_bg
    integer, intent(in) :: dist
    real(dp), intent(in) :: T_in, P_in, mu_in, grav_in, rho_d, r_c, sigma
    real(dp), dimension(:), intent(in) :: bg_VMR_in

    real(dp), intent(out) :: v_f

    integer :: n_gas
    real(dp) :: T, mu, nd_atm, rho, p, grav, mfp, eta, cT
    real(dp), allocatable, dimension(:) :: VMR_g
    real(dp) :: r_v, Kn, Kn_b, beta, vf_s, vf_e, fx

    real(dp) :: A, B

    !! Find the number density of the atmosphere
    T = T_in             ! Convert temperature to K
    p = P_in * 10.0_dp   ! Convert pascal to dyne cm-2

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

    !! Thermal velocity
    cT = sqrt((2.0_dp * kb * T) / (mu * amu))

    !! Calculate dynamical viscosity for this layer
    call eta_construct(n_gas, sp_bg, VMR_g, T, eta)

    !! Calculate mean free path for this layer
    mfp = (2.0_dp*eta/rho) * sqrt((pi * mu)/(8.0_dp*R_gas*T))

    !! Volume (or mass) weighted mean radius of particle assuming given distribution
    if (dist == 1) then
      ! Lognormal distribution - r_c is median particle size and sigma the geometric std. dev.
      r_v = max(r_c * exp(7.0_dp/2.0_dp * log(sigma)**2), r_seed)
      ! print*, sigma, r_c*1e4, r_v*1e4
      ! stop
    else if (dist == 2) then
      ! Gamma distribution - r_c is now the number weighted radius and sigma is related to the trigamma function
      A = inv_trigamma_pos(log(sigma)**2)
      B = A/r_c
      r_v = max((A + 3.0_dp)/B, r_seed)
      ! print*,A, B, r_c*1e4, r_v*1e4
      ! stop
    end if

    !! Knudsen number
    Kn = mfp/r_v

    !! Cunningham slip factor (Jung et al. 2012)
    beta = 1.0_dp + Kn*(1.165_dp + 0.480_dp * exp(-0.101_dp/Kn))

    !! Settling velocity (Stokes regime)
    vf_s = (2.0_dp * beta * grav * r_v**2 * (rho_d - rho))/(9.0_dp * eta) & 
     & * (1.0_dp &
     & + ((0.45_dp*grav*r_c**3*rho*rho_d)/(54.0_dp*eta**2))**(0.4_dp))**(-1.25_dp)

    v_f = max(vf_s, 1.0e-30_dp)

    deallocate(d_g, LJ_g, molg_g, eta_g)

  end subroutine mini_cloud_vf_sat_adj

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

  pure real(dp) function inv_trigamma_pos(y_target) result(x)
    real(dp), intent(in) :: y_target
    real(dp) :: x_lo, x_hi, x_mid
    integer :: i

    if (y_target <= 0.0_dp) then
      x = huge(1.0_dp)
      return
    end if

    ! trigamma(x) is strictly decreasing for x > 0.
    x_lo = epsilon(1.0_dp)
    x_hi = max(1.0_dp, 1.0_dp/y_target + 1.0_dp/sqrt(y_target))

    do while (trigamma_pos(x_hi) > y_target)
      x_hi = 2.0_dp*x_hi
    end do

    do i = 1, 100
      x_mid = 0.5_dp*(x_lo + x_hi)
      if (trigamma_pos(x_mid) > y_target) then
        x_lo = x_mid
      else
        x_hi = x_mid
      end if
    end do

    x = 0.5_dp*(x_lo + x_hi)

  end function inv_trigamma_pos

  pure real(dp) function trigamma_pos(x) result(y)
    real(dp), intent(in) :: x
    real(dp) :: z, inv, inv2

    z = x
    y = 0.0_dp

    do while (z < 8.0_dp)
      y = y + 1.0_dp / (z*z)
      z = z + 1.0_dp
    end do

    inv  = 1.0_dp / z
    inv2 = inv * inv

    y = y + inv + 0.5_dp*inv2 + inv2*inv/6.0_dp &
          - inv2*inv2*inv/30.0_dp &
          + inv2*inv2*inv2*inv/42.0_dp &
          - inv2*inv2*inv2*inv2*inv/30.0_dp
  end function trigamma_pos


end module mini_cloud_vf_sat_adj_mod
