module mini_cloud_vf_num_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: R_gas = 8.31446261815324e7_dp
  real(dp), parameter :: amu = 1.66053906892e-24_dp
  real(dp), parameter :: third = 1.0_dp/3.0_dp

  real(dp), parameter :: r_seed = 1.0e-7_dp
  real(dp), parameter :: V_seed = 4.0_dp/3.0_dp * pi * r_seed**3
  real(dp), parameter :: v_floor = 1.0e-30_dp
  real(dp), parameter :: tiny_pos = 1.0e-300_dp
  real(dp), parameter :: r_dyn_max = 1.0e-1_dp

  integer, parameter :: nq_gl = 16

  ! Gauss-Laguerre nodes and weights:
  ! integral_0^inf exp(-x) f(x) dx ~= sum_i lag_w(i) f(lag_x(i)).
  real(dp), parameter :: lag_x(nq_gl) = [ &
    8.76494104789278111e-02_dp, 4.62696328915081445e-01_dp, &
    1.14105777483122739e+00_dp, 2.12928364509838053e+00_dp, &
    3.43708663389320668e+00_dp, 5.07801861454976766e+00_dp, &
    7.07033853504823462e+00_dp, 9.43831433639193840e+00_dp, &
    1.22142233688661594e+01_dp, 1.54415273687816175e+01_dp, &
    1.91801568567531362e+01_dp, 2.35159056939919076e+01_dp, &
    2.85787297428821425e+01_dp, 3.45833987022866225e+01_dp, &
    4.19404526476883319e+01_dp, 5.17011603395433212e+01_dp ]

  real(dp), parameter :: lag_w(nq_gl) = [ &
    2.06151714957799909e-01_dp, 3.31057854950890795e-01_dp, &
    2.65795777644212639e-01_dp, 1.36296934296374711e-01_dp, &
    4.73289286941242782e-02_dp, 1.12999000803392299e-02_dp, &
    1.84907094352627899e-03_dp, 2.04271915308275255e-04_dp, &
    1.48445868739810618e-05_dp, 6.82831933087107655e-07_dp, &
    1.88102484107965474e-08_dp, 2.86235024297378986e-10_dp, &
    2.12707903322408582e-12_dp, 6.29796700251783083e-15_dp, &
    5.05047370003545801e-18_dp, 4.16146237037272733e-22_dp ]

  real(dp), parameter :: d_OH = 3.06e-8_dp, LJ_OH = 100.0_dp * kb, molg_OH = 17.00734_dp
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

  real(dp), allocatable, dimension(:) :: d_g, LJ_g, molg_g, eta_g

  !$omp threadprivate(d_g, LJ_g, molg_g, eta_g)

  public :: mini_cloud_vf
  private :: eta_construct, gamma_moment_velocities, terminal_velocity_particle

contains

  subroutine mini_cloud_vf(T_in, P_in, grav_in, mu_in, bg_VMR_in, rho_d, sp_bg, &
    & ndust, q_0, q_1, q_2, v_f)
    implicit none

    integer, intent(in) :: ndust
    character(len=20), dimension(:), intent(in) :: sp_bg
    real(dp), intent(in) :: T_in, P_in, mu_in, grav_in, q_0, q_2
    real(dp), dimension(:), intent(in) :: bg_VMR_in
    real(dp), dimension(ndust), intent(in) :: rho_d, q_1
    real(dp), dimension(:), intent(out) :: v_f

    integer :: n_bg
    real(dp) :: T, mu, nd_atm, rho, p, grav, mfp, eta
    real(dp), allocatable, dimension(:) :: VMR_bg
    real(dp) :: N_c, sig2, nu, lam, Z_c_t
    real(dp) :: m_c, m_seed, rho_c_t, rho_d_m, V_tot
    real(dp), dimension(ndust) :: rho_c

    v_f(:) = v_floor

    T = T_in
    p = P_in * 10.0_dp
    mu = mu_in
    grav = grav_in * 100.0_dp
    if (T <= 0.0_dp .or. p <= 0.0_dp .or. mu <= 0.0_dp .or. grav <= 0.0_dp) return

    nd_atm = p/(kb*T)
    rho = (p*mu*amu)/(kb*T)

    n_bg = size(bg_VMR_in)
    allocate(VMR_bg(n_bg))
    VMR_bg(:) = bg_VMR_in(:)

    call eta_construct(n_bg, sp_bg, VMR_bg, T, eta)
    if (eta <= 0.0_dp .or. rho <= 0.0_dp) then
      call cleanup_local(VMR_bg)
      return
    end if

    mfp = (2.0_dp*eta/rho) * sqrt((pi*mu)/(8.0_dp*R_gas*T))

    m_seed = V_seed * rho_d(1)
    N_c = max(q_0*nd_atm, tiny_pos)
    rho_c(:) = q_1(:)*rho
    rho_c_t = max(sum(rho_c(:)), N_c*m_seed)
    Z_c_t = max(q_2*rho**2, rho_c_t**2/N_c)

    m_c = max(rho_c_t/N_c, m_seed)

    V_tot = sum(rho_c(:)/rho_d(:))
    if (V_tot > 0.0_dp) then
      rho_d_m = rho_c_t/V_tot
    else
      rho_d_m = rho_d(1)
    end if
    rho_d_m = max(rho_d_m, tiny_pos)

    sig2 = max(Z_c_t/N_c - (rho_c_t/N_c)**2, m_seed**2)
    nu = max(m_c**2/sig2, 0.01_dp)
    nu = min(nu, 100.0_dp)
    lam = max(m_seed, m_c/nu)

    call gamma_moment_velocities(nu, lam, rho_d_m, rho, eta, mfp, grav, v_f)

    call cleanup_local(VMR_bg)

  contains

    subroutine cleanup_local(VMR_bg_local)
      real(dp), allocatable, dimension(:), intent(inout) :: VMR_bg_local

      if (allocated(d_g)) deallocate(d_g)
      if (allocated(LJ_g)) deallocate(LJ_g)
      if (allocated(molg_g)) deallocate(molg_g)
      if (allocated(eta_g)) deallocate(eta_g)
      if (allocated(VMR_bg_local)) deallocate(VMR_bg_local)
    end subroutine cleanup_local

  end subroutine mini_cloud_vf

  subroutine gamma_moment_velocities(nu, lam, rho_d_m, rho, eta, mfp, grav, v_f)
    implicit none

    real(dp), intent(in) :: nu, lam, rho_d_m, rho, eta, mfp, grav
    real(dp), dimension(:), intent(out) :: v_f

    integer :: iq
    real(dp) :: x, m, r, v
    real(dp) :: weight0, weight1, weight2
    real(dp) :: num0, num1, num2, den0, den1, den2

    num0 = 0.0_dp; num1 = 0.0_dp; num2 = 0.0_dp
    den0 = 0.0_dp; den1 = 0.0_dp; den2 = 0.0_dp

    do iq = 1, nq_gl
      x = lag_x(iq)
      m = max(lam*x, tiny_pos)
      r = max(((3.0_dp*m)/(4.0_dp*pi*rho_d_m))**third, r_seed)
      r = min(r, r_dyn_max)

      v = terminal_velocity_particle(r, rho_d_m, rho, eta, mfp, grav)

      weight0 = lag_w(iq) * exp((nu - 1.0_dp)*log(x))
      weight1 = weight0 * m
      weight2 = weight1 * m

      num0 = num0 + weight0*v
      den0 = den0 + weight0
      num1 = num1 + weight1*v
      den1 = den1 + weight1
      num2 = num2 + weight2*v
      den2 = den2 + weight2
    end do

    if (size(v_f) >= 1) v_f(1) = max(num0/max(den0, tiny_pos), v_floor)
    if (size(v_f) >= 2) v_f(2) = max(num1/max(den1, tiny_pos), v_floor)
    if (size(v_f) >= 3) v_f(3) = max(num2/max(den2, tiny_pos), v_floor)

  end subroutine gamma_moment_velocities

  function terminal_velocity_particle(r, rho_d_m, rho, eta, mfp, grav) result(v)
    implicit none

    real(dp), intent(in) :: r, rho_d_m, rho, eta, mfp, grav
    real(dp) :: v
    real(dp) :: Kn, beta, reynolds_fac, delta_rho

    Kn = mfp/max(r, r_seed)
    delta_rho = max(rho_d_m - rho, 0.0_dp)

    ! Cunningham slip factor, Jung et al. 2012
    beta = 1.0_dp + Kn*(1.165_dp + 0.480_dp*exp(-0.101_dp/max(Kn, tiny_pos)))
    reynolds_fac = (1.0_dp + &
      ((0.45_dp*grav*r**3*rho*rho_d_m)/max(54.0_dp*eta**2, tiny_pos))**0.4_dp)**(-1.25_dp)

    v = (2.0_dp*beta*grav*r**2*delta_rho)/max(9.0_dp*eta, tiny_pos) * reynolds_fac
    v = max(v, v_floor)

  end function terminal_velocity_particle

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
    character(len=20) :: sp

    allocate(d_g(n_bg), LJ_g(n_bg), molg_g(n_bg), eta_g(n_bg))

    do i = 1, n_bg
      sp = trim(adjustl(sp_bg(i)))
      select case(sp)
      case('OH')
        d_g(i) = d_OH; LJ_g(i) = LJ_OH; molg_g(i) = molg_OH
      case('H2')
        d_g(i) = d_H2; LJ_g(i) = LJ_H2; molg_g(i) = molg_H2
      case('H2O')
        d_g(i) = d_H2O; LJ_g(i) = LJ_H2O; molg_g(i) = molg_H2O
      case('H')
        d_g(i) = d_H; LJ_g(i) = LJ_H; molg_g(i) = molg_H
      case('CO')
        d_g(i) = d_CO; LJ_g(i) = LJ_CO; molg_g(i) = molg_CO
      case('CO2')
        d_g(i) = d_CO2; LJ_g(i) = LJ_CO2; molg_g(i) = molg_CO2
      case('O')
        d_g(i) = d_O; LJ_g(i) = LJ_O; molg_g(i) = molg_O
      case('CH4')
        d_g(i) = d_CH4; LJ_g(i) = LJ_CH4; molg_g(i) = molg_CH4
      case('C2H2')
        d_g(i) = d_C2H2; LJ_g(i) = LJ_C2H2; molg_g(i) = molg_C2H2
      case('NH3')
        d_g(i) = d_NH3; LJ_g(i) = LJ_NH3; molg_g(i) = molg_NH3
      case('N2')
        d_g(i) = d_N2; LJ_g(i) = LJ_N2; molg_g(i) = molg_N2
      case('HCN')
        d_g(i) = d_HCN; LJ_g(i) = LJ_HCN; molg_g(i) = molg_HCN
      case('He')
        d_g(i) = d_He; LJ_g(i) = LJ_He; molg_g(i) = molg_He
      case default
        print*, 'Background gas species data not found: ', trim(sp_bg(i)), 'STOP'
        stop
      end select
    end do

    do i = 1, n_bg
      eta_g(i) = (5.0_dp/16.0_dp) * (sqrt(pi*(molg_g(i)*amu)*kb*T)/(pi*d_g(i)**2)) &
        & * ((((kb*T)/LJ_g(i))**0.16_dp)/1.22_dp)
    end do

    bot = 0.0_dp
    do i = 1, n_bg
      bot = bot + VMR_bg(i)*sqrt(molg_g(i))
    end do
    if (bot <= 0.0_dp) then
      eta_out = 0.0_dp
      return
    end if

    y(:) = (VMR_bg(:)*sqrt(molg_g(:)))/bot

    eta_out = 0.0_dp
    do i = 1, n_bg
      do j = 1, n_bg
        Eij = ((2.0_dp*sqrt(molg_g(i)*molg_g(j)))/(molg_g(i) + molg_g(j)))**0.375_dp
        part = (y(i)*y(j))/(sqrt(eta_g(i)*eta_g(j))) * Eij
        eta_out = eta_out + part
      end do
    end do

    eta_out = 1.0_dp/eta_out

  end subroutine eta_construct

end module mini_cloud_vf_num_mod
