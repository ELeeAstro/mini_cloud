module mini_cloud_class_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  ! Fortran 2008 intrinsic precisions - reccomended if possible
  !integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  !integer, parameter :: qp = REAL128

  !! Constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: third = 1.0_dp/3.0_dp
  real(dp), parameter :: twothird = 2.0_dp/3.0_dp

  !! Common physical constants - CODATA 2018
  real(dp), parameter :: kb = 1.380649e-16_dp ! erg K^-1 - Boltzmann's constant
  real(dp), parameter :: sb_c = 5.670374419e-5_dp ! erg cm-2 s-1 K^-4 - Stefan-Boltzmann's constant
  real(dp), parameter :: hp = 6.62607015e-27_dp ! erg s - Planck's constant
  real(dp), parameter :: c_s = 2.99792458e10_dp ! cm s^-1 - Vacuum speed of light
  real(dp), parameter :: amu = 1.66053906660e-24_dp ! g - Atomic mass unit
  real(dp), parameter :: N_A = 6.02214076e23_dp ! mol^-1 - Avogadro's constant
  real(dp), parameter :: eV = 1.60217663e-12_dp
  real(dp), parameter :: gam = 1.4_dp
  real(dp), parameter :: Rgas = 8.31446261815324e7_dp, Rgas_si = 8.31446261815324_dp
  real(dp), parameter :: cal = 4.184e7_dp, cal_si = 4.184_dp ! Carlories
  real(dp), parameter :: mH = 1.007825 * amu  ! Mass of Hydrogen atom

  !! Common unit conversions
  real(dp), parameter :: bar = 1.0e6_dp ! bar to dyne
  real(dp), parameter :: atm = 1.01325e6_dp ! atm to dyne
  real(dp), parameter :: pa = 10.0_dp ! pa to dyne
  real(dp), parameter :: mmHg = 1333.2239_dp  ! mmHg to dyne


  !! Turbulent scales
  real(dp), parameter :: l_k = 1.0_dp

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

  type dust

    integer :: idx
    character(len=20) :: name

    real(dp) :: rho_s
    real(dp) :: mw

    real(dp) :: m0, r0, V0, d0

    real(dp) :: r_seed

    real(dp) :: p_vap
    real(dp) :: sat
    real(dp) :: vth, D

    real(dp) :: sig

    integer :: inuc
    real(dp) :: Nl, Nf, a0, alpha
    real(dp) :: Vl, Vl13, Vl23
    real(dp) :: Js

  end type dust

  type(dust), allocatable, dimension(:) :: d_sp
  !$omp threadprivate (d_sp)

contains

  subroutine mini_cloud_init(n_dust, sp)

    integer, intent(in) :: n_dust
    character(len=20), dimension(n_dust), intent(in) :: sp

    integer :: n

    allocate(d_sp(n_dust))

    !! Large loop to hardcoded values
    do n = 1, n_dust

      print*, trim(sp(n))

      select case(trim(sp(n)))


      case ('TiO2')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho_s = 4.23_dp
        d_sp(n)%mw = 79.866_dp
        d_sp(n)%m0 = d_sp(n)%mw * amu
        d_sp(n)%V0 = d_sp(n)%m0 / d_sp(n)%rho_s
        d_sp(n)%r0 = ((3.0_dp*d_sp(n)%V0)/(4.0_dp*pi))**(third)
        d_sp(n)%d0 = 2.0_dp * r0 

        d_sp(n)%inuc = 1
        d_sp(n)%r_seed = 1e-7_dp
        d_sp(n)%Nl = (4.0_dp/3.0_dp * pi * d_sp(n)%r_seed**3) / d_sp(n)%V0
        d_sp(n)%Nf = 0.0_dp

      case ('Al2O3')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho_s = 3.986_dp
        d_sp(n)%mw = 101.961_dp
        d_sp(n)%m0 = d_sp(n)%mw * amu
        d_sp(n)%V0 = d_sp(n)%m0 / d_sp(n)%rho_s
        d_sp(n)%r0 = ((3.0_dp*d_sp(n)%V0)/(4.0_dp*pi))**(third)
        d_sp(n)%d0 = 2.0_dp * r0 

        d_sp(n)%inuc = 0

      case ('Fe')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho_s = 7.87_dp
        d_sp(n)%mw = 55.845_dp
        d_sp(n)%m0 = d_sp(n)%mw * amu
        d_sp(n)%V0 = d_sp(n)%m0 / d_sp(n)%rho_s
        d_sp(n)%r0 = ((3.0_dp*d_sp(n)%V0)/(4.0_dp*pi))**(third)
        d_sp(n)%d0 = 2.0_dp * r0 

        d_sp(n)%inuc = 0

      case('MgSiO3')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho_s = 3.19_dp
        d_sp(n)%mw = 100.389_dp
        d_sp(n)%m0 = d_sp(n)%mw * amu
        d_sp(n)%V0 = d_sp(n)%m0 / d_sp(n)%rho_s
        d_sp(n)%r0 = ((3.0_dp*d_sp(n)%V0)/(4.0_dp*pi))**(third)
        d_sp(n)%d0 = 2.0_dp * r0 

        d_sp(n)%inuc = 0

      case default
        print*, 'Class: Species not included yet: ', d_sp(:)%name, 'stopping'
        stop
      end select

    end do

  end subroutine mini_cloud_init

  !! Vapour pressure for each species
  subroutine p_vap_sp(n_dust,T)
    implicit none

    integer, intent(in) :: n_dust
    real(dp), intent(in) :: T

    real(dp) :: TC

    TC = T - 273.15_dp

    do n = 1, n_dust

      ! Return vapour pressure in dyne
      select case(trim(d_sp(n)%name))
      case('TiO2')
        ! Woitke & Helling (2004)
        d_sp(n)%p_vap = exp(35.8027_dp - 74734.7_dp/T)
      case('Al2O3')
        ! Kozasa et al. (1989)
        d_sp(n)%p_vap = exp(-73503.0_dp/T + 22.01_dp) * atm
      case('Fe')
        ! Elspeth note: Changed to Ackerman & Marley et al. (2001) expression
        if (T > 1800.0_dp) then
          d_sp(n)%p_vap = exp(9.86_dp - 37120.0_dp/T) * bar
        else
          d_sp(n)%p_vap = exp(15.71_dp - 47664.0_dp/T) * bar
        end if
      case('MgSiO3','MgSiO3_amorph')
        ! Ackerman & Marley (2001)
        d_sp(n)%p_vap = exp(-58663.0_dp/T + 25.37_dp) * bar
      case default
        print*, 'Vapour pressure species not found: ', trim(d_sp(n)%name), 'STOP'
        stop
      end select

    end do

  end subroutine p_vap_sp

  !! Surface tension for each species
  subroutine surface_tension(n_dust, T)
    implicit none

    integer, intent(in) :: n_dust
    real(dp), intent(in) :: T

    real(dp) :: TC

    TC = T - 273.15_dp

    do n = 1, n_dust

      ! Return surface tension in erg cm-2
      select case (trim(d_sp(n)%name))
      case('TiO2')
        ! Sindel et al. (2022)
        d_sp(n)%sig = 589.79_dp - 0.0708_dp * T
      case('Fe')
        ! http://www.kayelaby.npl.co.uk/general_physics/2_2/2_2_5.html
        d_sp(n)%sig = 1862.0_dp - 0.39_dp * (TC - 1530.0_dp)
      case('Al2O3')
        ! Pradhan et al. (2009)
        d_sp(n)%sig = 1024.0_dp - 0.177_dp * T
      case('MgSiO3')
        ! Janz 1988
        d_sp(n)%sig = 197.3_dp + 0.098_dp * T
      case default
        print*, 'Species surface tension not found: ', trim(d_sp(n)%name), 'STOP'
        stop
      end select

      d_sp(n)%sig = max(10.0_dp, d_sp(n)%sig)

    end do

  end subroutine surface_tension

  !! eta for background gas
  subroutine eta_construct(n_bg, sp_bg)
    implicit none

    integer, intent(in) :: n_bg
    character(len=20), dimension(:), intent(in) :: sp_bg
    
    integer :: i

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

end module mini_cloud_class
