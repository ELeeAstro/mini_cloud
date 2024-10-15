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
  real(dp), parameter :: amu = 1.66053906892e-24_dp ! g - Atomic mass unit
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
  real(dp), parameter :: mmHg = 1333.22387415_dp  ! mmHg to dyne

  !! Seed particle size
  real(dp), parameter :: r_seed = 1e-7_dp
  real(dp) :: V_seed, m_seed

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
  real(dp), allocatable, dimension(:) :: d_g, LJ_g, molg_g, eta_g

  type dust

    integer :: idx
    character(len=20) :: name

    real(dp) :: rho
    real(dp) :: mw
    real(dp) :: chi_v

    real(dp) :: m0, r0, V0, d0

    real(dp) :: p_vap
    real(dp) :: sat
    real(dp) :: vth, D

    real(dp) :: sig, L

    integer :: inuc
    real(dp) :: Nl, Nf, a0, al

  end type dust

  type gas

    integer :: idx
    character(len=20) :: name

    real(dp) :: mw, chi

  end type gas

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

      case('C')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 2.5_dp
        d_sp(n)%mw = 12.01070_dp

        d_sp(n)%L = 41523.0_dp * log(10.0_dp) * Rgas / d_sp(n)%mw 

        d_sp(n)%inuc = 1
        d_sp(n)%Nf = 5.0_dp
        d_sp(n)%al = 0.39_dp

        V_seed = 4.0_dp/3.0_dp * pi * r_seed**3
        m_seed = V_seed * d_sp(n)%rho

      case('TiC')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 4.93_dp
        d_sp(n)%mw = 59.8777_dp

        d_sp(n)%L = 33600.0_dp * log(10.0_dp) * Rgas / d_sp(n)%mw 

        d_sp(n)%inuc = 1
        d_sp(n)%Nf = 5.0_dp
        d_sp(n)%al = 0.29_dp

        V_seed = 4.0_dp/3.0_dp * pi * r_seed**3
        m_seed = V_seed * d_sp(n)%rho

      case('SiC')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 3.21_dp
        d_sp(n)%mw = 40.0962_dp

        d_sp(n)%L = 9.51431385e4_dp * Rgas / d_sp(n)%mw 

        d_sp(n)%inuc = 0
        
      case('CaTiO3')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 3.91_dp
        d_sp(n)%mw = 135.9432_dp

        d_sp(n)%L = 79568.2_dp * Rgas / d_sp(n)%mw 

        d_sp(n)%inuc = 0
        
      case ('Al2O3')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 3.97_dp
        d_sp(n)%mw = 101.961_dp

        d_sp(n)%chi_v = 2.0_dp

        d_sp(n)%L = 73503.0_dp * Rgas / d_sp(n)%mw 

        d_sp(n)%inuc = 0

      case('TiO2')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 4.23_dp
        d_sp(n)%mw = 79.866_dp

        d_sp(n)%chi_v = 1.0_dp

        d_sp(n)%L = 7.70443e4_dp * Rgas / d_sp(n)%mw 

        d_sp(n)%inuc = 1 
        d_sp(n)%Nf = 0.0_dp
        d_sp(n)%al = 0.2_dp

        V_seed = 4.0_dp/3.0_dp * pi * r_seed**3
        m_seed = V_seed * d_sp(n)%rho

      case('VO')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 5.76_dp
        d_sp(n)%mw = 66.94090_dp

        d_sp(n)%L = 6.74603e4_dp * Rgas / d_sp(n)%mw 

        d_sp(n)%inuc = 0

      case('Fe')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 7.87_dp
        d_sp(n)%mw = 55.845_dp

        d_sp(n)%chi_v = 1.0_dp


        d_sp(n)%L = 37120.0_dp * Rgas / d_sp(n)%mw 

        d_sp(n)%inuc = 0

      case('FeS')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 4.83_dp
        d_sp(n)%mw = 87.91_dp

        d_sp(n)%L = 5.69922e4_dp * Rgas / d_sp(n)%mw 

        d_sp(n)%inuc = 0

      case('FeO')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 5.99_dp
        d_sp(n)%mw = 71.8444_dp

        d_sp(n)%chi_v = 1.0_dp

        d_sp(n)%L = 6.30018e4_dp * Rgas / d_sp(n)%mw 

        d_sp(n)%inuc = 0

      case('Mg2SiO4')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 3.21_dp
        d_sp(n)%mw = 140.6931_dp

        d_sp(n)%chi_v = 2.0_dp


        d_sp(n)%L = 62279.0_dp * Rgas / d_sp(n)%mw 

        d_sp(n)%inuc = 0

      case('MgSiO3')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 3.19_dp
        d_sp(n)%mw = 100.389_dp

        d_sp(n)%chi_v = 1.0_dp
        
        d_sp(n)%L = 58663.0_dp * Rgas / d_sp(n)%mw 

        d_sp(n)%inuc = 0

      case('MgS')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 2.68_dp
        d_sp(n)%mw = 56.3700_dp

        d_sp(n)%L = 5.92010440e4_dp * Rgas / d_sp(n)%mw 

        d_sp(n)%inuc = 0       

      case('MgO')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 3.58_dp
        d_sp(n)%mw = 40.30440_dp

        d_sp(n)%L = 7.91838e4_dp * Rgas / d_sp(n)%mw 

        d_sp(n)%inuc = 0

      case('SiO2')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 2.648_dp
        d_sp(n)%mw = 60.08430_dp

        d_sp(n)%L = 7.28086e4_dp * Rgas / d_sp(n)%mw 

        d_sp(n)%inuc = 0

      case('SiO')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 2.18_dp
        d_sp(n)%mw = 44.08490_dp

        d_sp(n)%L = 49520.0_dp * Rgas / d_sp(n)%mw 

        d_sp(n)%inuc = 1

        d_sp(n)%Nf = 5.0_dp
        d_sp(n)%al = 0.2_dp

        V_seed = 4.0_dp/3.0_dp * pi * r_seed**3
        m_seed = V_seed * d_sp(n)%rho

      case('Cr')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 7.19_dp
        d_sp(n)%mw = 51.99610_dp

        d_sp(n)%L = 4.78455e4_dp * Rgas / d_sp(n)%mw 

        d_sp(n)%inuc = 0

      case('MnS')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 4.08_dp
        d_sp(n)%mw = 87.0030_dp

        d_sp(n)%L = 23810.0_dp * log(10.0_dp) * Rgas / d_sp(n)%mw 

        d_sp(n)%inuc = 0

      case('Na2S')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 1.86_dp
        d_sp(n)%mw = 78.0445_dp

        d_sp(n)%L = 13889.0_dp * log(10.0_dp) * Rgas / d_sp(n)%mw 

        d_sp(n)%inuc = 0        
     
      case('ZnS')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 4.09_dp
        d_sp(n)%mw = 97.4450_dp

        d_sp(n)%L = 4.75507888e4_dp * Rgas / d_sp(n)%mw 

        d_sp(n)%inuc = 0  
        
      case('KCl')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 1.99_dp
        d_sp(n)%mw = 74.551_dp

        d_sp(n)%L = 2.69250e4_dp * Rgas / d_sp(n)%mw

        d_sp(n)%inuc = 1
        d_sp(n)%Nf = 5.0_dp
        d_sp(n)%al = 1.0_dp

        V_seed = 4.0_dp/3.0_dp * pi * r_seed**3
        m_seed = V_seed * d_sp(n)%rho

      case('NaCl')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 2.165_dp
        d_sp(n)%mw = 58.4428_dp

        d_sp(n)%L = 2.79146e4_dp * Rgas / d_sp(n)%mw

        d_sp(n)%inuc = 1
        d_sp(n)%Nf = 5.0_dp
        d_sp(n)%al = 1.0_dp

        V_seed = 4.0_dp/3.0_dp * pi * r_seed**3
        m_seed = V_seed * d_sp(n)%rho

      case('NH4Cl')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 1.52_dp
        d_sp(n)%mw = 53.4915_dp

        d_sp(n)%L = 4302.0_dp * log(10.0_dp) * Rgas / d_sp(n)%mw 

        d_sp(n)%inuc = 0  

      case('H2O')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 0.93_dp
        d_sp(n)%mw = 18.01528_dp

        d_sp(n)%L = 2257.0e7_dp ! Vapourisation assumed

        d_sp(n)%inuc = 0  

      case('NH3')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 0.87_dp
        d_sp(n)%mw = 17.03052_dp

        d_sp(n)%L = 1371.0e7_dp ! Vapourisation assumed

        d_sp(n)%inuc = 0 

      case('CH4')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 0.656_dp
        d_sp(n)%mw = 16.0425_dp

        d_sp(n)%L = 480.6e7_dp ! Vapourisation assumed

        d_sp(n)%inuc = 0

      case('NH4SH')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 1.17_dp
        d_sp(n)%mw = 51.1114_dp

        d_sp(n)%L = 2409.4_dp * log(10.0_dp) * Rgas / d_sp(n)%mw

        d_sp(n)%inuc = 0

      case('H2S')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 1.12_dp
        d_sp(n)%mw = 34.0809_dp

        d_sp(n)%L = 958.587_dp * log(10.0_dp) * Rgas / d_sp(n)%mw

        d_sp(n)%inuc = 0

      case('S2')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 2.07_dp
        d_sp(n)%mw = 64.1300_dp

        d_sp(n)%L = 14000.0_dp * Rgas / d_sp(n)%mw

        d_sp(n)%inuc = 0

      case('S8')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 2.07_dp
        d_sp(n)%mw = 256.5200_dp

        d_sp(n)%L = 7510.0_dp * Rgas / d_sp(n)%mw

        d_sp(n)%inuc = 0

      case('CO')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 1.14_dp
        d_sp(n)%mw = 28.0101_dp

        d_sp(n)%L = 7.8824e2_dp * log(10.0_dp) * Rgas / d_sp(n)%mw

        d_sp(n)%inuc = 0   

      case('CO2')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 1.98_dp
        d_sp(n)%mw = 44.0095_dp

        d_sp(n)%L = 1.5119e3_dp * log(10.0_dp) * Rgas / d_sp(n)%mw

        d_sp(n)%inuc = 0

      case('H2SO4')

        d_sp(n)%idx = n
        d_sp(n)%name = sp(n)

        d_sp(n)%rho = 1.6_dp
        d_sp(n)%mw = 98.0785_dp

        d_sp(n)%L = 1.01294e4_dp * Rgas / d_sp(n)%mw

        d_sp(n)%inuc = 0

      case default
        print*, 'Class: Species not included yet: ', trim(d_sp(n)%name), 'stopping'
        stop
      end select

      !! Variables common to all material
      d_sp(n)%m0 = d_sp(n)%mw * amu
      d_sp(n)%V0 = d_sp(n)%m0 / d_sp(n)%rho
      d_sp(n)%r0 = ((3.0_dp*d_sp(n)%V0)/(4.0_dp*pi))**(third)
      d_sp(n)%d0 = 2.0_dp * d_sp(n)%r0 

      if (d_sp(n)%inuc == 1) then
        d_sp(n)%Nl = (4.0_dp/3.0_dp * pi * r_seed**3) / d_sp(n)%V0
      end if

    end do

  end subroutine mini_cloud_init

  !! Vapour pressure for each species
  subroutine p_vap_sp(n_dust,T)
    implicit none

    integer, intent(in) :: n_dust
    real(dp), intent(in) :: T

    integer :: n
    real(dp) :: TC

    TC = T - 273.15_dp

    do n = 1, n_dust

      ! Return vapour pressure in dyne
      select case(trim(d_sp(n)%name))
      case('C')
        !d_sp(n)%p_vap = exp(3.27860e1_dp - 8.65139e4_dp/(T + 4.80395e-1_dp))
        ! Kimura et al. (2023)
        d_sp(n)%p_vap = 10.0_dp**(-41523.0_dp/T + 10.609_dp) * atm
      case('TiC')
        ! Kimura et al. (2023)
        d_sp(n)%p_vap = 10.0_dp**(-33600.0_dp/T + 7.652_dp) * atm
      case('SiC')
        ! Elspeth 5 polynomial JANAF-NIST fit
        d_sp(n)%p_vap = exp(-9.51431385e4_dp/T + 3.72019157e1_dp + 1.09809718e-3_dp*T &
          & -5.63629542e-7_dp*T**2 + 6.97886017e-11_dp*T**3)
      case('CaTiO3')
        ! Kozasa et al. (1987)
        d_sp(n)%p_vap = exp(-79568.2_dp/T + 42.0204_dp) * atm        
      case('Al2O3')
        ! Kozasa et al. (1989)
        d_sp(n)%p_vap = exp(-73503.0_dp/T + 22.01_dp) * atm
      case('TiO2')
        ! GGChem 5 polynomial NIST fit
        d_sp(n)%p_vap = exp(-7.70443e4_dp/T +  4.03144e1_dp - 2.59140e-3_dp*T &
          &  + 6.02422e-7_dp*T**2 - 6.86899e-11_dp*T**3)
      case('VO')
        ! NIST 5 param fit
        d_sp(n)%p_vap = exp(-6.74603e4_dp/T + 3.82717e1_dp - 2.78551e-3_dp*T &
          & + 5.72078e-7_dp*T**2 - 7.41840e-11_dp*T**3)        
      case('Fe')
        ! Elspeth note: Changed to Ackerman & Marley (2001) expression
        if (T > 1800.0_dp) then
          d_sp(n)%p_vap = exp(9.86_dp - 37120.0_dp/T) * bar
        else
          d_sp(n)%p_vap = exp(15.71_dp - 47664.0_dp/T) * bar
        end if
      case('FeS')
        ! GGChem 5 polynomial NIST fit
        d_sp(n)%p_vap = exp(-5.69922e4_dp/T + 3.86753e1_dp - 4.68301e-3_dp*T &
          & + 1.03559e-6_dp*T**2 - 8.42872e-11_dp*T**3)
      case('FeO')
        ! GGChem 5 polynomial NIST fit
        d_sp(n)%p_vap = exp(-6.30018e4_dp/T + 3.66364e1_dp - 2.42990e-3_dp*T &
          & + 3.18636e-7_dp*T**2)
      case('MgS')
        ! Elspeth 5 polynomial JANAF-NIST fit
        d_sp(n)%p_vap = exp(-5.92010440e4_dp/T +  3.58148138e1_dp - 1.93822353e-3_dp*T &
          &  + 9.77971406e-7_dp*T**2 - 11.74813601e-10_dp*T**3)
      case('Mg2SiO4')
        ! Kozasa et al. (1989)
        d_sp(n)%p_vap = exp(-62279.0_dp/T + 20.944_dp) * atm          
      case('MgSiO3','MgSiO3_amorph')
        ! Ackerman & Marley (2001)
        d_sp(n)%p_vap = exp(-58663.0_dp/T + 25.37_dp) * bar
      case('MgO')
        ! GGChem 5 polynomial NIST fit
        d_sp(n)%p_vap = exp(-7.91838e4_dp/T + 3.57312e1_dp + 1.45021e-4_dp*T &
          & - 8.47194e-8_dp*T**2 + 4.49221e-12_dp*T**3)
      case('SiO2','SiO2_amorph')
        ! GGChem 5 polynomial NIST fit
        d_sp(n)%p_vap = exp(-7.28086e4_dp/T + 3.65312e1_dp - 2.56109e-4_dp*T &
          & - 5.24980e-7_dp*T**2 + 1.53343E-10_dp*T**3) 
      case('SiO')
        ! Gail et al. (2013)
        d_sp(n)%p_vap = exp(-49520.0_dp/T + 32.52_dp)
      case('Cr')
        ! GGChem 5 polynomial NIST fit
        d_sp(n)%p_vap = exp(-4.78455e4_dp/T + 3.22423e1_dp - 5.28710e-4_dp*T & 
          &  - 6.17347e-8_dp*T**2 + 2.88469e-12_dp*T**3)
      case('MnS')
        ! Morley et al. (2012)
        d_sp(n)%p_vap = 10.0_dp**(11.532_dp - 23810.0_dp/T) * bar
      case('Na2S')
        ! Morley et al. (2012)
        d_sp(n)%p_vap =  10.0_dp**(8.550_dp - 13889.0_dp/T) * bar
      case('ZnS')
        ! Elspeth 5 polynomial Barin data fit
        d_sp(n)%p_vap = exp(-4.75507888e4_dp/T + 3.66993865e1_dp - 2.49490016e-3_dp*T &
          &  + 7.29116854e-7_dp*T**2 - 1.12734453e-10_dp*T**3)
        ! Morley et al. (2012)
        !d_sp(n)%p_vap = 10.0_dp**(12.812_dp - 15873.0_dp/T) * bar                   
      case('KCl')
        ! GGChem 5 polynomial NIST fit
        d_sp(n)%p_vap = exp(-2.69250e4_dp/T + 3.39574e1_dp - 2.04903e-3_dp*T &
          & -2.83957e-7_dp*T**2 + 1.82974e-10_dp*T**3)
      case('NaCl')
        ! GGChem 5 polynomial NIST fit
        d_sp(n)%p_vap = exp(-2.79146e4_dp/T + 3.46023e1_dp - 3.11287e3_dp*T & 
        & + 5.30965e-7_dp*T**2 -2.59584e-12_dp*T**3)
      case('NH4Cl')
        d_sp(n)%p_vap = 10.0_dp**(7.0220_dp - 4302.0_dp/T) * bar   
      case('H2O')
        ! Ackerman & Marley (2001) H2O liquid & ice vapour pressure expressions
        if (T > 1048.0_dp) then
          d_sp(n)%p_vap = 6.0e8_dp
        else if (T < 273.16_dp) then
          d_sp(n)%p_vap = 6111.5_dp * exp((23.036_dp * TC - TC**2/333.7_dp)/(TC + 279.82_dp))
        else
          d_sp(n)%p_vap = 6112.1_dp * exp((18.729_dp * TC - TC**2/227.3_dp)/(TC + 257.87_dp)) 
        end if
      case('NH3')
        ! Ackerman & Marley (2001) NH3 ice vapour pressure expression from Weast (1971)
        d_sp(n)%p_vap = exp(10.53_dp - 2161.0_dp/T - 86596.0_dp/T**2)  * bar
      case('CH4')
        !--- Prydz, R.; Goodwin, R.D., J. Chem. Thermodyn., 1972, 4,1 ---
        if (T < 0.5_dp) then ! Limiter for very very very cold T
          d_sp(n)%p_vap = 10.0_dp**(3.9895_dp - 443.028_dp/(0.5_dp-0.49_dp)) * bar       
        else
          d_sp(n)%p_vap = 10.0_dp**(3.9895_dp - 443.028_dp/(T-0.49_dp)) * bar
        end if
      case('NH4SH')
        !--- E.Lee's fit to Walker & Lumsden (1897) ---
        d_sp(n)%p_vap = 10.0_dp**(7.8974_dp - 2409.4_dp/T) * bar    
      case('H2S')
        !--- Stull (1947) ---
        if (T < 30.0_dp) then ! Limiter for very cold T
          d_sp(n)%p_vap = 10.0_dp**(4.43681_dp - 829.439_dp/(30.0_dp-25.412_dp)) * bar
        else if (T < 212.8_dp) then
          d_sp(n)%p_vap = 10.0_dp**(4.43681_dp - 829.439_dp/(T-25.412_dp)) * bar
        else
          d_sp(n)%p_vap = 10.0_dp**(4.52887_dp - 958.587_dp/(T-0.539_dp)) * bar
        end if
      case('S2')
        !--- Zahnle et al. (2016) ---
        if (T < 413.0_dp) then
          d_sp(n)%p_vap = exp(27.0_dp - 18500.0_dp/T) * bar
        else
          d_sp(n)%p_vap = exp(16.1_dp - 14000.0_dp/T) * bar
        end if
      case('S8')
        !--- Zahnle et al. (2016) ---
        if (T < 413.0_dp) then
          d_sp(n)%p_vap = exp(20.0_dp - 11800.0_dp/T) * bar
        else
          d_sp(n)%p_vap = exp(9.6_dp - 7510.0_dp/T) * bar
        end if
      case('CO')
        ! Yaws
        d_sp(n)%p_vap = 10.0_dp**(51.8145e0_dp - 7.8824e2_dp/T - 2.2734e1_dp*log10(T) &
          & + 5.1225e-2_dp*T + 4.6603e-11_dp*T**2) * mmHg
      case('CO2')
        ! Yaws
        d_sp(n)%p_vap = 10.0_dp**(35.0187e0_dp - 1.5119e3_dp/T - 1.1335e1_dp*log10(T) &
          & + 9.3383e-3_dp*T + 7.7626e-10_dp*T**2) * mmHg
      case('H2SO4')
        ! GGChem 5 polynomial NIST fit
        d_sp(n)%p_vap = exp(-1.01294e4_dp/T + 3.55465e1_dp - 8.34848e-3_dp*T)                               
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

    integer :: n
    real(dp) :: TC

    TC = T - 273.15_dp

    do n = 1, n_dust

      ! Return surface tension in erg cm-2
      select case (trim(d_sp(n)%name))
      case('C')
        ! Tabak et al. (1995)
        d_sp(n)%sig = 4200.0_dp
      case('TiC')
        ! Chigai et al. (1999)
        d_sp(n)%sig = 2700.0_dp
      case ('SiC')
        ! Nozawa et al. (2003)
        d_sp(n)%sig = 1800.0_dp   
      case('CaTiO3')
        ! Kozasa et al. (1987)
        d_sp(n)%sig = 494.0_dp
      case('Al2O3')
        ! Pradhan et al. (2009)
        d_sp(n)%sig = 1024.0_dp - 0.177_dp * T                     
      case('TiO2')
        ! Sindel et al. (2022)
        d_sp(n)%sig = 589.79_dp - 0.0708_dp * T
      case('Fe')
        ! http://www.kayelaby.npl.co.uk/general_physics/2_2/2_2_5.html
        d_sp(n)%sig = 1862.0_dp - 0.39_dp * (TC - 1530.0_dp)
      case('Fe2O3')
        d_sp(n)%sig = 410.0_dp
      case('FeO')
        ! Janz 1988
        d_sp(n)%sig = 585.0_dp
      case('FeS')
        d_sp(n)%sig = 380.0_dp
      case('MgS')
        ! Gail and Sedlmayr (1986)
        d_sp(n)%sig = 800.0_dp
      case('Mg2SiO4')
        ! Kozasa et al. (1989)
        d_sp(n)%sig = 436.0_dp           
      case('MgSiO3')
        ! Janz 1988
        d_sp(n)%sig = 197.3_dp + 0.098_dp * T
      case('MgO')
         ! Nozawa et al. (2003)
        d_sp(n)%sig = 1100.0_dp 
      case('SiO')
        ! Gail and Sedlmayr (1986)
        d_sp(n)%sig = 500.0_dp
      case('SiO2')
        ! Janz 1988
        d_sp(n)%sig = 243.2_dp + 0.031_dp * T 
      case('Cr')
        ! http://www.kayelaby.npl.co.uk/general_physics/2_2/2_2_5.html
        d_sp(n)%sig = 1642.0_dp - 0.20_dp * (TC - 1860.0_dp)      
      case('MnS')
        d_sp(n)%sig = 2326.0_dp
      case('Na2S')
        d_sp(n)%sig = 1033.0_dp  
      case('ZnS')
        d_sp(n)%sig = 860.0_dp                     
      case('KCl')
        ! Janz 1988
        d_sp(n)%sig = 175.57_dp - 0.07321_dp * T
      case('NaCl')
        ! Janz 1988
        d_sp(n)%sig = 191.16_dp - 0.07188_dp * T
      case('H2O')
        ! Hale and Plummer (1974)
        d_sp(n)%sig = 141.0_dp - 0.15_dp*TC
      case('NH3')
        ! Weast et al. (1988)
        d_sp(n)%sig = 23.4_dp
      case('NH4Cl')
        d_sp(n)%sig = 56.0_dp
      case('NH4SH')
        d_sp(n)%sig = 50.0_dp
      case('CH4')
        ! USCG (1984)
        d_sp(n)%sig = 14.0_dp
      case('H2S')
        ! Nehb and Vydra (2006)
        d_sp(n)%sig = 58.1_dp
      case('S2','S8')
        ! Fanelli 1950
        d_sp(n)%sig = 60.8_dp        
      case default
        print*, 'Species surface tension not found: ', trim(d_sp(n)%name), 'STOP'
        stop
      end select

      d_sp(n)%sig = max(10.0_dp, d_sp(n)%sig)

    end do

  end subroutine surface_tension

  !! eta for background gas
  subroutine eta_construct(n_bg, sp_bg, VMR_bg, T, eta_out)
    implicit none

    integer, intent(in) :: n_bg
    character(len=20), dimension(:), intent(in) :: sp_bg
    real(dp), dimension(n_bg), intent(in) :: VMR_bg
    real(dp), intent(in) :: T

    real(dp), intent(out) :: eta_out
    
    integer :: i, j, n
    real(dp) :: eta_sum, phi_ij_top, phi_ij_bot, phi_ij

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

    !! do Wilke (1950) classical mixing rule
    !! First calculate each species eta
    do n = 1, n_bg
      eta_g(n) = (5.0_dp/16.0_dp) * (sqrt(pi*(molg_g(n)*amu)*kb*T)/(pi*d_g(n)**2)) &
        & * ((((kb*T)/LJ_g(n))**(0.16_dp))/1.22_dp)
    end do

    !! Find weighting factor for each species
    eta_out = 0.0_dp
    do i = 1, n_bg
      eta_sum = 0.0_dp
      do j = 1, n_bg
        phi_ij_top = (1.0_dp + sqrt(eta_g(i)/eta_g(j)) * (molg_g(j)/molg_g(i))**(0.25_dp))**2
        phi_ij_bot = sqrt(8.0_dp*(1.0_dp + (molg_g(i)/molg_g(j))))
        phi_ij = phi_ij_top  / phi_ij_bot
        eta_sum = eta_sum + VMR_bg(j) * phi_ij
      end do
      eta_out = eta_out + (VMR_bg(i) * eta_g(i)) / eta_sum
    end do
    
  end subroutine eta_construct

  !! Thermal conductivity of background gass
  subroutine kappa_a_construct(n_bg, sp_bg, VMR_bg, T, kappa_out)
    implicit none

    integer, intent(in) :: n_bg
    character(len=20), dimension(:), intent(in) :: sp_bg
    real(dp), dimension(n_bg), intent(in) :: VMR_bg
    real(dp), intent(in) :: T

    real(dp), intent(out) :: kappa_out

    integer :: i, j, n
    real(dp), dimension(n_bg) :: kappa_g
    real(dp) :: top, bot, TTc

    real(dp) :: eta_sum, phi_ij_top, phi_ij_bot, phi_ij

    !! H2 parameters
    real(dp), dimension(7), parameter :: A1 = &
      & (/-3.40976e-1_dp, 4.58820e0_dp, -1.45080e0_dp, &
      &    3.26394e-1_dp, 3.16939e-3_dp, 1.90592e-4_dp, -1.13900e-6_dp/)
    real(dp), dimension(4), parameter :: A2 = &
      & (/1.38497e2_dp, -2.21878e1_dp, 4.57151e0_dp, 1.0e0_dp/)
    real(dp), parameter :: Tc = 33.145_dp

    !! He parameters
    real(dp), parameter :: A = 2.7870034e-3_dp, B = 7.034007057e-1_dp
    real(dp), dimension(4), parameter :: C = (/3.739232544_dp, &
      &  -2.620316969e1_dp, 5.982252246e1_dp, -4.26397634e1_dp/)

    !! H2O parameters
    real(dp), dimension(5), parameter :: L = (/2.443221e-3_dp, 1.323095e-2_dp, &
      & 6.770357e-3_dp, -3.454586e-3_dp, 4.096266e-4_dp /)
    real(dp), parameter :: Ts = 647.096_dp
    real(dp) :: Tbar

    !! CO parameters 
    real(dp), dimension(3), parameter :: AA = (/3.29558e-4_dp, 3.05976e-6_dp, -3.13222e-9_dp/)

    do n = 1, n_bg

      select case(sp_bg(n))

      case('H2')
        !! Pure Normal Hydrogen (H2) dilute limit fitting function from Assael et al. (2011):
        !! Correlation of the Thermal Conductivity of Normal and Parahydrogen from the Triple Point to 1000 K and up to 100 MPa

        TTc = T/Tc

        top = 0.0_dp
        do i = 1, 7
          top = top + A1(i)*TTc**(i-1)
        end do

        bot = 0.0_dp
        do i = 1, 4
          bot = bot + A2(i)*TTc**(i-1)
        end do

        !! Fitting function is in W m-1 K-1, convert to erg s-1 cm-1 K-1
        kappa_g(n) = top/bot * 1e5_dp 

      case('He')

        !! A correlation of thermal conductivity data for helium (Hands & Arp 1981)
        kappa_g(n) = 0.0_dp
        do i = 1, 4
          kappa_g(n) = kappa_g(n) + C(i)/T**i
        end do

        !! Fitting function is in W m-1 K-1, convert to erg s-1 cm-1 K-1
        kappa_g(n) = A * T**B * exp(kappa_g(n)) * 1e5_dp

      case('H2O')

        !! New International Formulation for the Thermal Conductivity of H2O (Huber et al. 2012)
        Tbar = T/Ts

        kappa_g(n) = 0.0_dp
        do i = 1, 5
          kappa_g(n) = kappa_g(n) + L(i)/Tbar**(i-1)
        end do

        !! Fitting function is in W m-1 K-1, convert to erg s-1 cm-1 K-1
        kappa_g(n) = sqrt(Tbar)/kappa_g(n) * 1e5_dp

      case('CO')

        !! Models for Viscosity, Thermal Conductivity, and Surface Tension of Selected Pure Fluids as Implemented in REFPROP v10.0
        !! Huber (2018) - NIST

        kappa_g(n) = 0.0_dp
        do i = 1, 3
          kappa_g(n) = kappa_g(n) + AA(i)*T**(i-1)
        end do

        !! Fitting function is in W m-1 K-1, convert to erg s-1 cm-1 K-1
        kappa_g(n) = kappa_g(n) * 1e5_dp

      case default
        !! No conductivity data avilable for background species, assume zero
        kappa_g(n) = 0.0_dp
      end select

    end do

    !! Find weighting factor for each species - following Wilke (1950) mixing rule
    kappa_out = 0.0_dp
    do i = 1, n_bg
      eta_sum = 0.0_dp
      do j = 1, n_bg
        phi_ij_top = (1.0_dp + sqrt(eta_g(i)/eta_g(j)) * (molg_g(j)/molg_g(i))**(0.25_dp))**2
        phi_ij_bot = sqrt(8.0_dp*(1.0_dp + (molg_g(i)/molg_g(j))))
        phi_ij = phi_ij_top  / phi_ij_bot
        eta_sum = eta_sum + VMR_bg(j) * phi_ij
      end do
      kappa_out = kappa_out + (VMR_bg(i) * kappa_g(i)) / eta_sum
    end do

  end subroutine kappa_a_construct

end module mini_cloud_class_mod
