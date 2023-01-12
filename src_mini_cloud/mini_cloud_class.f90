module mini_cloud_class
  use mini_cloud_precision
  implicit none

  !! Common constants and geometric conversions
  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp)
  real(dp), parameter :: twopi = 2.0_dp*pi, fourpi = 4.0_dp*pi
  real(dp), parameter :: half = 1.0_dp/2.0_dp
  real(dp), parameter :: third = 1.0_dp/3.0_dp, twothird = 2.0_dp/3.0_dp
  real(dp), parameter :: fourpi3 = fourpi/3.0_dp

  !! Common physical constants - CODATA 2018
  real(dp), parameter :: kb = 1.380649e-16_dp ! erg K^-1 - Boltzmann's constant
  real(dp), parameter :: sb_c = 5.670374419e-5_dp ! erg cm-2 s-1 K^-4 - Stefan-Boltzmann's constant
  real(dp), parameter :: hp = 6.62607015e-27_dp ! erg s - Planck's constant
  real(dp), parameter :: c_s = 2.99792458e10_dp ! cm s^-1 - Vacuum speed of light
  real(dp), parameter :: amu = 1.66053906660e-24_dp ! g - Atomic mass unit
  real(dp), parameter :: N_A = 6.02214076e23_dp ! mol^-1 - Avogadro's constant
  real(dp), parameter :: gam = 1.4_dp
  real(dp), parameter :: Rgas = 8.31446261815324e7_dp, Rgas_si = 8.31446261815324_dp
  real(dp), parameter :: cal = 4.184e7_dp, cal_si = 4.184_dp ! Carlories
  real(dp), parameter :: mH = 1.007825 * amu  ! Mass of Hydrogen atom


  !! Common unit conversions
  real(dp), parameter :: bar = 1.0e6_dp ! bar to dyne
  real(dp), parameter :: atm = 1.01325e6_dp ! atm to dyne
  real(dp), parameter :: pa = 10.0_dp ! pa to dyne
  real(dp), parameter :: mmHg = 1333.2239_dp  ! mmHg to dyne



  real(dp) :: P_cgs, T, nd_atm
  real(dp) :: a_seed, Ar_seed, V_seed

  type dust
    integer :: idx
    character(len=11) :: name
    real(dp) :: bulk_den
    real(dp) :: mol_wght
    real(dp) :: mass
    real(dp) :: dV ! Volume of one unit

    real(dp) :: sat

    integer :: isig
    real(dp) :: sig

    integer :: inuc
    real(dp) :: Nl, Nf, a0, alpha
    real(dp) :: Vl, Vl13, Vl23
    real(dp) :: Js

    real(dp) :: v_stoi

    real(dp) :: chis, sevap

  end type dust

  type(dust), allocatable, dimension(:) :: d_sp



contains


  subroutine mini_cloud_init(n_dust, sp)

    integer, intent(in) :: n_dust
    character(len=11), dimension(n_dust), intent(in) :: sp

    integer :: n

    allocate(d_sp(n_dust))

    a_seed = 1.0e-7_dp 
    Ar_seed = fourpi * a_seed**2
    V_seed = fourpi3 * a_seed**3

    d_sp(:)%name = sp(:)
    

    !! Large loop to hardcoded values
    do n = 1, n_dust

      print*, trim(d_sp(n)%name)

      select case(trim(d_sp(n)%name))
      case('C')

        d_sp(n)%bulk_den = 2.27_dp
        d_sp(n)%mol_wght = 12.0107_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 1
        d_sp(n)%Nl = (4.0_dp/3.0_dp * pi * a_seed**3) / d_sp(n)%dV !1000.0_dp
        d_sp(n)%Nf = 5.0_dp
        d_sp(n)%a0 = 1.553_dp * 1e-8_dp
        d_sp(n)%alpha = 1.0_dp
        d_sp(n)%sig = 1200.0_dp
        d_sp(n)%v_stoi = 1.0_dp

      case('TiC')

        d_sp(n)%bulk_den = 4.93_dp
        d_sp(n)%mol_wght = 59.8777_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 0
        d_sp(n)%sig = 1000.0_dp
        d_sp(n)%v_stoi = 1.0_dp
        d_sp(n)%alpha = 1.0_dp

      case('SiC')

        d_sp(n)%bulk_den = 3.21_dp
        d_sp(n)%mol_wght = 40.0962_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 0
        d_sp(n)%sig = 1800.0_dp
        d_sp(n)%v_stoi = 1.0_dp
        d_sp(n)%alpha = 1.0_dp

      case('CaTiO3')

        d_sp(n)%bulk_den = 3.98_dp
        d_sp(n)%mol_wght = 135.943_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 0
        d_sp(n)%sig = 915.0_dp
        d_sp(n)%v_stoi = 1.0_dp
        d_sp(n)%alpha = 1.0_dp

      case ('TiO2')

        d_sp(n)%bulk_den = 4.23_dp
        d_sp(n)%mol_wght = 79.866_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 1
        d_sp(n)%Nl = (4.0_dp/3.0_dp * pi * a_seed**3) / d_sp(n)%dV ! 1000.0_dp
        d_sp(n)%Nf = 0.0_dp
        d_sp(n)%a0 = ((3.0_dp*d_sp(n)%dV)/(4.0_dp*pi))**(third)  !1.956_dp * 1e-8_dp !
        d_sp(n)%alpha = 1.0_dp
        d_sp(n)%v_stoi = 1.0_dp

      case ('VO')

        d_sp(n)%bulk_den = 5.76_dp
        d_sp(n)%mol_wght = 66.94090_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 0
        d_sp(n)%sig = 600.0_dp
        d_sp(n)%v_stoi = 1.0_dp
        d_sp(n)%alpha = 1.0_dp

      case ('Al2O3')

        d_sp(n)%bulk_den = 3.986_dp
        d_sp(n)%mol_wght = 101.961_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 0
        d_sp(n)%v_stoi = 2.0_dp
        d_sp(n)%alpha = 1.0_dp

      case ('Fe')

        d_sp(n)%bulk_den = 7.87_dp
        d_sp(n)%mol_wght = 55.845_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 0
        d_sp(n)%v_stoi = 1.0_dp
        d_sp(n)%alpha = 1.0_dp

      case('Mg2SiO4')

        d_sp(n)%bulk_den = 3.21_dp
        d_sp(n)%mol_wght = 140.693_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 0
        d_sp(n)%v_stoi = 2.0_dp
        d_sp(n)%alpha = 1.0_dp

      case('MgSiO3')

        d_sp(n)%bulk_den = 3.19_dp
        d_sp(n)%mol_wght = 100.389_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 0
        d_sp(n)%sig = 400.0_dp
        d_sp(n)%v_stoi = 1.0_dp
        d_sp(n)%alpha = 1.0_dp

      case ('SiO2')

        d_sp(n)%bulk_den = 2.648_dp
        d_sp(n)%mol_wght = 60.084_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 0
        d_sp(n)%v_stoi = 1.0_dp
        d_sp(n)%alpha = 1.0_dp

      case ('SiO')

        d_sp(n)%bulk_den = 2.18_dp
        d_sp(n)%mol_wght = 44.085_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 1
        d_sp(n)%Nl = (4.0_dp/3.0_dp * pi * a_seed**3) / d_sp(n)%dV !1000.0_dp
        d_sp(n)%Nf = 1.0_dp
        d_sp(n)%a0 = ((3.0_dp*d_sp(n)%dV)/(4.0_dp*pi))**(third) !2.001_dp * 1e-8_dp
        d_sp(n)%alpha = 1.0_dp
        d_sp(n)%sig = 500.0_dp
        d_sp(n)%v_stoi = 1.0_dp
        d_sp(n)%alpha = 1.0_dp

      case ('Cr')

        d_sp(n)%bulk_den = 7.19_dp
        d_sp(n)%mol_wght = 51.996_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 1
        d_sp(n)%Nl = (4.0_dp/3.0_dp * pi * a_seed**3) / d_sp(n)%dV !1000.0_dp
        d_sp(n)%Nf = 1.0_dp
        d_sp(n)%a0 = ((3.0_dp*d_sp(n)%dV)/(4.0_dp*pi))**(third) !1.421_dp * 1e-8_dp
        d_sp(n)%alpha = 1.0_dp
        d_sp(n)%v_stoi = 1.0_dp
        d_sp(n)%alpha = 1.0_dp

      case ('MnS')

        d_sp(n)%bulk_den = 4.08_dp
        d_sp(n)%mol_wght = 87.003_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 0
        d_sp(n)%sig = 2326.0_dp
        d_sp(n)%v_stoi = 1.0_dp
        d_sp(n)%alpha = 1.0_dp

      case ('Na2S')

        d_sp(n)%bulk_den = 1.856_dp
        d_sp(n)%mol_wght = 78.0445_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 0
        d_sp(n)%sig = 1033.0_dp
        d_sp(n)%v_stoi = 2.0_dp
        d_sp(n)%alpha = 1.0_dp

      case ('ZnS')

        d_sp(n)%bulk_den = 4.09_dp
        d_sp(n)%mol_wght = 97.445_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 0
        d_sp(n)%sig = 860.0_dp
        d_sp(n)%v_stoi = 1.0_dp
        d_sp(n)%alpha = 1.0_dp

      case ('KCl')

        d_sp(n)%bulk_den = 1.99_dp
        d_sp(n)%mol_wght = 74.551_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 1
        d_sp(n)%Nl = (4.0_dp/3.0_dp * pi * a_seed**3) / d_sp(n)%dV !1000.0_dp
        d_sp(n)%Nf = 1.0_dp
        d_sp(n)%a0 = ((3.0_dp*d_sp(n)%dV)/(4.0_dp*pi))**(third) !2.458_dp * 1e-8_dp
        d_sp(n)%alpha = 1.0_dp
        d_sp(n)%v_stoi = 1.0_dp
        d_sp(n)%alpha = 1.0_dp

      case ('NaCl')

        d_sp(n)%bulk_den = 2.165_dp
        d_sp(n)%mol_wght = 58.443_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 1
        d_sp(n)%Nl = (4.0_dp/3.0_dp * pi * a_seed**3) / d_sp(n)%dV !1000.0_dp
        d_sp(n)%Nf = 1.0_dp
        d_sp(n)%a0 = ((3.0_dp*d_sp(n)%dV)/(4.0_dp*pi))**(third) !2.205_dp * 1e-8_dp
        d_sp(n)%alpha = 1.0_dp
        d_sp(n)%v_stoi = 1.0_dp
        d_sp(n)%alpha = 1.0_dp

      case ('NH4Cl')

        d_sp(n)%bulk_den = 1.52_dp
        d_sp(n)%mol_wght = 53.4915_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 0
        d_sp(n)%sig = 56.0_dp
        d_sp(n)%v_stoi = 1.0_dp
        d_sp(n)%alpha = 1.0_dp

      case ('H2O')

        d_sp(n)%bulk_den = 0.93_dp
        d_sp(n)%mol_wght = 18.015_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 0
        d_sp(n)%v_stoi = 1.0_dp
        d_sp(n)%alpha = 1.0_dp

      case ('NH3')

        d_sp(n)%bulk_den = 0.87_dp
        d_sp(n)%mol_wght = 17.031_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 0
        d_sp(n)%sig = 23.4_dp
        d_sp(n)%v_stoi = 1.0_dp
        d_sp(n)%alpha = 1.0_dp

      case ('CH4')

        d_sp(n)%bulk_den = 0.656_dp
        d_sp(n)%mol_wght = 16.043_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 0
        d_sp(n)%sig = 14.0_dp
        d_sp(n)%v_stoi = 1.0_dp
        d_sp(n)%alpha = 1.0_dp

      case ('NH4SH')

        d_sp(n)%bulk_den = 1.17_dp
        d_sp(n)%mol_wght = 51.111_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 0
        d_sp(n)%sig = 50.0_dp
        d_sp(n)%v_stoi = 1.0_dp
        d_sp(n)%alpha = 1.0_dp

      case ('H2S')

        d_sp(n)%bulk_den = 1.12_dp
        d_sp(n)%mol_wght = 34.08_dp
        d_sp(n)%mass = d_sp(n)%mol_wght * amu
        d_sp(n)%dV = d_sp(n)%mass / d_sp(n)%bulk_den
        d_sp(n)%inuc = 0
        d_sp(n)%sig = 58.1_dp
        d_sp(n)%v_stoi = 1.0_dp
        d_sp(n)%alpha = 1.0_dp

      case default
        print*, 'Class: Species not included yet: ', d_sp(:)%name, 'stopping'
        stop
      end select

    end do

  end subroutine mini_cloud_init

end module mini_cloud_class
