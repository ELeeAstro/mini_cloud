module mini_cloud_opac_mie_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  use lxmie_mod, only : lxmie
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp) ! value of pi
  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: amu = 1.66053906660e-24_dp ! g mol-1 (note, has to be cgs g mol-1 !!!)

  type nk_table

    character(len=11) :: name
    character(len=50) :: fname

    real(dp), allocatable, dimension(:) :: n, k

  end type nk_table  

  character(len=50) :: p_2_nk

  type(nk_table), allocatable, dimension(:)  :: nk
  !$omp threadprivate (nk)

  logical :: first_call = .True.
  !$omp threadprivate (first_call)

  public :: opac_mie
  private :: locate, linear_log_interp, m2e, e2m, &
    & read_nk_tables, read_and_interp_nk_table

contains

  subroutine opac_mie(n_dust, sp, T_in, mu_in, P_in, q_c, rm, rho_d, sigma, n_wl, wl, k_ext, alb, gg)
    implicit none

    integer, intent(in) :: n_dust, n_wl
    character(len=11), dimension(n_dust), intent(in) :: sp
    real(dp), intent(in) :: T_in, mu_in, P_in
    real(dp), intent(in) :: q_c, rm, sigma, rho_d
    real(dp), dimension(n_wl), intent(in) :: wl

    real(dp), dimension(n_wl), intent(out) :: k_ext, alb, gg

    integer :: l, n
    real(dp) :: rho, amean, n_d
    real(dp) :: x, xsec, q_ext, q_sca, q_abs, g
    real(dp), dimension(n_dust) :: b_mix
    complex(dp) :: e_eff, e_eff0
    complex(dp) :: N_eff
    complex(dp), dimension(n_dust) :: N_inc, e_inc


    if (first_call .eqv. .True.) then
      p_2_nk = 'nk_tables/'
      call read_nk_tables(n_dust, sp, n_wl, wl)
      first_call = .False.
    end if

    if (q_c < 1e-10_dp) then
      k_ext(:) = 0.0_dp
      alb(:) = 0.0
      gg(:) = 0.0_dp
      return
    end if

    rho = (P_in * 10.0_dp * mu_in * amu)/(kb * T_in)

    amean = rm * exp(5.0_dp/2.0_dp * log(sigma)**2)! Find effective particle size [cm]

    n_d = (3.0_dp * q_c * rho)/(4.0_dp*pi*rho_d*rm**3) * exp(-9.0_dp/2.0_dp * log(sigma)**2)! Find particle number density

    xsec = pi * amean**2
    
    b_mix(:) = 1.0_dp

    do l = 1, n_wl

      !! Size parameter 
      x = (2.0_dp * pi * amean) / (wl(l) * 1e-4_dp)

      !! Refractive index
      !! Use Landau-Lifshitz-Looyenga (LLL) method
      e_eff0 = cmplx(0.0_dp,0.0_dp,dp)
      do n = 1, n_dust
        N_inc(n) = cmplx(nk(n)%n(l),nk(n)%k(l),dp)
        e_inc(n) = m2e(N_inc(n))
        e_eff0 = e_eff0 + b_mix(n)*e_inc(n)**(1.0_dp/3.0_dp)
      end do
      e_eff = e_eff0**3
      N_eff = e2m(e_eff)

      if (x < 0.01_dp) then
        !! Use Rayleigh approximation
        call rayleigh(x, N_eff, q_abs, q_sca, q_ext)
        g = 0.0_dp
      else if (x > 10.0_dp) then
        call madt(x, N_eff, q_abs, q_sca, q_ext, g)
      else
        !! Call LX-MIE with negative k value
        N_eff = cmplx(real(N_eff,dp),-aimag(N_eff),dp)
        call lxmie(N_eff, x, q_ext, q_sca, q_abs, g)
      end if

      !! Calculate the opacity, abledo and mean cosine angle (asymmetry factor)
      k_ext(l) = (q_ext * xsec * n_d)/rho
      alb(l) = max(q_sca/q_ext, 0.95_dp)
      gg(l) = max(g, 0.0_dp)

    end do

  end subroutine opac_mie

  ! Small dielectric sphere approximation, small particle limit x << 1
  subroutine rayleigh(x, ri, q_abs, q_sca, q_ext)
    implicit none

    real(dp), intent(in) :: x
    complex(dp), intent(in) :: ri

    real(dp), intent(out) :: q_abs, q_sca, q_ext

    complex(dp) :: alp

    alp = (ri**2 - 1.0_dp)/(ri**2 + 2.0_dp)

    q_sca = 8.0_dp/3.0_dp * x**4 * abs(alp)**2
    q_ext = 4.0_dp * x * &
     & aimag(alp * (1.0_dp + x**2/15.0_dp*alp * ((ri**4+27.0_dp*ri**2+38.0_dp)/(2.0_dp*ri**2+3.0_dp)))) + &
     & 8.0_dp/3.0_dp * x**4 * real(alp**2,dp)

    q_abs = q_ext - q_sca

  end subroutine rayleigh
  
  subroutine madt(x, ri, q_abs, q_sca, q_ext, g)
    implicit none

    real(dp), intent(in) :: x
    complex(dp), intent(in) :: ri

    real(dp), intent(out) :: q_abs, q_sca, q_ext, g

    real(dp) :: n, k, rho, beta, tan_b
    real(dp) :: C1, C2, eps, q_edge, Cm

    n = real(ri,dp)
    k = aimag(ri)

    rho = 2.0_dp*x*(n - 1.0_dp)
    beta = atan(k/(n - 1.0_dp))
    tan_b = tan(beta)

    if (k < 1.0e-20_dp) then
      q_sca = 2.0_dp - (4.0_dp/rho)*sin(rho) - (4.0_dp/rho**2)*(cos(rho) - 1.0_dp)
      q_abs = 0.0_dp
      q_ext = q_sca
    else
      q_ext = 2.0_dp - 4.0_dp * exp(-rho*tan_b)*(cos(beta)/rho)*sin(rho-beta) &
        & - 4.0_dp*exp(-rho*tan_b)*(cos(beta)/rho)**2*cos(rho-2.0_dp*beta) & 
        & + 4.0_dp*(cos(beta)/rho)**2*cos(2.0_dp*beta)
      q_abs = 1.0_dp + 2.0_dp*(exp(-4.0_dp*k*x)/(4.0_dp*k*x)) &
        & + 2.0_dp*((exp(-4.0_dp*k*x)-1.0_dp)/(4.0_dp*k*x)**2)

      C1 = 0.25_dp * (1.0_dp + exp(-1167.0_dp*k))*(1.0_dp - q_abs)

      eps = 0.25_dp + 0.61_dp*(1.0_dp - exp(-(8.0_dp*pi/3.0_dp)*k))**2
      C2 = sqrt(2.0_dp*eps*(x/pi))*exp(0.5_dp - eps*(x/pi))*(0.79393_dp*n - 0.6069_dp)

      q_abs = (1.0_dp + C1 + C2)*q_abs 

      q_edge = (1.0_dp - exp(-0.06_dp*x))*x**(-2.0_dp/3.0_dp)

      q_ext = (1.0_dp + 0.5_dp*C2)*q_ext + q_edge

      q_sca = q_ext - q_abs

    end if

    !! Estimate g
    if (k < 1.0e-20_dp) then
      Cm = (6.0_dp + 5.0_dp*n**2 + n**4)/(45.0_dp*30.0_dp*n**2)
    else
      Cm = (-2.0_dp*k**6+k**4*(13.0_dp-2.0_dp*n**2)+k**2*(2.0_dp*n**4+2.0_dp*n**2-27.0_dp) & 
        & + 2.0_dp*n**6 + 13.0_dp*n**4 + 27.0_dp*n**2 + 18.0_dp) & 
        & /(15.0_dp * (4.0_dp*k**4 + 4.0_dp*k**2*(2.0_dp*n**2 - 3.0_dp) + (2.0_dp*n**2 + 3.0_dp)**2))
    end if

    !! Rayleigh regime, but limit to 0.9 to try get the constant region
    g = min(Cm * x**2, 0.9_dp)

  end subroutine madt

  subroutine read_nk_tables(n_dust, sp, n_wl, wl)
    implicit none

    integer, intent(in) :: n_dust, n_wl
    real(dp), dimension(n_wl), intent(in) :: wl
    character(len=11), dimension(n_dust), intent(in) :: sp

    integer :: n

    allocate(nk(n_dust))

    do n = 1, n_dust

      allocate(nk(n)%n(n_wl), nk(n)%k(n_wl))

      select case(trim(sp(n)))
      case('CaTiO3')
        nk(n)%name = sp(n)
        nk(n)%fname = 'CaTiO3[s].dat'
      case('TiO2')
        nk(n)%name = sp(n)
        nk(n)%fname = 'TiO2[s].dat'
      case('Al2O3')
        nk(n)%name = sp(n)
        nk(n)%fname = 'Al2O3[s].dat'
      case('Fe')
        nk(n)%name = sp(n)
        nk(n)%fname = 'Fe[s].dat'
      case('FeO')
        nk(n)%name = sp(n)
        nk(n)%fname = 'FeO[s].dat'        
      case('Mg2SiO4') 
        nk(n)%name = sp(n)
        nk(n)%fname = 'Mg2SiO4_amorph[s].dat'
      case('MgSiO3') 
        nk(n)%name = sp(n)
        nk(n)%fname = 'MgSiO3_amorph[s].dat'
      case('MnS') 
        nk(n)%name = sp(n)
        nk(n)%fname = 'MnS[s].dat'     
      case('Na2S') 
        nk(n)%name = sp(n)
        nk(n)%fname = 'Na2S[s].dat'
      case('ZnS') 
        nk(n)%name = sp(n)
        nk(n)%fname = 'ZnS[s].dat'
      case('KCl') 
        nk(n)%name = sp(n)
        nk(n)%fname = 'KCl[s].dat'   
      case('NaCl') 
        nk(n)%name = sp(n)
        nk(n)%fname = 'NaCl[s].dat'              
      case('C') 
        nk(n)%name = sp(n)
        nk(n)%fname = 'C[s].dat'
      case('H2O') 
        nk(n)%name = sp(n)
        nk(n)%fname = 'H2O[s].dat'
      case('NH3') 
        nk(n)%name = sp(n)
        nk(n)%fname = 'NH3[s].dat'
      case default
        print*,  'No availible n,k data for species: ', trim(sp(n)), 'STOP'
        stop
      end select

      call read_and_interp_nk_table(n, n_wl, wl)

    end do

  end subroutine read_nk_tables

  subroutine read_and_interp_nk_table(nn, n_wl, wl_in)
    implicit none

    integer, intent(in) :: nn, n_wl
    real(dp), dimension(n_wl), intent(in) :: wl_in

    integer :: u, l, nlines
    logical :: con_flag
    real(dp), allocatable, dimension(:) :: wl, n, k

    integer :: iwl1, iwl2
    real(dp) :: wl1, wl2, n1, n2, k1, k2

    open(newunit=u,file=trim(p_2_nk)//trim(nk(nn)%fname),action='read')
    print*, 'reading nk table @: ', trim(p_2_nk)//trim(nk(nn)%fname)

    ! Read number of lines and conducting flag
    read(u,*) nlines, con_flag

    allocate(wl(nlines),n(nlines),k(nlines))

    ! Read 4 blank lines
    read(u,*) ; read(u,*); read(u,*) ; read(u,*)

    do l = 1, nlines
      read(u,*) wl(l), n(l), k(l)
      n(l) = max(n(l),1e-30_dp)
      k(l) = max(k(l),1e-30_dp)
      !print*, l, wl(l), n(l), k(l)
    end do

    close(u)

    !! Perform 1D log-linear interpolation to get n,k at specific wavelengths
    do l = 1, n_wl

      ! Find wavelength array triplet and check bounds
      call locate(wl(:), nlines, wl_in(l), iwl1)

      iwl2 = iwl1 + 1

      if (iwl1 <= 0) then
        ! Use lowest wavelength n,k values in table
        nk(nn)%n(l) = n(1)
        nk(nn)%k(l) = k(1)
        cycle
      else if (iwl2 > nlines) then
        ! Use greatest wavelength n,k values in table
        nk(nn)%n(l) = n(nlines)
        nk(nn)%k(l) = k(nlines)
        cycle
      end if

      wl1 = wl(iwl1)
      wl2 = wl(iwl2)

      !! Interpolate n values
      n1 = n(iwl1)
      n2 = n(iwl2)
      call linear_log_interp(wl_in(l), wl1 , wl2, n1, n2, nk(nn)%n(l))

      !! Interpolate k values
      k1 = k(iwl1)
      k2 = k(iwl2)
      call linear_log_interp(wl_in(l), wl1 , wl2, k1, k2, nk(nn)%k(l))

    end do

    deallocate(wl, n, k)

  end subroutine read_and_interp_nk_table

  ! ------------- Functions for LLL theory ------------- !!
  complex(dp) function e2m(e)
    implicit none

    complex(dp), intent(in) :: e

    real(dp) :: ereal, eimag, n, k
    real(dp) :: sqrte2

    ereal = real(e,dp)
    eimag = aimag(e)
    sqrte2 = sqrt(ereal*ereal + eimag*eimag)
    n = sqrt(0.5_dp * ( ereal + sqrte2))
    k = sqrt(0.5_dp * (-ereal + sqrte2))
    e2m = cmplx(n,k,dp)

  end function e2m

  complex(dp) function m2e(m)
    implicit none

    complex(dp), intent(in) :: m

    real(dp) :: ereal, eimag, n, k

    n = real(m,dp)
    k = aimag(m)
    ereal = n*n - k*k
    eimag = 2.0_dp * n * k
    m2e = cmplx(ereal,eimag,dp)

  end function m2e  

  ! Perform linear interpolation in log10 space
  subroutine linear_log_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(dp), intent(in) :: xval, y1, y2, x1, x2

    real(dp), intent(out) :: yval

    real(dp) :: ly1, ly2
    real(dp) :: norm

    ly1 = log10(y1); ly2 = log10(y2)

    norm = 1.0_dp / log10(x2/x1)

    yval = 10.0_dp**((ly1 * log10(x2/xval) + ly2 * log10(xval/x1)) * norm)

  end subroutine linear_log_interp

  subroutine locate(arr, n, var, idx)
    implicit none

    integer, intent(in) :: n
    integer, intent(out) :: idx
    real(dp), dimension(n), intent(in) :: arr
    real(dp),intent(in) ::  var
    integer :: jl, jm, ju

    ! Search an array using bi-section (numerical methods)
    ! Then return array index that is lower than var in arr

    jl = 0
    ju = n+1
    do while (ju-jl > 1)
      jm = (ju+jl)/2
      if ((arr(n) > arr(1)).eqv.(var > arr(jm))) then
        jl=jm
      else
        ju=jm
      end if
    end do

    idx = jl

  end subroutine locate

end module mini_cloud_opac_mie_mod