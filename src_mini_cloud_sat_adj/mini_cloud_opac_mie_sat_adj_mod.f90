module mini_cloud_opac_mie_sat_adj_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
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
  !$omp threadprivate (p_2_nk)

  type(nk_table), allocatable, dimension(:)  :: nk
  !$omp threadprivate (nk)

  logical :: first_call = .True.
  !$omp threadprivate (first_call)

  real(dp), parameter :: r_seed = 1e-7_dp

  public :: opac_mie_sat_adj
  private :: locate, linear_log_interp, &
    & load_nk_table, set_nk_filename, read_and_interp_nk_table

contains

  subroutine opac_mie_sat_adj(n_dust, sp, T_in, mu_in, P_in, q_c, r_c, rho_d, sigma, n_wl, wl, k_ext, alb, gg, dist)
    implicit none

    integer, intent(in) :: n_dust, n_wl, dist
    character(len=11), dimension(n_dust), intent(in) :: sp
    real(dp), intent(in) :: T_in, mu_in, P_in
    real(dp), intent(in) :: q_c, r_c, sigma, rho_d
    real(dp), dimension(n_wl), intent(in) :: wl

    real(dp), dimension(n_wl), intent(out) :: k_ext, alb, gg

    integer :: l, nk_idx
    real(dp) :: rho, r_eff, n_d
    real(dp) :: x, xsec, q_ext, q_sca, q_abs, g
    complex(dp) :: N_eff

    real(dp) :: A, B

    if (n_dust /= 1) then
      print*, 'opac_mie_sat_adj expects one condensate species when effective medium theory is disabled.'
      stop
    end if

    if (first_call .eqv. .True.) then
      p_2_nk = 'nk_tables/'
      first_call = .False.
    end if
    call load_nk_table(sp(1), n_wl, wl, nk_idx)

    if (q_c < 1e-10_dp) then
      k_ext(:) = 0.0_dp
      alb(:) = 0.0
      gg(:) = 0.0_dp
      return
    end if

    rho = (P_in * 10.0_dp * mu_in * amu)/(kb * T_in)


    if (dist == 1) then
      ! lognormal total number density - r_c is the median particle size
      n_d = (3.0_dp * q_c * rho)/(4.0_dp*pi*rho_d*r_c**3) * exp(-9.0_dp/2.0_dp * log(sigma)**2)
      ! Find effective particle size [cm]
      r_eff = max(r_c * exp(5.0_dp/2.0_dp * log(sigma)**2), r_seed)
    else if (dist == 2) then
      ! gamma total number density - r_c is the number weighted particle size
      A = inv_trigamma_pos(log(sigma)**2)
      B = A/r_c
      n_d = (3.0_dp * q_c * rho)/(4.0_dp*pi*rho_d)
      ! Find effective particle size [cm]
      r_eff = max((A + 2.0_dp)/B, r_seed)
    else
      print*, 'opac_mie_sat_adj invalid dist: ', dist
      stop
    end if
    

    xsec = pi * r_eff**2
    
    do l = 1, n_wl

      !! Size parameter 
      x = (2.0_dp * pi * r_eff) / (wl(l) * 1e-4_dp)

      !! Refractive index for the single condensate species.
      N_eff = cmplx(nk(nk_idx)%n(l), nk(nk_idx)%k(l), dp)

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
      alb(l) = min(q_sca/q_ext, 0.99_dp)
      gg(l) = max(g, 0.0_dp)

      if (.not. ieee_is_finite(k_ext(l))) k_ext(l) = 0.0_dp
      if (.not. ieee_is_finite(alb(l))) alb(l) = 0.0_dp
      if (.not. ieee_is_finite(gg(l))) gg(l) = 0.0_dp

    end do

  end subroutine opac_mie_sat_adj

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
      Cm = (6.0_dp + 5.0_dp*n**2 + n**4)/(45.0_dp + 30.0_dp*n**2)
    else
      Cm = (-2.0_dp*k**6+k**4*(13.0_dp-2.0_dp*n**2)+k**2*(2.0_dp*n**4+2.0_dp*n**2-27.0_dp) & 
        & + 2.0_dp*n**6 + 13.0_dp*n**4 + 27.0_dp*n**2 + 18.0_dp) & 
        & /(15.0_dp * (4.0_dp*k**4 + 4.0_dp*k**2*(2.0_dp*n**2 - 3.0_dp) + (2.0_dp*n**2 + 3.0_dp)**2))
    end if

    !! Rayleigh regime, but limit to 0.9 to try get the constant region
    g = min(Cm * x**2, 0.9_dp)

  end subroutine madt

  subroutine load_nk_table(sp, n_wl, wl, idx)
    implicit none

    character(len=11), intent(in) :: sp
    integer, intent(in) :: n_wl
    real(dp), dimension(n_wl), intent(in) :: wl
    integer, intent(out) :: idx

    integer :: n, n_old
    type(nk_table), allocatable, dimension(:) :: nk_old

    idx = 0
    if (allocated(nk)) then
      do n = 1, size(nk)
        if (trim(nk(n)%name) == trim(sp)) then
          idx = n
          return
        end if
      end do

      n_old = size(nk)
      call move_alloc(nk, nk_old)
      allocate(nk(n_old+1))

      do n = 1, n_old
        nk(n)%name = nk_old(n)%name
        nk(n)%fname = nk_old(n)%fname
        call move_alloc(nk_old(n)%n, nk(n)%n)
        call move_alloc(nk_old(n)%k, nk(n)%k)
      end do
    else
      n_old = 0
      allocate(nk(1))
    end if

    idx = n_old + 1
    allocate(nk(idx)%n(n_wl), nk(idx)%k(n_wl))
    call set_nk_filename(idx, sp)
    call read_and_interp_nk_table(idx, n_wl, wl)

  end subroutine load_nk_table

  subroutine set_nk_filename(idx, sp)
    implicit none

    integer, intent(in) :: idx
    character(len=11), intent(in) :: sp

    select case(trim(sp))
    case('CaTiO3')
      nk(idx)%name = sp
      nk(idx)%fname = 'CaTiO3[s].dat'
    case('TiO2')
      nk(idx)%name = sp
      nk(idx)%fname = 'TiO2[s].dat'
    case('Al2O3')
      nk(idx)%name = sp
      nk(idx)%fname = 'Al2O3[s].dat'
    case('Fe')
      nk(idx)%name = sp
      nk(idx)%fname = 'Fe[s].dat'
    case('FeO')
      nk(idx)%name = sp
      nk(idx)%fname = 'FeO[s].dat'
    case('Mg2SiO4')
      nk(idx)%name = sp
      nk(idx)%fname = 'Mg2SiO4_amorph[s].dat'
    case('MgSiO3')
      nk(idx)%name = sp
      nk(idx)%fname = 'MgSiO3_amorph[s].dat'
    case('MnS')
      nk(idx)%name = sp
      nk(idx)%fname = 'MnS[s].dat'
    case('Na2S')
      nk(idx)%name = sp
      nk(idx)%fname = 'Na2S[s].dat'
    case('ZnS')
      nk(idx)%name = sp
      nk(idx)%fname = 'ZnS[s].dat'
    case('KCl')
      nk(idx)%name = sp
      nk(idx)%fname = 'KCl[s].dat'
    case('NaCl')
      nk(idx)%name = sp
      nk(idx)%fname = 'NaCl[s].dat'
    case('C')
      nk(idx)%name = sp
      nk(idx)%fname = 'C[s].dat'
    case('H2O')
      nk(idx)%name = sp
      nk(idx)%fname = 'H2O[s].dat'
    case('NH3')
      nk(idx)%name = sp
      nk(idx)%fname = 'NH3[s].dat'
    case('Soot_Lavvas')
      nk(idx)%name = sp
      nk(idx)%fname = 'Soot_Lavvas[s].dat'
    case default
      print*,  'No availible n,k data for species: ', trim(sp), 'STOP'
      stop
    end select

  end subroutine set_nk_filename

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

end module mini_cloud_opac_mie_sat_adj_mod
