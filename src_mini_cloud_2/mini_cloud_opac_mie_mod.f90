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

  subroutine opac_mie(n_dust, sp, T_in, mu_in, P_in, q_0, q_1, rho_d, n_wl, wl, k_ext, alb, gg)
    implicit none

    integer, intent(in) :: n_dust, n_wl
    character(len=11), dimension(n_dust), intent(in) :: sp
    real(dp), intent(in) :: T_in, mu_in, P_in
    real(dp), intent(in) :: q_0, q_1, rho_d
    real(dp), dimension(n_wl), intent(in) :: wl

    real(dp), dimension(n_wl), intent(out) :: k_ext, alb, gg

    integer :: l, n
    real(dp) :: rho, nd_atm, amean, n_d, m_c
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

    rho = (P_in * 10.0_dp * mu_in * amu)/(kb * T_in)
    nd_atm = (P_in * 10.0_dp)/(kb * T_in)

    if (q_0*nd_atm < 1e-10_dp) then
      k_ext(:) = 0.0_dp
      alb(:) = 0.0_dp
      gg(:) = 0.0_dp
      return
    end if

    m_c = (q_1*rho)/(q_0*nd_atm)

    amean =  max(((3.0_dp*m_c)/(4.0_dp*pi*rho_d))**(1.0_dp/3.0_dp), 1e-7_dp)

    n_d = q_0 * nd_atm

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
        call rayleigh(x, N_eff, q_abs, q_sca, q_ext, g)
      else if (x > 100.0_dp) then
        call adt(x, N_eff, q_abs, q_sca, q_ext, g)
      else
        !! Call LX-MIE with negative k value
        N_eff = cmplx(real(N_eff,dp),-aimag(N_eff),dp)
        call lxmie(N_eff, x, q_ext, q_sca, q_abs, g)
      end if

      !! Calculate the opacity, abledo and mean cosine angle (asymmetry factor)
      k_ext(l) = (q_ext * xsec * n_d)/rho
      alb(l) = min(q_sca/q_ext, 0.95_dp)
      gg(l) = max(g, 0.0_dp)

    end do

  end subroutine opac_mie

  ! Small dielectric sphere approximation, small particle limit x << 1
  subroutine rayleigh(x, ri, q_abs, q_sca, q_ext, g)
    implicit none

    real(dp), intent(in) :: x
    complex(dp), intent(in) :: ri

    real(dp), intent(out) :: q_abs, q_sca, q_ext, g

    complex(dp) :: alp

    alp = (ri**2 - 1.0_dp)/(ri**2 + 2.0_dp)

    q_sca = 8.0_dp/3.0_dp * x**4 * abs(alp)**2
    q_ext = 4.0_dp * x * aimag(alp * (1.0_dp + x**2/15.0_dp*alp * ((ri**4+27.0_dp*ri**2+38.0_dp)/(2.0_dp*ri**2+3.0_dp)))) + &
     & 8.0_dp/3.0_dp * x**4 * real(alp**2,dp)

    q_abs = q_ext - q_sca

    g = 0.0_dp

  end subroutine rayleigh

  subroutine adt(x, ri, q_abs, q_sca, q_ext, g)
    implicit none

    real(dp), intent(in) :: x
    complex(dp), intent(in) :: ri

    real(dp), intent(out) :: q_abs, q_sca, q_ext, g

    real(dp) :: rho1, rho2, rho0, beta0, beta, fac, fac2
    complex(dp) ::  rho

    rho = 2.0_dp * x * (ri - 1.0_dp)
    rho0 = abs(rho)
    rho1 = real(rho, dp)
    rho2 = aimag(rho)

    if (abs(rho1) > 0.0_dp) then
      beta0 = atan(abs(rho2)/abs(rho1))
      if (rho1 < 0.0_dp .and. rho2 > 0.0_dp) then
        beta = pi - beta0
      else if (rho1 < 0.0_dp .and. rho2 < 0.0_dp) then
        beta = pi + beta0
      else if (rho1 > 0.0_dp .and. rho2 < 0.0_dp) then 
        beta = 2.0_dp*pi - beta0
      else 
        beta = beta0 
      endif     
    else
      if (rho2 > 0.0_dp) then
        beta = 0.5_dp*pi 
      else
        beta = 1.5_dp*pi
      endif
    end if

    if (rho0 < 1.0e-3_dp) then
      q_ext = (4.0_dp/3.0_dp)*rho2 + 0.5_dp*(rho1**2 - rho2**2)
      q_abs = (4.0_dp/3.0_dp)*rho2 - rho2**2
      q_sca = 0.5_dp*rho0**2
    else
      fac = exp(-rho2)
      fac2 = fac**2
      q_ext = 2.0_dp + (4.0_dp/rho0**2)*(cos(2.0_dp*beta) - fac*(cos(rho1 - 2*beta) + rho0*sin(rho1 - beta)))
      q_abs = 1.0_dp + fac2/rho2 + (fac2 - 1.0_dp)/(2.0_dp*rho2**2)
      q_sca = q_ext - q_abs
    end if

    if (x >= 100.0_dp) then
      q_ext = q_ext + 2.0_dp * x**(-2.0_dp/3.0_dp)
    end if

    if (aimag(ri) < 1) then

      if (x <= 3.0_dp) then
        g = 0.7_dp*(x/3.0_dp)**2
      else
        g = 0.7_dp
      end if

    else if (aimag(ri) >= 1) then

      if (x <= 3.0_dp) then
        g = -0.2_dp
      else
        g = 0.5_dp
      end if

    end if

  end subroutine adt

  subroutine read_nk_tables(n_dust, sp, n_wl, wl)
    implicit none

    integer, intent(in) :: n_dust, n_wl
    real(dp), dimension(n_wl), intent(in) :: wl
    character(len=11), dimension(n_dust), intent(in) :: sp

    integer :: n

    allocate(nk(n_dust))

    do n = 1, n_dust

      allocate(nk(n)%n(n_wl))
      allocate(nk(n)%k(n_wl))

      select case(trim(sp(n)))
      case('CaTiO3')
        nk(n)%name = sp(n)
        nk(n)%fname = 'CaTiO3[s].dat'
      case('TiO2')
        nk(n)%name = sp(n)
        nk(n)%fname = 'TiO2_anatase[s].dat'
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
      case('KCl') 
        nk(n)%name = sp(n)
        nk(n)%fname = 'KCl[s].dat'        
      case('C') 
        nk(n)%name = sp(n)
        nk(n)%fname = 'C[s].dat'
      case('H2O') 
        nk(n)%name = sp(n)
        nk(n)%fname = 'H2O[s].dat'
       case('NH3') 
        nk(n)%name = sp(n)
        nk(n)%fname = 'NH3[s].dat'
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