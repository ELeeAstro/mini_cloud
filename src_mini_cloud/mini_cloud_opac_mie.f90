module mini_cloud_opac_mie
  use mini_cloud_precision
  use mini_cloud_class
  use lxmie_mod, only : lxmie
  implicit none

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
  private :: locate, bezier_interp, m2e, e2m, &
    & read_nk_tables, read_and_interp_nk_table

contains

  subroutine opac_mie(n_dust, sp, T_in, mu_in, P_in, k, k3, n_wl, wl, k_ext, alb, gg)
    implicit none

     integer, intent(in) :: n_dust, n_wl
    character(len=11), dimension(n_dust), intent(in) :: sp
    real(dp), intent(in) :: T_in, mu_in, P_in
    real(dp), dimension(4), intent(in) :: k
    real(dp), dimension(n_dust), intent(in) :: k3
    real(dp), dimension(n_wl), intent(in) :: wl

    real(dp), dimension(n_wl), intent(out) :: k_ext, alb, gg

    integer :: l, n
    real(dp) :: rho, T, mu, P_cgs
    real(dp) :: amean, x, xsec, q_ext, q_sca, q_abs, g
    real(dp), dimension(n_dust) :: b_mix
    complex(dp) :: e_eff, e_eff0
    complex(dp) :: N_eff, N_eff0
    complex(dp), dimension(n_dust) :: N_inc, e_inc


    if (first_call .eqv. .True.) then
      p_2_nk = 'nk_tables/'
      call read_nk_tables(n_dust, sp, n_wl, wl)
      first_call = .False.
    end if

    if (k(1) < 1e-6_dp) then
      k_ext(:) = 0.0_dp
      alb(:) = 0.0
      gg(:) = 0.0_dp
      return
    end if

    T = T_in
    mu = mu_in
    P_cgs = P_in * 10.0_dp
    rho = (P_cgs * mu * amu)/(kb * T)

    amean = k(2)/k(1)
    xsec = pi * amean**2
    b_mix(:) = max(k3(:)/k(4),1e-30_dp)
    b_mix(:) = min(k3(:)/k(4),1.0_dp)

    do l = 1, n_wl

      !! Size parameter 
      x = (2.0_dp * pi * (amean*1e4_dp)**2) / wl(l)

      if (x > 100.0_dp) then
        !! Use large size paramater approximation - guess alb and gg
        k_ext(l) = 2.0_dp * xsec * k(1)
        alb(l) = 0.5_dp
        gg(l) = 0.5_dp
        cycle 
      end if

      !! Refractive index
      !! Use LLL method
      do n = 1, n_dust
        N_inc(n) = cmplx(nk(n)%n(l),nk(n)%k(l),dp)
      end do

      N_eff = cmplx(0.0_dp, 0.0_dp,dp)
      e_eff = cmplx(0.0_dp, 0.0_dp,dp)

      N_eff0 = cmplx(0.0_dp, 0.0_dp,dp)
      do n = 1, n_dust
        N_eff0 = N_eff0 + b_mix(n) * N_inc(n)
        e_inc(n) = m2e(N_inc(n))
      end do

      e_eff0 = cmplx(0.0_dp,0.0_dp,dp)
      do n = 1, n_dust
        e_eff0 = e_eff0 + b_mix(n)*e_inc(n)**(1.0_dp/3.0_dp)
      end do
      e_eff = e_eff0**3
      N_eff = e2m(e_eff)

      !! Call LX-MIE with negative k value
      N_eff = cmplx(real(N_eff,dp),-aimag(N_eff),dp)
      call lxmie(N_eff, x, q_ext, q_sca, q_abs, g)

      k_ext(l) = (q_ext * xsec * k(1))/rho
      alb(l) = q_sca/q_ext
      gg(l) = g

    end do

  end subroutine opac_mie

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
      case('Mg2SiO4') 
        nk(n)%name = sp(n)
        nk(n)%fname = 'Mg2SiO4_amorph[s].dat'
      case('MgSiO3') 
        nk(n)%name = sp(n)
        nk(n)%fname = 'MgSiO3_amorph[s].dat'
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

    integer :: iwl1, iwl2, iwl3
    real(dp), dimension(3) :: lwl, ln, lk

    open(newunit=u,file=trim(p_2_nk)//trim(nk(nn)%fname),action='read')
    print*, 'reading nk table @: ', trim(p_2_nk)//trim(nk(nn)%fname)

    ! Read number of lines and conducting flag
    read(u,*) nlines, con_flag

    allocate(wl(nlines),n(nlines),k(nlines))

    ! Read 4 blank lines
    read(u,*) ; read(u,*); read(u,*) ; read(u,*)

    do l = 1, nlines
      read(u,*) wl(l), n(l), k(l)
    end do


    !! Perform 1D Bezier interpolation to get n,k at specific wavelengths
    do l = 1, n_wl

      ! Find wavelength array triplet and check bounds
      call locate(wl(:), nlines, wl_in(l), iwl2)

      iwl1 = iwl2 - 1
      iwl3 = iwl2 + 1

      if (iwl1 <= 0) then
        ! Use lowest wavelength n,k values in table
        nk(nn)%n(l) = n(1)
        nk(nn)%k(l) = k(1)
        cycle
      else if (iwl3 > nlines) then
        ! Use greatest wavelength n,k values in table
        nk(nn)%n(l) = n(nlines)
        nk(nn)%k(l) = k(nlines)
        cycle
      end if

      lwl(1) = wl(iwl1)
      lwl(2) = wl(iwl2)
      lwl(3) = wl(iwl3)

      !! Interpolate n and k values
      ln(1) = n(iwl1)
      ln(2) = n(iwl2)
      ln(3) = n(iwl3)
      call Bezier_interp(lwl(:), ln(:), 3, wl_in(l), nk(nn)%n(l))

      lk(1) = k(iwl1)
      lk(2) = k(iwl2)
      lk(3) = k(iwl3)
      call Bezier_interp(lwl(:), lk(:), 3, wl_in(l), nk(nn)%k(l))

    end do

    deallocate(wl, n, k)

    close(u)
  
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

  subroutine bezier_interp(xi, yi, ni, x, y)
    implicit none

    integer, intent(in) :: ni
    real(dp), dimension(ni), intent(in) :: xi, yi
    real(dp), intent(in) :: x
    real(dp), intent(out) :: y

    real(dp) :: dx, dx1, dy, dy1, w, yc, t, wlim, wlim1

    !xc = (xi(1) + xi(2))/2.0_dp ! Control point (no needed here, implicitly included)
    dx = xi(2) - xi(1)
    dx1 = xi(3) - xi(2)
    dy = yi(2) - yi(1)
    dy1 = yi(3) - yi(2)

    if (x > xi(1) .and. x < xi(2)) then
      ! left hand side interpolation
      !print*,'left'
      w = dx1/(dx + dx1)
      wlim = 1.0_dp + 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      wlim1 = 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      if (w <= min(wlim,wlim1) .or. w >= max(wlim,wlim1)) then
        w = 1.0_dp
      end if
      yc = yi(2) - dx/2.0_dp * (w*dy/dx + (1.0_dp - w)*dy1/dx1)
      t = (x - xi(1))/dx
      y = (1.0_dp - t)**2 * yi(1) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(2)
    else ! (x > xi(2) and x < xi(3)) then
      ! right hand side interpolation
      !print*,'right'
      w = dx/(dx + dx1)
      wlim = 1.0_dp/(1.0_dp - (dy1/dy) * (dx/dx1))
      wlim1 = 1.0_dp + 1.0_dp/(1.0_dp - (dy/dy1) * (dx1/dx))
      if (w <= min(wlim,wlim1) .or. w >= max(wlim,wlim1)) then
        w = 1.0_dp
      end if
      yc = yi(2) + dx1/2.0_dp * (w*dy1/dx1 + (1.0_dp - w)*dy/dx)
      t = (x - xi(2))/(dx1)
      y = (1.0_dp - t)**2 * yi(2) + 2.0_dp*t*(1.0_dp - t)*yc + t**2*yi(3)
    end if

  end subroutine bezier_interp

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

end module mini_cloud_opac_mie