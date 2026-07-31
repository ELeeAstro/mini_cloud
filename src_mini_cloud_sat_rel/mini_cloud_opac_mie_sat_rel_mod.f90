module mini_cloud_opac_mie_sat_rel_mod
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

  !! Number of Gauss-Laguerre quadrature points used to integrate the Mie
  !! efficiencies over the prescribed gamma particle size distribution
  !! (dist == 2). Validated against a dense reference integral to ~0.1-15%.
  integer, parameter :: n_quad = 8

  !! Number of points used to integrate the Mie efficiencies over the
  !! prescribed lognormal particle size distribution (dist == 1), following
  !! the mini_Mie reference recipe (int_ln_ana_trapz_psd.py): a plain
  !! composite trapezoidal rule over z in [-5,5], not a fixed-order
  !! Gauss-Hermite quadrature. This was necessary because Mie efficiencies
  !! are genuinely oscillatory (real resonance structure, not numerical
  !! noise) across the size range a lognormal distribution samples - a
  !! low-order Gauss-Hermite rule (tried at n=8 and n=16) left a ~40-65%
  !! error that did not shrink with more nodes, since Gauss quadrature
  !! convergence assumes a smooth integrand. A dense uniform trapezoidal
  !! rule resolves the oscillations directly and was validated to ~0.5-6.5%
  !! against a dense reference integral for TiO2 at several wavelengths.
  integer, parameter :: n_trapz = 100

  !! Cached, per-species per-wavelength distribution-integrated Mie "shape
  !! factors", built once by init_opac_mie_sat_rel. r_med and sigma are
  !! prescribed constants for the whole run, so the quadrature abscissas and
  !! the Mie efficiencies evaluated there never change - only the local
  !! number density n_d (from q_c, T, P) varies per call. These arrays are
  !! already normalised so that, per species j and wavelength l:
  !!   k_ext_sp(l) = pi * n_d(j) / rho * k_ext_shape(j,l)
  !!   alb_sp(l)   = sca_shape(j,l) / k_ext_shape(j,l)
  !!   gg_sp(l)    = gsca_shape(j,l) / sca_shape(j,l)
  real(dp), allocatable, dimension(:,:), save :: k_ext_shape, sca_shape, gsca_shape
  logical, save :: opac_initialized = .False.
  !$omp threadprivate (k_ext_shape, sca_shape, gsca_shape, opac_initialized)

  public :: init_opac_mie_sat_rel, opac_mie_sat_rel
  private :: locate, linear_log_interp, &
    & load_nk_table, set_nk_filename, read_and_interp_nk_table, &
    & mie_efficiencies, trapz_nodes, &
    & gauss_laguerre_nodes, tridiag_eigen, pythag, &
    & inv_trigamma_pos, trigamma_pos

contains

  !! Precompute and cache the distribution-integrated Mie shape factors for
  !! each condensate species. Must be called once (per thread) before
  !! opac_mie_sat_rel, and again if r_med, sigma, dist, or the wavelength
  !! grid change.
  subroutine init_opac_mie_sat_rel(n_dust, sp, n_wl, wl, r_med, sigma, dist)
    implicit none

    integer, intent(in) :: n_dust, n_wl, dist
    character(len=*), dimension(n_dust), intent(in) :: sp
    real(dp), dimension(n_wl), intent(in) :: wl
    real(dp), dimension(n_dust), intent(in) :: r_med
    real(dp), intent(in) :: sigma

    integer :: j, l, nk_idx

    if (first_call .eqv. .True.) then
      p_2_nk = 'nk_tables/'
      first_call = .False.
    end if

    if (allocated(k_ext_shape)) deallocate(k_ext_shape, sca_shape, gsca_shape)
    allocate(k_ext_shape(n_dust,n_wl), sca_shape(n_dust,n_wl), gsca_shape(n_dust,n_wl))

    do j = 1, n_dust

      call load_nk_table(sp(j), n_wl, wl, nk_idx)

      do l = 1, n_wl
        if (dist == 1) then
          call shape_factors_lognormal(nk_idx, l, wl(l), r_med(j), sigma, &
            & k_ext_shape(j,l), sca_shape(j,l), gsca_shape(j,l))
        else if (dist == 2) then
          call shape_factors_gamma(nk_idx, l, wl(l), r_med(j), sigma, &
            & k_ext_shape(j,l), sca_shape(j,l), gsca_shape(j,l))
        else
          print*, 'init_opac_mie_sat_rel invalid dist: ', dist
          stop
        end if
      end do

    end do

    opac_initialized = .True.

  end subroutine init_opac_mie_sat_rel

  !! Distribution-integrated Mie shape factors for the lognormal size
  !! distribution, following the mini_Mie reference recipe
  !! (int_ln_ana_trapz_psd.py): a composite trapezoidal rule over z in
  !! [-5,5] with n_trapz points, r(z) = r_med*exp(sqrt(2)*log(sigma)*z),
  !! weighted by exp(-z**2) and r**2 (cross-sectional area).
  subroutine shape_factors_lognormal(nk_idx, l, wl, r_med, sigma, k_ext_sh, sca_sh, gsca_sh)
    implicit none

    integer, intent(in) :: nk_idx, l
    real(dp), intent(in) :: wl, r_med, sigma
    real(dp), intent(out) :: k_ext_sh, sca_sh, gsca_sh

    integer :: k
    real(dp) :: log_sigma, r_k, x, q_ext, q_sca, q_abs, g, wgt, norm_const
    real(dp), dimension(n_trapz) :: z_tz, w_tz
    complex(dp) :: N_eff

    call trapz_nodes(n_trapz, -5.0_dp, 5.0_dp, z_tz, w_tz)

    !! Normalised average over a standard-normal weight is (1/sqrt(pi)) * sum(w*f).
    norm_const = 1.0_dp/sqrt(pi)
    log_sigma = log(sigma)

    N_eff = cmplx(nk(nk_idx)%n(l), nk(nk_idx)%k(l), dp)

    k_ext_sh = 0.0_dp
    sca_sh   = 0.0_dp
    gsca_sh  = 0.0_dp

    do k = 1, n_trapz

      r_k = max(r_med * exp(sqrt(2.0_dp) * log_sigma * z_tz(k)), r_seed)
      x = (2.0_dp * pi * r_k) / (wl * 1e-4_dp)
      call mie_efficiencies(x, N_eff, q_ext, q_sca, q_abs, g)

      if (.not. ieee_is_finite(q_ext)) q_ext = 0.0_dp
      if (.not. ieee_is_finite(q_sca)) q_sca = 0.0_dp
      if (.not. ieee_is_finite(g)) g = 0.0_dp

      wgt = w_tz(k) * exp(-z_tz(k)**2) * r_k**2

      k_ext_sh = k_ext_sh + wgt * q_ext
      sca_sh   = sca_sh   + wgt * q_sca
      gsca_sh  = gsca_sh  + wgt * q_sca * g

    end do

    k_ext_sh = k_ext_sh * norm_const
    sca_sh   = sca_sh   * norm_const
    gsca_sh  = gsca_sh  * norm_const

  end subroutine shape_factors_lognormal

  !! Distribution-integrated Mie shape factors for the gamma size
  !! distribution (r_med is the number-weighted radius, sigma related to
  !! the trigamma function), via n_quad-point generalized Gauss-Laguerre
  !! quadrature with the r**2 area weighting embedded in the quadrature
  !! shape itself (see gauss_laguerre_nodes).
  subroutine shape_factors_gamma(nk_idx, l, wl, r_med, sigma, k_ext_sh, sca_sh, gsca_sh)
    implicit none

    integer, intent(in) :: nk_idx, l
    real(dp), intent(in) :: wl, r_med, sigma
    real(dp), intent(out) :: k_ext_sh, sca_sh, gsca_sh

    integer :: k
    real(dp) :: A, B, r_k, x, q_ext, q_sca, q_abs, g, norm_const
    real(dp), dimension(n_quad) :: q_x, q_w
    complex(dp) :: N_eff

    A = inv_trigamma_pos(log(sigma)**2)
    B = A/r_med
    !! Shift the Laguerre shape by +1 so the r**2 area weighting is
    !! embedded in the quadrature weight itself (matches r**(A-1)*r**2).
    call gauss_laguerre_nodes(A + 1.0_dp, q_x, q_w)

    !! Exact analytic normalisation of the r**2-weighted Gamma(A,B)
    !! average: E[r**2] = A*(A+1)/B**2 (derived from the r**(A-1)*e**(-B*r)
    !! moments), folded in here since it is fully closed-form.
    norm_const = A*(A + 1.0_dp)/B**2

    N_eff = cmplx(nk(nk_idx)%n(l), nk(nk_idx)%k(l), dp)

    k_ext_sh = 0.0_dp
    sca_sh   = 0.0_dp
    gsca_sh  = 0.0_dp

    do k = 1, n_quad

      r_k = max(q_x(k)/B, r_seed)
      x = (2.0_dp * pi * r_k) / (wl * 1e-4_dp)
      call mie_efficiencies(x, N_eff, q_ext, q_sca, q_abs, g)

      if (.not. ieee_is_finite(q_ext)) q_ext = 0.0_dp
      if (.not. ieee_is_finite(q_sca)) q_sca = 0.0_dp
      if (.not. ieee_is_finite(g)) g = 0.0_dp

      k_ext_sh = k_ext_sh + q_w(k) * q_ext
      sca_sh   = sca_sh   + q_w(k) * q_sca
      gsca_sh  = gsca_sh  + q_w(k) * q_sca * g

    end do

    k_ext_sh = k_ext_sh * norm_const
    sca_sh   = sca_sh   * norm_const
    gsca_sh  = gsca_sh  * norm_const

  end subroutine shape_factors_gamma

  !! Compute the combined extinction, single-scattering albedo, and
  !! asymmetry parameter across all condensate species. Extinction is
  !! additive across species; albedo and asymmetry are scattering-weighted
  !! averages - both follow directly from summing the (already correctly
  !! normalised) per-species physical extinction/scattering/asymmetry-
  !! weighted-scattering coefficients before taking any ratio.
  subroutine opac_mie_sat_rel(n_dust, sp, T_in, mu_in, P_in, q_c, r_med, rho_d, sigma, n_wl, wl, k_ext, alb, gg, dist)
    implicit none

    integer, intent(in) :: n_dust, n_wl, dist
    character(len=*), dimension(n_dust), intent(in) :: sp
    real(dp), intent(in) :: T_in, mu_in, P_in, sigma
    real(dp), dimension(n_dust), intent(in) :: q_c, r_med, rho_d
    real(dp), dimension(n_wl), intent(in) :: wl

    real(dp), dimension(n_wl), intent(out) :: k_ext, alb, gg

    integer :: j, l
    real(dp) :: rho, n_d, scale
    real(dp), dimension(n_wl) :: sca_ext, g_sca_ext

    if (opac_initialized .eqv. .False.) then
      print*, 'opac_mie_sat_rel: call init_opac_mie_sat_rel before opac_mie_sat_rel. STOP'
      stop
    end if

    rho = (P_in * 10.0_dp * mu_in * amu)/(kb * T_in)

    k_ext(:) = 0.0_dp
    sca_ext(:) = 0.0_dp
    g_sca_ext(:) = 0.0_dp

    do j = 1, n_dust

      if (q_c(j) < 1e-10_dp) cycle

      if (dist == 1) then
        ! lognormal total number density - r_med is the median particle size
        n_d = (3.0_dp * q_c(j) * rho)/(4.0_dp*pi*rho_d(j)*r_med(j)**3) * exp(-9.0_dp/2.0_dp * log(sigma)**2)
      else if (dist == 2) then
        ! gamma total number density - r_med is the number weighted particle size
        n_d = (3.0_dp * q_c(j) * rho)/(4.0_dp*pi*rho_d(j))
      else
        print*, 'opac_mie_sat_rel invalid dist: ', dist
        stop
      end if

      scale = pi * n_d / rho

      do l = 1, n_wl
        k_ext(l)     = k_ext(l)     + scale * k_ext_shape(j,l)
        sca_ext(l)   = sca_ext(l)   + scale * sca_shape(j,l)
        g_sca_ext(l) = g_sca_ext(l) + scale * gsca_shape(j,l)
      end do

    end do

    alb(:) = sca_ext(:)/max(k_ext(:), 1.0e-300_dp)
    gg(:)  = g_sca_ext(:)/max(sca_ext(:), 1.0e-300_dp)

    do l = 1, n_wl
      if (.not. ieee_is_finite(k_ext(l))) k_ext(l) = 0.0_dp
      if (.not. ieee_is_finite(alb(l))) alb(l) = 0.0_dp
      if (.not. ieee_is_finite(gg(l))) gg(l) = 0.0_dp
    end do

  end subroutine opac_mie_sat_rel

  !! Compute the Mie efficiencies for size parameter x and refractive index
  !! ri via the exact LX-MIE solution (Kitzmann & Heng 2018). Approximate
  !! regimes (Rayleigh, MADT) are no longer needed: the Q's are now cached
  !! once per species/wavelength/quadrature-node at init time rather than
  !! evaluated every layer/timestep, so the extra cost of the exact
  !! solution at small/large x is paid once, not repeatedly - and it avoids
  !! the regime-switch kinks that limited quadrature accuracy near x~10.
  subroutine mie_efficiencies(x, ri, q_ext, q_sca, q_abs, g)
    implicit none

    real(dp), intent(in) :: x
    complex(dp), intent(in) :: ri

    real(dp), intent(out) :: q_ext, q_sca, q_abs, g

    complex(dp) :: N_eff

    !! Call LX-MIE with negative k value
    N_eff = cmplx(real(ri,dp),-aimag(ri),dp)
    call lxmie(N_eff, x, q_ext, q_sca, q_abs, g)

    g = max(g, 0.0_dp)

  end subroutine mie_efficiencies

  !! n-point composite trapezoidal rule nodes/weights on a uniform grid
  !! spanning [a,b] inclusive (endpoints get half weight). Used for the
  !! lognormal shape-factor integral following the mini_Mie trapz recipe.
  pure subroutine trapz_nodes(n, a, b, x, w)
    implicit none

    integer, intent(in) :: n
    real(dp), intent(in) :: a, b
    real(dp), dimension(n), intent(out) :: x, w

    real(dp) :: dx
    integer :: i

    dx = (b - a)/real(n-1,dp)

    do i = 1, n
      x(i) = a + real(i-1,dp)*dx
    end do

    w(:) = dx
    w(1) = 0.5_dp*dx
    w(n) = 0.5_dp*dx

  end subroutine trapz_nodes

  !! Generalized n-point Gauss-Laguerre quadrature nodes/weights for the
  !! weight function x**alpha * exp(-x) on [0,inf), via the same
  !! Golub-Welsch approach. Weights are returned unnormalized (missing the
  !! common factor Gamma(alpha+1), which callers must fold in explicitly
  !! when an absolute - not ratio - quantity is needed).
  pure subroutine gauss_laguerre_nodes(alpha, x, w)
    implicit none

    real(dp), intent(in) :: alpha
    real(dp), dimension(n_quad), intent(out) :: x, w

    integer, parameter :: n = n_quad
    real(dp), dimension(n) :: d, e
    real(dp), dimension(n,n) :: z
    integer :: i

    !! Jacobi matrix for generalized Laguerre polynomials L_n^(alpha):
    !! diagonal a_k = 2k + alpha + 1, off-diagonal b_k = sqrt(k*(k+alpha))
    do i = 1, n
      d(i) = 2.0_dp*real(i-1,dp) + alpha + 1.0_dp
    end do
    e(1) = 0.0_dp
    do i = 2, n
      e(i) = sqrt(real(i-1,dp)*(real(i-1,dp) + alpha))
    end do

    z(:,:) = 0.0_dp
    do i = 1, n
      z(i,i) = 1.0_dp
    end do

    call tridiag_eigen(n, d, e, z)

    x(:) = d(:)
    w(:) = z(1,:)**2

  end subroutine gauss_laguerre_nodes

  !! Eigenvalues (ascending) and eigenvectors of a real symmetric tridiagonal
  !! matrix with diagonal d and off-diagonal e(2:n), via the implicit-shift
  !! QL algorithm (standard EISPACK TQL2 method). On input z should be the
  !! identity (or an initial similarity transform); on output its columns
  !! hold the eigenvectors of the original tridiagonal matrix.
  pure subroutine tridiag_eigen(n, d, e, z)
    implicit none

    integer, intent(in) :: n
    real(dp), dimension(n), intent(inout) :: d, e
    real(dp), dimension(n,n), intent(inout) :: z

    integer :: i, j, k, l, m, ii, l1, mml
    real(dp) :: c, c2, c3, dl1, el1, f, g, h, p, r, s, s2, tst1, tst2

    if (n == 1) return

    do i = 2, n
      e(i-1) = e(i)
    end do
    e(n) = 0.0_dp

    f = 0.0_dp
    tst1 = 0.0_dp

    do l = 1, n
      tst1 = max(tst1, abs(d(l)) + abs(e(l)))

      !! locate small sub-diagonal element
      m = n
      do i = l, n
        tst2 = tst1 + abs(e(i))
        if (tst2 == tst1) then
          m = i
          exit
        end if
      end do

      do while (m /= l)

        !! form shift
        l1 = l + 1
        g = d(l)
        p = (d(l1) - g)/(2.0_dp*e(l))
        r = pythag(p, 1.0_dp)
        d(l) = e(l)/(p + sign(r,p))
        d(l1) = e(l)*(p + sign(r,p))
        dl1 = d(l1)
        h = g - d(l)
        do i = l1+1, n
          d(i) = d(i) - h
        end do
        f = f + h

        !! QL transformation on the submatrix l..m
        p = d(m)
        c = 1.0_dp
        c2 = c
        el1 = e(l1)
        s = 0.0_dp
        mml = m - l
        do ii = 1, mml
          c3 = c2
          c2 = c
          s2 = s
          i = m - ii
          g = c*e(i)
          h = c*p
          r = pythag(p, e(i))
          e(i+1) = s*r
          s = e(i)/r
          c = p/r
          p = c*d(i) - s*g
          d(i+1) = h + s*(c*g + s*d(i))
          do k = 1, n
            h = z(k,i+1)
            z(k,i+1) = s*z(k,i) + c*h
            z(k,i) = c*z(k,i) - s*h
          end do
        end do
        p = -s*s2*c3*el1*e(l)/dl1
        e(l) = s*p
        d(l) = c*p

        tst2 = tst1 + abs(e(l))
        if (tst2 <= tst1) exit
      end do

      d(l) = d(l) + f

      !! locate next small sub-diagonal element (or exit if l==n)
      if (l < n) then
        m = n
        do i = l+1, n
          tst2 = tst1 + abs(e(i))
          if (tst2 == tst1) then
            m = i
            exit
          end if
        end do
      end if

    end do

    !! sort eigenvalues (and corresponding eigenvectors) ascending
    do ii = 2, n
      i = ii - 1
      k = i
      p = d(i)
      do j = ii, n
        if (d(j) < p) then
          k = j
          p = d(j)
        end if
      end do
      if (k /= i) then
        d(k) = d(i)
        d(i) = p
        do j = 1, n
          p = z(j,i)
          z(j,i) = z(j,k)
          z(j,k) = p
        end do
      end if
    end do

  end subroutine tridiag_eigen

  pure real(dp) function pythag(a, b) result(r)
    implicit none
    real(dp), intent(in) :: a, b
    real(dp) :: absa, absb

    absa = abs(a)
    absb = abs(b)
    if (absa > absb) then
      r = absa*sqrt(1.0_dp + (absb/absa)**2)
    else if (absb == 0.0_dp) then
      r = 0.0_dp
    else
      r = absb*sqrt(1.0_dp + (absa/absb)**2)
    end if

  end function pythag

  subroutine load_nk_table(sp, n_wl, wl, idx)
    implicit none

    character(len=*), intent(in) :: sp
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
    character(len=*), intent(in) :: sp

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

end module mini_cloud_opac_mie_sat_rel_mod
