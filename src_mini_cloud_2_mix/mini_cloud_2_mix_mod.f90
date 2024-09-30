module mini_cloud_2_mix_mod
  use  mini_cloud_class_mod
  implicit none

  !! Global variables
  real(dp) :: T, mu, nd_atm, rho, p, grav, Rd

  real(dp) :: mfp, eta, nu, eps_d

  logical :: first_call = .True.

  integer :: n_dust

  public :: mini_cloud_2_mix, RHS_mom, jac_dum
  private :: calc_coal, calc_coag, calc_cond, calc_hom_nuc, calc_seed_evap, &
    calc_turb_acc, calc_turb_shear, calc_turb_couple

  contains

  subroutine mini_cloud_2_mix(T_in, P_in, grav_in, mu_in, VMR_bg_in, t_end, sp, sp_bg, q_v, q_0, q_1s)
    implicit none

    ! Input variables
    character(len=20), dimension(:), intent(in) :: sp
    character(len=20), dimension(:), intent(in) :: sp_bg
    real(dp), intent(in) :: T_in, P_in, mu_in, grav_in, t_end
    real(dp), dimension(:), intent(in) :: VMR_bg_in

    ! Input/Output tracer values
    real(dp), intent(inout) :: q_0
    real(dp), dimension(:), intent(inout) :: q_v, q_1s

    integer :: ncall, n

    ! Time controls
    real(dp) :: t_now

    ! DLSODE variables
    integer :: n_eq
    real(dp), allocatable, dimension(:) :: y
    real(dp), allocatable, dimension(:) :: rwork
    integer, allocatable, dimension(:) :: iwork
    integer :: itol, itask, istate, iopt, mf
    integer :: rworkdim, iworkdim
    real(dp) :: rtol, atol

    !! Work variables
    integer :: n_bg
    real(dp), allocatable, dimension(:) :: VMR_bg
    real(dp) :: eta_g, bot, top

    n_dust = size(sp)

    if (first_call .eqv. .True.) then
      call mini_cloud_init(n_dust, sp)
      first_call = .False.
    end if


    !! Find the number density of the atmosphere
    T = T_in             ! Convert temperature to K
    p = P_in * 10.0_dp   ! Convert pascal to dyne cm-2

    !! Change mu_in to mu
    mu = mu_in ! Convert mean molecular weight to mu [g mol-1]

    !! Change gravity to cgs [cm s-2]
    grav = grav_in * 100.0_dp

    !! Number density [cm-3] of layer
    nd_atm = p/(kb*T)  

    !! Mass density of layer
    rho = (p*mu*amu)/(kb * T) ! Mass density [g cm-3]

    !! Specific gas constant of layer [erg g-1 K-1]
    Rd = Rgas/mu

    !! Allocate bg gas VMR
    n_bg = size(VMR_bg_in)
    allocate(VMR_bg(n_bg))
    VMR_bg(:) = VMR_bg_in(:)

    !! Calculate dynamical viscosity for this layer - do square root mixing law from Rosner 2012
    call eta_construct(n_bg, sp_bg)
    top = 0.0_dp
    bot = 0.0_dp
    do n = 1, n_bg
      eta_g = (5.0_dp/16.0_dp) * (sqrt(pi*(molg_g(n)*amu)*kb*T)/(pi*d_g(n)**2)) &
        & * ((((kb*T)/LJ_g(n))**(0.16_dp))/1.22_dp)
      top = top + sqrt(molg_g(n))*VMR_bg(n)*eta_g
      bot = bot + sqrt(molg_g(n))*VMR_bg(n)
    end do

    !! Mixture dynamical viscosity
    eta = top/bot

    !! Mixture kinematic viscosity
    nu = eta/rho

    !! dissipation of turbulent kinetic energy
    eps_d = nu**3/l_k**4

    !! Calculate mean free path for this layer
    mfp = (2.0_dp*eta/rho) * sqrt((pi * mu)/(8.0_dp*Rgas*T))

    !! Find saturation vapour pressure for each species
    call p_vap_sp(n_dust, T)

    !! Thermal velocity of gas species
    d_sp(:)%vth = sqrt((kb*T)/(2.0_dp*pi*d_sp(:)%m0))

    !! Gaseous diffusion constant of gas species
    d_sp(:)%D = 5.0_dp/(16.0_dp*N_A*d_sp(:)%d0**2*rho) * &
      & sqrt((Rgas*T*mu)/(2.0_dp*pi) * (d_sp(:)%mw + mu)/d_sp(:)%mw)

    !! Surface tension of each material
    call surface_tension(n_dust, T)

    ! -----------------------------------------
    ! ***  parameters for the DLSODE solver  ***
    ! -----------------------------------------

    itask = 1
    istate = 1
    iopt = 1

    n_eq = 1 + n_dust + n_dust

    ! Problem is stiff (usual)
    ! mf = 21 - full jacobian matrix with jacobian save
    ! mf = 22 - internal calculated jacobian
    mf = 22
    rworkdim = 22 + 9*n_eq + n_eq**2
    iworkdim = 20 + n_eq
    allocate(rwork(rworkdim), iwork(iworkdim))

    itol = 1
    rtol = 1.0e-3_dp           ! Relative tolerances for each scalar
    atol = 1.0e-30_dp               ! Absolute tolerance for each scalar (floor value)

    rwork(:) = 0.0_dp
    iwork(:) = 0

    rwork(1) = 0.0_dp               ! Critical T value (don't integrate past time here)
    rwork(5) = 0.0_dp              ! Initial starting timestep (start low, will adapt in DVODE)
    rwork(6) = 0.0_dp       ! Maximum timestep

    iwork(5) = 0               ! Max order required
    iwork(6) = 100000               ! Max number of internal steps
    iwork(7) = 1                ! Number of error messages

    allocate(y(n_eq))

    !! Give tracer values to y
    y(1) = q_0
    y(2:2+n_dust-1) = q_1s(:)
    y(2+n_dust:2+n_dust+n_dust-1) = q_v(:)

    !! Limit y values
    y(:) = max(y(:),1e-30_dp)

    t_now = 0.0_dp

    ! Set the printing flag
    ! 0 = no printing, 1 = printing
    call xsetf(1)

    ncall = 0

    do while ((t_now < t_end) .and. (ncall < 100))

      call DLSODE (RHS_mom, n_eq, y, t_now, t_end, itol, rtol, atol, itask, &
        & istate, iopt, rwork, rworkdim, iwork, iworkdim, jac_dum, mf)

      ncall = ncall + 1

      if (mod(ncall,10) == 0) then
        istate = 1
      else  if (istate == -1) then
        istate = 2
      else if (istate < -1) then
        print*, 'dlsode: ', istate
        exit
      end if

    end do

    !! Limit y values
    y(:) = max(y(:),1e-30_dp)

    !! Give y values to tracers
    q_0 = y(1)
    q_1s(:) = y(2:2+n_dust-1)
    q_v(:) = y(2+n_dust:2+n_dust+n_dust-1)

    deallocate(y, rwork, iwork, d_g, LJ_g, molg_g)

  end subroutine mini_cloud_2_mix

  subroutine RHS_mom(n_eq, time, y, f)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), intent(inout) :: time
    real(dp), dimension(n_eq), intent(inout) :: y
    real(dp), dimension(n_eq), intent(inout) :: f

    real(dp) :: w_sh2, w_acc2, w_co2
    real(dp) :: f_coal, f_coag, f_turb
    real(dp) :: m_c, r_c, Kn, beta, vf

    real(dp) :: n_d, rho_t, rho_d
    real(dp), dimension(n_dust) :: m_c_s
    real(dp), dimension(n_dust) :: J_hom, J_het, J_evap
    real(dp), dimension(n_dust) :: rho_s, rho_v, p_v, sat, n_v, dmdt
    real(dp), dimension(n_dust) :: V_frac, V_c_s

    !! In this routine, you calculate the instanenous new fluxes (f) for each moment
    !! The current values of each moment (y) are typically kept constant
    !! Basically, you solve for the RHS of the ODE for each moment

    !! Limit y values
    !y(:) = max(y(:),1e-99_dp)

    !! Convert y to real physical numbers to calculate f
    n_d = y(1) * nd_atm ! Convert to real number density
    rho_s(:) = y(2:2+n_dust-1) * rho   ! Convert to real mass density
    rho_v(:)= y(2+n_dust:2+n_dust+n_dust-1) * rho   ! Convert to real number density

    !! Find the true vapour VMR
    p_v(:) = rho_v(:) * Rd * T     !! Pressure of vapour
    n_v(:) = p_v(:)/(kb*T)        !! Number density of vapour

    !! Mean mass of particle
    rho_t = sum(rho_s(:))
    m_c = max(rho_t/n_d, m_seed)
    !! Find weighted bulk density of mixed grain composition using volume fraction
    V_frac(:) = max(rho_s(:)/rho_t,1e-99_dp)
    rho_d = sum(V_frac(:)*d_sp(:)%rho)

    ! !! Mean mass of particle per species
    ! m_c_s(:) = rho_s(:)/n_d
    ! !! Mean volumes of particle per species
    ! V_c_s(:) = m_c_s(:)/d_sp(:)%rho
    ! !! Volume fraction of species in grain
    ! V_frac(:) = V_c_s(:)/sum(V_c_s(:))
    ! !! Volume weighted bulk density
    ! rho_d = sum(V_frac(:)*d_sp(:)%rho)
    ! !! Mean particle mass [g]
    ! rho_t = sum(rho_s(:))
    ! m_c = max(rho_t/n_d,m_seed)

    !! Mass weighted mean radius of particle
    r_c = max(((3.0_dp*m_c)/(4.0_dp*pi*rho_d))**(third), r_seed)

    !! Knudsen number
    Kn = mfp/r_c

    !! Cunningham slip factor
    beta = 1.0_dp + Kn*(1.257_dp + 0.4_dp * exp(-1.1_dp/Kn))

    !! Settling velocity
    vf = (2.0_dp * beta * grav * r_c**2 * rho_d)/(9.0_dp * eta) & 
      & * (1.0_dp &
      & + ((0.45_dp*grav*r_c**3*rho*rho_d)/(54.0_dp*eta**2))**(0.4_dp))**(-1.25_dp)

    !! Find supersaturation ratio
    sat(:) = p_v(:)/d_sp(:)%p_vap

    !! Calculate condensation rate
    call calc_cond(n_dust, r_c, Kn, n_v, sat, V_frac, dmdt)

    !! Calculate homogenous nucleation rate
    call calc_hom_nuc(n_dust, sat, n_v, J_hom)

    !! Calculate hetrogenous nucleation rate
    call calc_het_nuc(n_dust, n_d, sat, r_c, J_het)
    J_het(:) = 0.0_dp

    !! Calculate seed particle evaporation rate
    call calc_seed_evap(n_dust, n_d, m_c, dmdt, J_evap)

    !! Calculate the coagulation rate
    call calc_coag(n_d, m_c, r_c, beta, f_coag)

    !! Calculate the coalesence rate
    call calc_coal(n_d, r_c, Kn, vf, f_coal)

    !! Calculate the turbulent shear collision velocity
    call calc_turb_shear(r_c, w_sh2)

    !! Calculate the turbulent acceleration collision velocity
    call calc_turb_acc(vf, rho_d, w_acc2)

    !! Calculate the turbulent fluid coupling collision velocity
    call calc_turb_couple(r_c, vf, rho_d, w_co2)

    !! Combine turbulent velocities into collision rate using total kernel
    f_turb = -pi * (2.0_dp*r_c)**2 * sqrt(2.0_dp/pi) * sqrt(w_sh2 + w_acc2 + w_co2) * n_d**2

    !! Calculate final net flux rate for each moment and vapour
    f(1) = sum(J_hom(:) + J_het(:) + J_evap(:)) + f_coag + f_coal + f_turb
    f(2:2+n_dust-1) = m_seed*(J_hom(:) + J_het(:) + J_evap(:)) + dmdt(:)*n_d
    f(2+n_dust:2+n_dust+n_dust-1) = -f(2:2+n_dust-1)

    !! Convert f to ratios
    f(1) = f(1)/nd_atm
    f(2:2+n_dust-1) = f(2:2+n_dust-1)/rho
    f(2+n_dust:2+n_dust+n_dust-1) = f(2+n_dust:2+n_dust+n_dust-1)/rho

  end subroutine RHS_mom

  !! Condensation and evaporation
  subroutine calc_cond(n_dust, r_c, Kn, n_v, sat, V_frac, dmdt)
    implicit none

    integer, intent(in) :: n_dust
    real(dp), intent(in) :: r_c, Kn
    real(dp), dimension(n_dust) :: n_v, sat, V_frac

    real(dp), dimension(n_dust), intent(out) :: dmdt

    integer :: n
    real(dp) :: Ak

    do n = 1, n_dust

      Ak = exp((2.0_dp*d_sp(n)%V0*d_sp(n)%sig)/(kb * T * r_c))

      if (Kn >= 1.0_dp) then
        !! Kinetic regime [g s-1]
        dmdt(n) = 4.0_dp * pi * r_c**2 * d_sp(n)%vth * d_sp(n)%m0 * n_v(n) * (1.0_dp - Ak/sat(n)) * V_frac(n)
      else
        !! Diffusive limited regime [g s-1]
        dmdt(n) = 4.0_dp * pi * r_c * d_sp(n)%D * d_sp(n)%m0 * n_v(n) * (1.0_dp - Ak/sat(n)) * V_frac(n)
      end if

    end do

  end subroutine calc_cond

  !! Classical nucleation theory (CNT)
  subroutine calc_hom_nuc(n_dust, sat, n_v, J_hom)
    implicit none

    integer, intent(in) :: n_dust
    real(dp), dimension(n_dust), intent(in) :: sat, n_v

    real(dp), dimension(n_dust), intent(out) :: J_hom

    integer :: n
    real(dp), parameter :: alpha = 1.0_dp
    real(dp) :: ln_ss, theta_inf, N_inf, N_star, N_star_1, dg_rt, Zel, tau_gr
    real(dp) :: f0, kbT


    do n = 1, n_dust

      if (d_sp(n)%inuc == 0) then
        J_hom(n) = 0.0_dp
        cycle
      end if

      if (sat(n) > 1.0_dp) then

        ! Efficency Variables
        ln_ss = log(sat(n)) ! Natural log of saturation ratio
        f0 = 4.0_dp * pi * d_sp(n)%r0**2 ! Monomer Area
        kbT = kb * T         ! kb * T

        !! Find Nstar, theta_inf -> Dg/RT (Eq. 11 Lee et al. (2015a))
        theta_inf = (f0 * d_sp(n)%sig)/(kbT)  !Theta_infty eq.8 (Lee et al. (2015a)
        N_inf = (((twothird) * theta_inf) / ln_ss)**3

        !! Gail et al. (2014) ! note, typo in Lee et al. (2015a)
        N_star = 1.0_dp + (N_inf / 8.0_dp) &
          & * (1.0_dp + sqrt(1.0_dp + 2.0_dp*(d_sp(n)%Nf/N_inf)**third) &
          & - 2.0_dp*(d_sp(n)%Nf/N_inf)**third)**3
        N_star = max(1.00001_dp, N_star) ! Ensure no div 0
        N_star_1 = N_star - 1.0_dp

        !! delta G (N_star) / RT (Eq. 11 Lee et al. (2015a))
        dg_rt = theta_inf * (N_star_1 / (N_star_1**third + d_sp(n)%Nf**third))

        !! Calculate Zeldovich factor at N_star (Eq. 9 Lee et al. (2015a))
        Zel = sqrt((theta_inf / (9.0_dp * pi * (N_star_1)**(4.0_dp/3.0_dp))) &
          & * ((1.0_dp + 2.0_dp*(d_sp(n)%Nf/N_star_1)**third)/(1.0_dp + (d_sp(n)%Nf/N_star_1)**third)**3))

        !! Calculate !!inverse!! of tau_gr
        tau_gr = (f0 * N_star**(twothird)) * alpha * sqrt(kbT &
          & / (2.0_dp * pi * d_sp(n)%mw * amu)) * n_v(n)

        !! Finally calculate J_star [cm-3 s-1] ! Note underfloat limiter here
        J_hom(n) = n_v(n) * tau_gr * Zel * exp(max(-300.0_dp, N_star_1*ln_ss - dg_rt))
      else 
        !! Unsaturated, zero nucleation
        J_hom(n) = 0.0_dp
      end if

    end do
    
  end subroutine calc_hom_nuc

  !! Classical hetrogenous nucleation theory
  subroutine calc_het_nuc(n_dust, n_d, sat, r_c, J_het)
    implicit none

    integer, intent(in) :: n_dust
    real(dp), intent(in) :: r_c, n_d
    real(dp), dimension(n_dust), intent(in) :: sat

    real(dp), dimension(n_dust), intent(out) :: J_het

    integer :: n
    real(dp) :: r_crit, FG, Theta, Zel, Nl
    real(dp) :: f, mu_con, phi, f0, x
    real(dp) :: c_surf, nu, F_des

    do n = 1, n_dust

      if (d_sp(n)%inuc == 0) then
        J_het(n) = 0.0_dp
        cycle
      end if

      if (sat(n) > 1.0_dp) then

        !! Critical particle size
        r_crit = (2.0_dp * d_sp(n)%mw * d_sp(n)%sig)/(d_sp(n)%rho*Rgas*T*log(sat(n)))

        !! Formation energy of critical cluster - classical approximation
        FG = 4.0_dp/3.0_dp * pi * d_sp(n)%sig * r_crit**2

        !! Cluster-cluster diffusive interaction rate 
        Theta = p/(sqrt(2.0_dp*d_sp(n)%m0*kb*T))

        !! Number of vapour units in critical cluster size
        Nl = max((4.0_dp/3.0_dp * pi * r_crit**3)/d_sp(n)%V0, 1.0_dp)

        !! Zeldovich factor
        Zel = sqrt(FG/(3.0_dp*pi*kb*T*Nl**2))

        !! Cosine contact angle from surface tension
        mu_con = 0.5_dp

        !! Begin shape factor calculation
        x = r_c/r_crit
        phi = sqrt(1.0_dp - 2.0_dp*mu_con*x + x**2)
        f0 = (x - mu_con)/phi

        !! Shape factor
        f = 0.5_dp * (1.0_dp + ((1.0_dp-mu_con*x)/phi)**3 + x**3*(2.0_dp - 3.0_dp*f0 + f0**3) &
          & + 3.0_dp*mu_con*x**2*(f0 - 1.0_dp))

        !! Desorption energy
        F_des = 0.5_dp * eV !~ 0.5 eV (Gao & Powell 2021)

        !! Oscillation frequency
        nu = 10e11_dp  * sqrt(F_des/(kb*d_sp(n)%m0)) ! Gao per. comm.

        !! Number density of condensate molecules on the nucleating surface
        c_surf = Theta/nu * exp(F_des/(kb*T))
        
        !! Heterogenous nucleation rate [cm-3 s-1]
        J_het(n) = 4.0_dp * pi**2 * r_c**2 * r_crit**2 * Theta * c_surf * Zel * exp(-(FG*f)/(kb*T)) * n_d

      else
        !! Unsaturated, zero nucleation
        J_het(n) = 0.0_dp
      end if

    end do

  end subroutine calc_het_nuc

  !! Seed particle evaporation
  subroutine calc_seed_evap(n_dust, n_d, m_c, dmdt, J_evap)
    implicit none

    integer, intent(in) :: n_dust
    real(dp), intent(in) :: n_d, m_c
    real(dp), dimension(n_dust), intent(in) :: dmdt

    real(dp), dimension(n_dust), intent(out) :: J_evap

    integer :: n
    real(dp) :: tau_evap

    do n = 1, n_dust

      if (d_sp(n)%inuc == 0) then
        J_evap(n) = 0.0_dp
        cycle
      end if

      if ((dmdt(n) >= 0.0_dp) .or. (n_d/nd_atm <= 1e-29_dp)) then

        !! If growing or too little number density then evaporation can't take place
        J_evap(n) = 0.0_dp

      else 

        !! Check if average mass is around 0.1% the seed particle mass
        !! This means the core is (probably) exposed to the air and can evaporate freely
        if (m_c <= (1.001_dp * m_seed)) then
          tau_evap = 1.0_dp !max(m_c/abs(f_cond),1.0_dp)
          !! Seed particle evaporation rate [cm-3 s-1]
          J_evap(n) = -n_d/tau_evap
        else
          !! There is still some mantle to evaporate from
          J_evap(n) = 0.0_dp
        end if

      end if

    end do

  end subroutine calc_seed_evap

  !! Particle-particle Brownian coagulation
  subroutine calc_coag(n_d, m_c, r_c, beta, f_coag)
    implicit none

    real(dp), intent(in) :: n_d
    real(dp), intent(in) :: m_c, r_c, beta

    real(dp), intent(out) :: f_coag

    real(dp) :: phi, del_r, D_r, V_r, lam_r

    !! Particle diffusion rate
    D_r = (kb*T*beta)/(6.0_dp*pi*eta*r_c)

    !! Thermal velocity limit rate
    V_r = sqrt((8.0_dp*kb*T)/(pi*m_c))

    !! Ratio fraction
    lam_r = (8.0_dp*D_r)/(pi*V_r)

    !! Interpolation function
    del_r = ((2.0_dp*r_c + lam_r)**3 - (4.0_dp*r_c**2 + lam_r**2)**(1.5_dp))/(6.0_dp*r_c*lam_r) &
      & - 2.0_dp*r_c 

    !! Correction factor
    phi = 2.0_dp*r_c/(2.0_dp*r_c + sqrt(2.0_dp)*del_r) + (4.0_dp*D_r)/(r_c*sqrt(2.0_dp)*V_r) 

    !! Coagulation flux (Zeroth moment) [cm-3 s-1]
    f_coag = -(4.0_dp * kb * T * beta)/(3.0_dp * eta * phi)  * n_d**2

  end subroutine calc_coag

  !! Particle-particle gravitational coalesence
  subroutine calc_coal(n_d, r_c, Kn, vf, f_coal)
    implicit none

    real(dp), intent(in) :: n_d
    real(dp), intent(in) :: r_c, Kn, vf

    real(dp), intent(out) :: f_coal

    real(dp) :: d_vf,  Stk, E
    real(dp), parameter :: eps = 0.5_dp

    !! Estimate differential velocity
    d_vf = eps * vf

    !! Calculate E
    if (Kn >= 1.0_dp) then
      !! E = 1 when Kn > 1
      E = 1.0_dp
    else
      !! Calculate Stokes number
      Stk = (vf * d_vf)/(grav * r_c)
      E = max(0.0_dp,1.0_dp - 0.42_dp*Stk**(-0.75_dp))
    end if

    !! Coalesence flux (Zeroth moment) [cm-3 s-1]
    f_coal = -2.0_dp*pi*r_c**2*d_vf*E*n_d**2

  end subroutine calc_coal

  !! Particle-particle turbulent shear collisions
  subroutine calc_turb_shear(r_c, w_sh2)
    implicit none

    real(dp), intent(in) :: r_c

    real(dp), intent(out) :: w_sh2

    !! Shear relative velocity (Wang et al. 1998)
    w_sh2 = 1.0_dp/15.0_dp * (2.0_dp*r_c)**2 * eps_d/nu

  end subroutine calc_turb_shear

  !! Particle-particle turbulent acceleration collisions
  subroutine calc_turb_acc(vf, rho_d, w_acc2)
    implicit none

    real(dp), intent(in) :: vf, rho_d

    real(dp), intent(out) :: w_acc2

    real(dp), parameter :: gam = 10.0_dp
    real(dp) :: u2, tau_i, b, T_l, th_i, c1, c2

    !! Follows Wang et al. (2001) 
    !! - see also Park et el. (2002) and Kruis & Kusters (1997)

    u2 = (gam * sqrt(eps_d*nu))/0.183_dp
    b = (3.0_dp*rho)/(2.0_dp*rho_d + rho)
    T_L = (0.4_dp*u2)/eps_d

    tau_i = vf/grav
    th_i = tau_i/T_L

    c1 = sqrt((1.0_dp + th_i + th_i)/((1.0_dp + th_i)*(1.0_dp + th_i)))
    c2 = (1.0_dp/(((1.0_dp + th_i)*(1.0_dp + th_i))) - 1.0_dp/(((1.0_dp + gam*th_i)*(1.0_dp + gam*th_i))))

    w_acc2 = 3.0_dp*(1.0_dp-b)**2*u2*(gam/(gam-1.0_dp)) & 
      & * (((th_i + th_i)**2 - 4.0_dp*th_i*th_i*c1)/(th_i + th_i)) * c2

  end subroutine calc_turb_acc

  !! Particle-particle turbulent fluid coupling collisions
  subroutine calc_turb_couple(r_c, vf, rho_d, w_co2)
    implicit none

    real(dp), intent(in) :: r_c, vf, rho_d

    real(dp), intent(out) :: w_co2

    real(dp), parameter :: lam_d = 10.0_dp
    real(dp) :: dudt2, tau

    !! Fluid coupling turbulent term (Wang et al. 1998)
    tau = vf/grav
    dudt2 = 1.16_dp * eps_d**(1.5_dp)/sqrt(nu)
    w_co2 = 2.0_dp*(1.0_dp - rho/rho_d)**2 * tau**2 * dudt2 * ((2.0_dp*r_c)**2/lam_d**2)

  end subroutine calc_turb_couple

  !! Dummy jacobian subroutine required for dlsode
  subroutine jac_dum (NEQ, X, Y, ML, MU, PD, NROWPD)
    integer, intent(in) :: NEQ, ML, MU, NROWPD
    real(dp), intent(in) :: X
    real(dp), dimension(NEQ), intent(in) :: Y
    real(dp), dimension(NROWPD, NEQ), intent(inout) :: PD
  end subroutine jac_dum

end module mini_cloud_2_mix_mod