program test_mini_cloud_2
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  use mini_cloud_2_exp_mix_mod, only : mini_cloud_2_exp_mix
  use mini_cloud_vf_mod, only : mini_cloud_vf
  use mini_cloud_opac_mie_mod, only : opac_mie
  use vert_diff_exp_mod, only : vert_diff_exp
  use vert_adv_exp_mod, only : vert_adv_exp
  use vert_diff_imp_mod, only : vert_diff_imp
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: kb = 1.380649e-16_dp ! erg K^-1 - Boltzmann's constant
  real(dp), parameter :: amu = 1.66053906660e-24_dp ! g - Atomic mass unit
  real(dp), parameter :: r_seed = 1e-7_dp
  real(dp), parameter :: V_seed = 4.0_dp/3.0_dp * pi * r_seed**3

  integer :: example, n_it
  real(dp) :: t_step, time

  integer :: nlay, nlev, i, n, u
  character(len=20), allocatable, dimension(:) :: sp
  character(len=20), allocatable, dimension(:) :: sp_bg
  real(dp), allocatable, dimension(:) :: Tl, pl, mu, Kzz, pe, nd_atm, rho, rho_d, mol_w_sp, mol_w_v, cp, dTdt
  real(dp), allocatable, dimension(:) :: q_0, r_c, m_c, q0, r_c_old, del
  real(dp), allocatable, dimension(:,:) :: VMR, q, q_v, q_1, vf
  real(dp) :: grav, met

  integer :: n_wl
  real(dp), allocatable, dimension(:) :: wl_e, wl
  real(dp), allocatable, dimension(:,:) :: k_ext, ssa, g

  integer :: nlines, io, idx, idx1
  real(dp) :: p_bot, p_top
  real(dp), allocatable, dimension(:) :: T_f, p_f
  real(dp) :: m_seed, Nc, rho_c_tot, rho_d_mean

  integer :: j, nsp
  real(dp) :: T_eff, k_ir, tau

  logical :: end

  !! time step
  t_step = 100.0_dp

  !! Number of iterations
  n_it = 10000

  !! Start time
  time = 6840.0_dp

  !! example select
  example = 3

  select case (example)

    case(2)

      !! In this example we perform a 1D test of mini_cloud with diffusion and settling included

      nsp = 2

      nlay = 192
      nlev = nlay + 1

      n_wl = 11
      allocate(wl_e(n_wl+1), wl(n_wl), k_ext(nlay,n_wl), ssa(nlay,n_wl), g(nlay,n_wl))

      !! Wavelengths to calculate opacity
      wl_e = (/0.260, 0.420, 0.610, 0.850, 1.320, 2.020,2.500,3.500,4.400,8.70,20.00,324.68 /)
      wl(:) = (wl_e(2:n_wl+1) +  wl_e(1:n_wl))/ 2.0_dp

      !! Allocate all variables, set constant values
      allocate(Tl(nlay), pl(nlay), pe(nlev), mu(nlay), Kzz(nlay), nd_atm(nlay), rho(nlay), cp(nlay), dTdt(nlay))

      !! Find pressure level grid - logspaced between p_top and p_bot
      p_top = 3e-3_dp * 1e5_dp
      p_bot = 300.0_dp * 1e5_dp

      p_top = log10(p_top)
      p_bot = log10(p_bot)
      do i = 1, nlev
        pe(i) = 10.0_dp**((p_bot-p_top) * real(i-1,dp) / real(nlev-1,dp) + p_top) 
      end do
      do i = 1, nlay
      pl(i) = (pe(i+1) - pe(i)) / log(pe(i+1)/pe(i))
      end do
      p_top = 10.0_dp**p_top
      p_bot = 10.0_dp**p_bot

      !! Read T-p file and interpolate T
      open(newunit=u,file='Y_400K_paper/Gao_2018_400_425.txt',action='read')
      ! Read header
      read(u,*) ; read(u,*)
    ! Find number of lines in file
      nlines = 0
      do
        read(u,*,iostat=io)
        if (io /= 0) then 
            exit
        end if
        nlines = nlines + 1
      end do
      ! Allocate values for T-p profile
      allocate(T_f(nlines),p_f(nlines))
      ! Rewind file
      rewind(u)
      ! Read header again
      read(u,*); read(u,*)
      do i = 1, nlines
        read(u,*) T_f(i), p_f(i)
      end do
      p_f(:) = p_f(:)*1e5_dp ! Convert to pa
      close(u)
      ! Interpolate to pressure grid
      do i = 1, nlay
        call locate(p_f(:), nlines, pl(i), idx)
        if (idx < 1) then
          Tl(i) = T_f(1)
        else if (idx >= nlines) then
          Tl(i) = T_f(nlines)
        else
          idx1 = idx + 1
          call linear_interp(pl(i), p_f(idx), p_f(idx1), T_f(idx), T_f(idx1), Tl(i))
        end if
      end do

      close(u)

      !! Assume constant Kzz [cm2 s-1]
      Kzz(:) = 1e8_dp

      !! Print T-p-Kzz profile
      print*, 'i, pl [bar], T[k], Kzz [cm2 s-1]'
      do i = 1, nlay
        print*, i, pl(i)/1e5_dp, Tl(i), Kzz(i)
      end do

      !! Assume constant background gas mean molecular weight [g mol-1] @ approx solar
      mu(:) = 2.33_dp

      !! Assume constant gravity [m s-2]
      grav = (10.0_dp**(4.25_dp))/100.0_dp

      !! Number density [cm-3] of layer
      nd_atm(:) = (pl(:)*10.0_dp)/(kb*Tl(:))  

      !! Mass density of layer
      rho(:) = (pl(:)*10.0_dp*mu(:)*amu)/(kb * Tl(:)) ! Mass density [g cm-3]

      !! Assume constant H2, He and H background VMR @ approx solar
      allocate(VMR(nlay,2),sp_bg(2))
      sp_bg = (/'H2','He'/)
      VMR(:,1) = 0.85_dp
      VMR(:,2) = 0.15_dp

      !! Assumed condensate species
      allocate(sp(nsp), rho_d(nsp), mol_w_sp(nsp), mol_w_v(nsp))
      sp = (/'KCl', 'ZnS'/)
      rho_d = (/1.99_dp, 3.9_dp/)
      mol_w_sp = (/74.5513_dp, 97.4450_dp/)
      mol_w_v = (/74.5513_dp, 97.4450_dp/)

      m_seed = V_seed * rho_d(1)

      allocate(q_v(nlay,nsp), q_0(nlay), q_1(nlay,nsp), q0(nsp*2+1), q(nlay,nsp*2+1))
      allocate(r_c(nlay), m_c(nlay), vf(nlay,2), r_c_old(nlay), del(nlay))

      q_v(:,:) = 1e-30_dp
      q_0(:) = 1e-30_dp
      q_1(:,:) = 1e-30_dp

      q0(1) = 1.17e-7_dp * mol_w_v(1)/mu(nlay)
      q0(2) = 3.63e-8_dp * mol_w_v(2)/mu(nlay)
      q0(3:) = 1e-30_dp

      time = 0.0_dp
      n = 0

      r_c_old(:) = r_seed * 1e4_dp

   
      do n = 1, n_it

        !$omp parallel do default(shared), private(i), schedule(dynamic)
        do i = 1, nlay

          !! Call mini-cloud and perform integrations for a single layer
          call mini_cloud_2_exp_mix(i, Tl(i), pl(i), grav, mu(i), met, cp(i), VMR(i,:), t_step, sp, sp_bg, & 
            & nsp, q_v(i,:), q_0(i), q_1(i,:), dTdt(i))

          !! Calculate settling velocity for this layer
          call mini_cloud_vf(Tl(i), pl(i), grav, mu(i), VMR(i,:), rho_d(:), sp_bg, & 
            &  nsp, q_0(i), q_1(i,:), vf(i,:))

          !! Calculate the opacity at the wavelength grid
          call opac_mie(nsp, sp, Tl(i), mu(i), pl(i), q_0(i), q_1(i,:), rho_d(:), n_wl, wl, k_ext(i,:), ssa(i,:), g(i,:))
        end do
        !$omp end parallel do

        q(:,1:nsp) = q_v(:,:)
        q(:,nsp+1) = q_0(:) * nd_atm(:) / rho(:) ! Make mass ratio for vertical transport
        q(:,nsp+2:) = q_1(:,:)

        call vert_adv_exp(nlay, nlev, t_step, mu, grav, Tl, pl, pe, vf, nsp+1, q(:,nsp+1:))

        call vert_diff_exp(nlay, nlev, t_step, mu, grav, Tl, pl, pe, Kzz, nsp*2+1, q(:,:), q0(:))

        q_v(:,:) = q(:,1:nsp)
        q_0(:) = q(:,nsp+1) * rho(:) / nd_atm(:) ! Return to number density after vertical transport
        q_1(:,:) = q(:,nsp+2:)

        do i = 1, nlay
          rho_c_tot = sum(q_1(i,:))*rho(i)
          Nc = q_0(i)*nd_atm(i)

          !! Mean mass of particle [g]
          m_c(i) = max(rho_c_tot/Nc,m_seed)

          rho_d_mean = 0.0_dp
          do j = 1, nsp
            rho_d_mean = rho_d_mean + (q_1(i,j)*rho(i))/rho_c_tot * rho_d(j)
          end do

          !! Mass weighted mean radius of particle [um]
          r_c(i) = max(((3.0_dp*m_c(i))/(4.0_dp*pi*rho_d_mean))**(1.0_dp/3.0_dp),r_seed) * 1e4_dp
        end do

        !! increment time
        time = time + t_step

        print*, n, time

        end = .True.

        do i = 1, nlay
          del(i) = abs(r_c_old(i) - r_c(i))/r_c_old(i)
          if ((del(i) > 1e-4_dp) .or.  (del(i)/t_step > 1e-4_dp)) then
            end = .False.
            exit
          end if
        end do
        r_c_old(:) = r_c(:)

        !! increment time
        time = time + t_step

        print*, n, time, maxval(del(:)/t_step)

        if ((end .eqv. .True.) .and. (n > int(1e5))) then
          print*, 'exit: ', n, n_it, end
          exit
        end if


      end do

      print*, del(:)
      print*, del(:)/t_step

      do i = 1, nlay
        print*, i, pl(i)/1e5_dp, Tl(i), Kzz(i), q_v(i,:), q_0(i), q_1(i,:), r_c(i), vf(i,1), q_0(i)*nd_atm(i)
      end do

      !! mini-cloud test output
      call output(n, time, nlay)

    case(3)

      !! L-T transition dwarf case (TiO2, Al2O3, Fe + Mg2SiO4/MgSiO3) 

      nsp = 4

      nlay = 192
      nlev = nlay + 1

      n_wl = 11
      allocate(wl_e(n_wl+1), wl(n_wl), k_ext(nlay,n_wl), ssa(nlay,n_wl), g(nlay,n_wl))

      !! Wavelengths to calculate opacity
      wl_e = (/0.260, 0.420, 0.610, 0.850, 1.320, 2.020,2.500,3.500,4.400,8.70,20.00,324.68 /)
      wl(:) = (wl_e(2:n_wl+1) +  wl_e(1:n_wl))/ 2.0_dp

      !! Allocate all variables, set constant values
      allocate(Tl(nlay), pl(nlay), pe(nlev), mu(nlay), Kzz(nlay), nd_atm(nlay), rho(nlay), cp(nlay), dTdt(nlay))

      !! Find pressure level grid (pa) - logspaced between p_top and p_bot
      p_top = 1e-4_dp * 1e5_dp
      p_bot = 1000.0_dp * 1e5_dp

      p_top = log10(p_top)
      p_bot = log10(p_bot)
      do i = 1, nlev
        pe(i) = 10.0_dp**((p_bot-p_top) * real(i-1,dp) / real(nlev-1,dp) + p_top) 
      end do
      do i = 1, nlay
      pl(i) = (pe(i+1) - pe(i)) / log(pe(i+1)/pe(i))
      end do
      p_top = 10.0_dp**p_top
      p_bot = 10.0_dp**p_bot

      !! Assume constant background gas mean molecular weight [g mol-1] @ approx solar
      mu(:) = 2.33_dp

      !! Approx [M/H] log10 Solar metallicity of atmosphere
      met = 0.0_dp

      !! Assume constant gravity [m s-2]
      grav = (10.0_dp**(4.5_dp))/100.0_dp

      !! Find T-p profile - assume semi-grey atmosphere and use Eddington approximation
      T_eff = 1300.0_dp!1300.0_dp ! [K]
      k_ir = 1e-2_dp ! [cm^2 g-1]
      do i = 1, nlay
        tau = (k_ir * (pl(i) * 10.0_dp))/(grav*100.0_dp)
        Tl(i) = 3.0_dp/4.0_dp * T_eff**4 * (2.0_dp/3.0_dp + tau)
        Tl(i) = Tl(i)**(1.0_dp/4.0_dp)
      end do 

      !! Number density [cm-3] of layer
      nd_atm(:) = (pl(:)*10.0_dp)/(kb*Tl(:))  

      !! Mass density of layer
      rho(:) = (pl(:)*10.0_dp*mu(:)*amu)/(kb * Tl(:)) ! Mass density [g cm-3]

      !! Assume constant Kzz [cm2 s-1]
      Kzz(:) = 1e8_dp

      !! Heat capacity of atmosphere [J kg-1 K-1]
      cp(:) = 1.3e4_dp

      !! Change in temperature due to latent heat
      dTdt(:) = 0.0_dp

      !! Print T-p-Kzz profile
      print*, 'i, pl [bar], T[k], Kzz [cm2 s-1]'
      do i = 1, nlay
        print*, i, pl(i)/1e5_dp, Tl(i), Kzz(i)
      end do

      !! Assume constant H2, He and H background VMR @ approx solar
      allocate(VMR(nlay,2),sp_bg(2))
      sp_bg = (/'H2','He'/)
      VMR(:,1) = 0.85_dp
      VMR(:,2) = 0.15_dp

      !! Assumed condensate species
      allocate(sp(nsp), rho_d(nsp), mol_w_sp(nsp), mol_w_v(nsp))
      sp = (/'TiO2   ', 'Al2O3  ', 'Fe     ', 'Mg2SiO4'/)
      rho_d = (/4.23_dp, 3.986_dp, 7.874_dp, 3.21_dp/)
      mol_w_sp = (/79.8658_dp, 101.961_dp, 55.8450_dp, 140.693_dp/)
      mol_w_v = (/79.8658_dp, 26.98153860_dp, 55.8450_dp, 24.305_dp/)

      !! Seed particle mass (assume rho_d at first index)
      m_seed = V_seed * rho_d(1)

      allocate(q_v(nlay,nsp), q_0(nlay), q_1(nlay,nsp), q0(nsp*2+1), q(nlay,nsp*2+1))
      allocate(r_c(nlay), m_c(nlay), vf(nlay,nsp+1), r_c_old(nlay), del(nlay))

      !! Set everything to zero first
      do j = 1, nsp
        q_v(:,j) = 1e-30_dp
        q_0(:) = 1e-30_dp
        q_1(:,j) = 1e-30_dp
      end do

      !! Lower mass mixing ratio boundary conditions for vapour + cloud (VMR taken from Asplund et al. 2001)
      !! Assume cloud = zero boundary condition (all evaporated)
      q0(1) = 9.33e-8_dp * mol_w_v(1)/mu(nlay)
      q0(2) = 2.69e-6_dp * mol_w_v(2)/mu(nlay)
      q0(3) = 2.88e-5_dp * mol_w_v(3)/mu(nlay)
      q0(4) = 3.55e-5_dp * mol_w_v(4)/mu(nlay)
      q0(5:) = 1e-30_dp

      time = 0.0_dp
      n = 0

      r_c_old(:) = r_seed * 1e4_dp

      !! Start main iteration loops
      do n = 1, n_it

        !$omp parallel do default(shared), private(i), schedule(dynamic)
        do i = 1, nlay
          !! Calculate settling velocity for this layer
          call mini_cloud_vf(Tl(i), pl(i), grav, mu(i), VMR(i,:), rho_d(:), sp_bg, & 
            &  nsp, q_0(i), q_1(i,:), vf(i,1:2))
            vf(i,3:) = vf(i,2)
        end do
        !$omp end parallel do

        !! Combine everything q into a single 2D array for advection and diffusion
        q(:,1:nsp) = q_v(:,:)
        q(:,nsp+1) = q_0(:) 
        q(:,nsp+2:) = q_1(:,:)

        !! Vertical advection (settling tracers)
        call vert_adv_exp(nlay, nlev, t_step/2.0_dp, mu, grav, Tl, pl, pe, vf(:,:), nsp+1, q(:,nsp+1:))

        !! Vertical diffusion (diffused tracers)
        call vert_diff_imp(nlay, nlev, t_step/2.0_dp, mu, grav, Tl, pl, pe, Kzz, nsp*2+1, q(:,:), q0(:))

        q_v(:,:) = q(:,1:nsp)
        q_0(:) = q(:,nsp+1)
        q_1(:,:) = q(:,nsp+2:)

        !$omp parallel do default(shared), private(i), schedule(dynamic)
        do i = 1, nlay
          !! Call mini-cloud and perform integrations for a single layer
          call mini_cloud_2_exp_mix(i, Tl(i), pl(i), grav, mu(i), met, cp(i), VMR(i,:), t_step, sp, sp_bg, & 
            & nsp, q_v(i,:), q_0(i), q_1(i,:), dTdt(i))
        end do
        !$omp end parallel do

        q(:,1:nsp) = q_v(:,:)
        q(:,nsp+1) = q_0(:) 
        q(:,nsp+2:) = q_1(:,:)

        call vert_diff_imp(nlay, nlev, t_step/2.0_dp, mu, grav, Tl, pl, pe, Kzz, nsp*2+1, q(:,:), q0(:))

        q_v(:,:) = q(:,1:nsp)
        q_0(:) = q(:,nsp+1)
        q_1(:,:) = q(:,nsp+2:)

        !$omp parallel do default(shared), private(i), schedule(dynamic)
        do i = 1, nlay
          !! Calculate settling velocity for this layer
          call mini_cloud_vf(Tl(i), pl(i), grav, mu(i), VMR(i,:), rho_d(:), sp_bg, & 
            &  nsp, q_0(i), q_1(i,:), vf(i,1:2))
          vf(i,3:) = vf(i,2)
        end do
        !$omp end parallel do

        q(:,1:nsp) = q_v(:,:)
        q(:,nsp+1) = q_0(:) 
        q(:,nsp+2:) = q_1(:,:)

        call vert_adv_exp(nlay, nlev, t_step/2.0_dp, mu, grav, Tl, pl, pe, vf(:,:), nsp+1, q(:,nsp+1:))

        q_v(:,:) = max(q(:,1:nsp),1e-30_dp)
        q_0(:) = max(q(:,nsp+1),1e-30_dp)
        q_1(:,:) = max(q(:,nsp+2:),1e-30_dp)


        !$omp parallel do default(shared), private(i), schedule(dynamic)
        do i = 1, nlay
            !! Calculate the opacity at the wavelength grid
          call opac_mie(nsp, sp, Tl(i), mu(i), pl(i), q_0(i), q_1(i,:), rho_d(:), n_wl, wl, k_ext(i,:), ssa(i,:), g(i,:))
        end do
        !$omp end parallel do

        do i = 1, nlay
          !! Total condensed mass
          rho_c_tot = sum(q_1(i,:))*rho(i)

          !! Total number density of clouds
          Nc = q_0(i)*nd_atm(i)

          !! Mean mass of particle [g]
          m_c(i) = max(rho_c_tot/Nc,m_seed)

          !! Average bulk density
          rho_d_mean = 0.0_dp
          do j = 1, nsp
            rho_d_mean = rho_d_mean + (q_1(i,j)*rho(i))/rho_c_tot * rho_d(j)
          end do

          !! Mass weighted mean radius of particle [um]
          r_c(i) = max(((3.0_dp*m_c(i))/(4.0_dp*pi*rho_d_mean))**(1.0_dp/3.0_dp),r_seed) * 1e4_dp
        end do

        !! Check for convergence
        end = .True.

        do i = 1, nlay
          del(i) = abs(r_c_old(i) - r_c(i))/r_c_old(i)
          if ((del(i) > 1e-4_dp) .or.  (del(i)/t_step > 1e-4_dp)) then
            end = .False.
            !exit
          end if
        end do
        r_c_old(:) = r_c(:)

        !! increment time
        time = time + t_step

        print*, n, time, maxval(del(:)/t_step)

        if ((end .eqv. .True.) .and. (n > int(1e5))) then
          print*, 'exit: ', n, n_it, end
          !exit
        end if

      end do

      print*, del(:)
      print*, del(:)/t_step

      do i = 1, nlay
        print*, i, pl(i)/1e5_dp, Tl(i), Kzz(i), q_v(i,:), q_0(i), q_1(i,:), r_c(i), vf(i,1), q_0(i)*nd_atm(i)
      end do

      !! mini-cloud test output
      call output(n, time, nlay)

    case default 
      print*, 'Invalid test case example: ', example
      stop
    end select

contains

  subroutine output(t, time, nlay)
    implicit none
    integer, intent(in) :: t, nlay
    double precision, intent(in) :: time
    integer, save :: u1, u2, u3, u4
    logical, save :: first_call = .True.

    if (first_call .eqv. .True.) then
      open(newunit=u1,file='results_2_exp_mix/tracers.txt',action='readwrite')
      write(u1,*) nsp
      write(u1,*) rho_d(:)
      write(u1,*) mol_w_sp(:)
      open(newunit=u2,file='results_2_exp_mix/opac_k.txt',action='readwrite')
      write(u2,*) wl(:)
      open(newunit=u3,file='results_2_exp_mix/opac_a.txt',action='readwrite')
      write(u3,*) wl(:)
      open(newunit=u4,file='results_2_exp_mix/opac_g.txt',action='readwrite')
      write(u4,*) wl(:)            
      first_call = .False.
    end if

    do i = 1, nlay
      write(u1,*) t, time, Tl(i), pl(i), grav, mu(i), VMR(i,:), q_v(i,:), q_0(i), q_1(i,:), vf(i,1), dTdt(i)
      write(u2,*) t, time, pl(i), k_ext(i,:)
      write(u3,*) t, time, pl(i), ssa(i,:)
      write(u4,*) t, time, pl(i), g(i,:)
    end do

  end subroutine output

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

  subroutine linear_interp(xval, x1, x2, y1, y2, yval)
    implicit none

    real(kind=dp), intent(in) :: xval, y1, y2, x1, x2
    real(kind=dp), intent(out) :: yval
    real(kind=dp) :: norm

    if (x1 >= x2) then
      print*, 'Error in linear_interp: x1 >= x2 - STOPPING', x1, x2
      stop
    end if

    norm = 1.0_dp / (x2 - x1)

    yval = (y1 * (x2 - xval) + y2 * (xval - x1)) * norm

  end subroutine linear_interp

end program test_mini_cloud_2
