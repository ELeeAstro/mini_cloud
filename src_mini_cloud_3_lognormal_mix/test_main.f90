program test_mini_cloud_3
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  use mini_cloud_3_lognormal_mix_mod, only : mini_cloud_3_lognormal_mix
  use mini_cloud_vf_num_mod, only : mini_cloud_vf
  use mini_cloud_opac_mie_mod, only : opac_mie
  use vert_adv_exp_mod, only : vert_adv_exp
  use vert_diff_imp_mod, only : vert_diff_imp
  use cli_progress
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: amu = 1.66053906660e-24_dp
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
  real(dp), allocatable, dimension(:) :: q_2
  real(dp) :: grav, met

  integer :: n_wl
  real(dp), allocatable, dimension(:) :: wl_e, wl
  real(dp), allocatable, dimension(:,:) :: k_ext, ssa, g

  integer :: nlines, io, idx, idx1
  real(dp) :: p_bot, p_top
  real(dp), allocatable, dimension(:) :: T_f, p_f
  real(dp) :: m_seed, Nc, rho_c_tot, rho_d_mean

  integer :: nsp
  real(dp) :: T_eff, k_ir, tau
  real(dp), dimension(3) :: vf_tmp

  logical :: end

  call progress_begin()

  !! time step
  t_step = 100.0_dp

  !! Number of iterations
  n_it = 10000

  !! Start time
  time = 0.0_dp

  !! example select
  example = 3

  select case (example)

    case(2)

      !! Y-dwarf test case: KCl + ZnS

      nsp = 2

      nlay = 192
      nlev = nlay + 1

      n_wl = 11
      allocate(wl_e(n_wl+1), wl(n_wl), k_ext(nlay,n_wl), ssa(nlay,n_wl), g(nlay,n_wl))

      wl_e = (/0.260, 0.420, 0.610, 0.850, 1.320, 2.020,2.500,3.500,4.400,8.70,20.00,324.68 /)
      wl(:) = (wl_e(2:n_wl+1) + wl_e(1:n_wl))/ 2.0_dp

      allocate(Tl(nlay), pl(nlay), pe(nlev), mu(nlay), Kzz(nlay), nd_atm(nlay), rho(nlay), cp(nlay), dTdt(nlay))

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

      open(newunit=u,file='Y_400K_paper/Gao_2018_400_425.txt',action='read')
      read(u,*) ; read(u,*)
      nlines = 0
      do
        read(u,*,iostat=io)
        if (io /= 0) exit
        nlines = nlines + 1
      end do
      allocate(T_f(nlines),p_f(nlines))
      rewind(u)
      read(u,*); read(u,*)
      do i = 1, nlines
        read(u,*) T_f(i), p_f(i)
      end do
      p_f(:) = p_f(:)*1e5_dp
      close(u)
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

      Kzz(:) = 1e8_dp
      mu(:) = 2.33_dp
      grav = (10.0_dp**(4.25_dp))/100.0_dp
      met = 0.0_dp
      cp(:) = 1.3e4_dp
      dTdt(:) = 0.0_dp

      nd_atm(:) = (pl(:)*10.0_dp)/(kb*Tl(:))
      rho(:) = (pl(:)*10.0_dp*mu(:)*amu)/(kb * Tl(:))

      print*, 'i, pl [bar], T[k], Kzz [cm2 s-1]'
      do i = 1, nlay
        print*, i, pl(i)/1e5_dp, Tl(i), Kzz(i)
      end do

      allocate(VMR(nlay,2),sp_bg(2))
      sp_bg = (/'H2','He'/)
      VMR(:,1) = 0.85_dp
      VMR(:,2) = 0.15_dp

      allocate(sp(nsp), rho_d(nsp), mol_w_sp(nsp), mol_w_v(nsp))
      sp = (/'KCl', 'ZnS'/)
      rho_d = (/1.99_dp, 4.09_dp/)
      mol_w_sp = (/74.551_dp, 97.445_dp/)
      mol_w_v = (/74.551_dp, 65.38_dp/)

      m_seed = V_seed * rho_d(1)

      allocate(q_v(nlay,nsp), q_0(nlay), q_1(nlay,nsp), q_2(nlay), q0(2*nsp+2), q(nlay,2*nsp+2))
      allocate(r_c(nlay), m_c(nlay), vf(nlay,nsp+2), r_c_old(nlay), del(nlay))
      
      q_v(:,:) = 1e-30_dp
      q_0(:) = 1e-30_dp
      q_1(:,:) = 1e-30_dp
      q_2(:) = 1e-30_dp

      q0(:) = 1e-30_dp
      q0(1) = 1.17e-7_dp * mol_w_v(1)/mu(nlay)
      q0(2) = 3.63e-8_dp * mol_w_v(2)/mu(nlay)

      time = 0.0_dp
      n = 0
      r_c_old(:) = r_seed * 1e4_dp

      do n = 1, n_it

        !$omp parallel do default(shared), private(i, vf_tmp), schedule(dynamic)
        do i = 1, nlay
          call mini_cloud_vf(Tl(i), pl(i), grav, mu(i), VMR(i,:), rho_d(:), sp_bg, &
            & nsp, q_0(i), q_1(i,:), q_2(i), vf_tmp)
          vf(i,1) = vf_tmp(1)
          vf(i,2:nsp+1) = vf_tmp(2)
          vf(i,nsp+2) = vf_tmp(3)
        end do
        !$omp end parallel do

        q(:,1:nsp) = q_v(:,:)
        q(:,nsp+1) = q_0(:)
        q(:,nsp+2:2*nsp+1) = q_1(:,:)
        q(:,2*nsp+2) = q_2(:)

        call vert_adv_exp(nlay, nlev, t_step/2.0_dp, mu, grav, Tl, pl, pe, vf(:,:), nsp+2, q(:,nsp+1:))
        call vert_diff_imp(nlay, nlev, t_step/2.0_dp, mu, grav, Tl, pl, pe, Kzz, 2*nsp+2, q(:,:), q0(:))

        q_v(:,:) = q(:,1:nsp)
        q_0(:) = q(:,nsp+1)
        q_1(:,:) = q(:,nsp+2:2*nsp+1)
        q_2(:) = q(:,2*nsp+2)

        !$omp parallel do default(shared), private(i), schedule(dynamic)
        do i = 1, nlay
          call mini_cloud_3_lognormal_mix(i, Tl(i), pl(i), grav, mu(i), met, cp(i), VMR(i,:), t_step, sp, sp_bg, &
            & nsp, q_v(i,:), q_0(i), q_1(i,:), q_2(i), dTdt(i))
        end do
        !$omp end parallel do

        q(:,1:nsp) = q_v(:,:)
        q(:,nsp+1) = q_0(:)
        q(:,nsp+2:2*nsp+1) = q_1(:,:)
        q(:,2*nsp+2) = q_2(:)

        call vert_diff_imp(nlay, nlev, t_step/2.0_dp, mu, grav, Tl, pl, pe, Kzz, 2*nsp+2, q(:,:), q0(:))

        q_v(:,:) = q(:,1:nsp)
        q_0(:) = q(:,nsp+1)
        q_1(:,:) = q(:,nsp+2:2*nsp+1)
        q_2(:) = q(:,2*nsp+2)

        !$omp parallel do default(shared), private(i, vf_tmp), schedule(dynamic)
        do i = 1, nlay
          call mini_cloud_vf(Tl(i), pl(i), grav, mu(i), VMR(i,:), rho_d(:), sp_bg, &
            & nsp, q_0(i), q_1(i,:), q_2(i), vf_tmp)
          vf(i,1) = vf_tmp(1)
          vf(i,2:nsp+1) = vf_tmp(2)
          vf(i,nsp+2) = vf_tmp(3)
        end do
        !$omp end parallel do

        q(:,1:nsp) = q_v(:,:)
        q(:,nsp+1) = q_0(:)
        q(:,nsp+2:2*nsp+1) = q_1(:,:)
        q(:,2*nsp+2) = q_2(:)

        call vert_adv_exp(nlay, nlev, t_step/2.0_dp, mu, grav, Tl, pl, pe, vf(:,:), nsp+2, q(:,nsp+1:))

        q_v(:,:) = max(q(:,1:nsp), 1e-30_dp)
        q_0(:) = max(q(:,nsp+1), 1e-30_dp)
        q_1(:,:) = max(q(:,nsp+2:2*nsp+1), 1e-30_dp)
        q_2(:) = max(q(:,2*nsp+2), 1e-30_dp)

        !$omp parallel do default(shared), private(i), schedule(dynamic)
        do i = 1, nlay
          call opac_mie(nsp, sp, Tl(i), mu(i), pl(i), q_0(i), q_1(i,:), rho_d(:), n_wl, wl, k_ext(i,:), ssa(i,:), g(i,:))
        end do
        !$omp end parallel do

        do i = 1, nlay
          rho_c_tot = sum(q_1(i,:))*rho(i)
          Nc = q_0(i)*nd_atm(i)
          m_c(i) = max(rho_c_tot/max(Nc, 1.0e-300_dp), m_seed)
          if (rho_c_tot > 1e-30_dp) then
            rho_d_mean = rho_c_tot / sum((q_1(i,:) * rho(i)) / rho_d(:))
          else
            rho_d_mean = rho_d(1)
          end if
          r_c(i) = max(((3.0_dp*m_c(i))/(4.0_dp*pi*rho_d_mean))**(1.0_dp/3.0_dp),r_seed) * 1e4_dp
        end do

        end = .True.
        do i = 1, nlay
          del(i) = abs(r_c_old(i) - r_c(i))/r_c_old(i)
          if ((del(i) > 1e-4_dp) .or. (del(i)/t_step > 1e-4_dp)) end = .False.
        end do
        r_c_old(:) = r_c(:)

        time = time + t_step

        call progress_update(n, n_it, time, maxval(del(:)/t_step), width=100)

        if ((end .eqv. .True.) .and. (n > int(1e5))) then
          print*, 'exit: ', n, n_it, end
          exit
        end if

      end do

      do i = 1, nlay
        print*, i, pl(i)/1e5_dp, Tl(i), q_v(i,:), q_0(i), q_1(i,:), r_c(i), vf(i,1), q_0(i)*nd_atm(i)
      end do

      call output(n, time, nlay)

    case(3)

      !! L-T transition dwarf case: TiO2, Al2O3, Fe, Mg2SiO4

      nsp = 4

      nlay = 192
      nlev = nlay + 1

      n_wl = 11
      allocate(wl_e(n_wl+1), wl(n_wl), k_ext(nlay,n_wl), ssa(nlay,n_wl), g(nlay,n_wl))

      wl_e = (/0.260, 0.420, 0.610, 0.850, 1.320, 2.020,2.500,3.500,4.400,8.70,20.00,324.68 /)
      wl(:) = (wl_e(2:n_wl+1) + wl_e(1:n_wl))/ 2.0_dp

      allocate(Tl(nlay), pl(nlay), pe(nlev), mu(nlay), Kzz(nlay), nd_atm(nlay), rho(nlay), cp(nlay), dTdt(nlay))

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

      mu(:) = 2.33_dp
      met = 0.0_dp
      grav = (10.0_dp**(4.5_dp))/100.0_dp

      T_eff = 1300.0_dp
      k_ir = 1e-2_dp
      do i = 1, nlay
        tau = (k_ir * (pl(i) * 10.0_dp))/(grav*100.0_dp)
        Tl(i) = (3.0_dp/4.0_dp * T_eff**4 * (2.0_dp/3.0_dp + tau))**(1.0_dp/4.0_dp)
      end do

      nd_atm(:) = (pl(:)*10.0_dp)/(kb*Tl(:))
      rho(:) = (pl(:)*10.0_dp*mu(:)*amu)/(kb * Tl(:))

      Kzz(:) = 1e8_dp
      cp(:) = 1.3e4_dp
      dTdt(:) = 0.0_dp

      print*, 'i, pl [bar], T[k], Kzz [cm2 s-1]'
      do i = 1, nlay
        print*, i, pl(i)/1e5_dp, Tl(i), Kzz(i)
      end do

      allocate(VMR(nlay,2),sp_bg(2))
      sp_bg = (/'H2','He'/)
      VMR(:,1) = 0.85_dp
      VMR(:,2) = 0.15_dp

      allocate(sp(nsp), rho_d(nsp), mol_w_sp(nsp), mol_w_v(nsp))
      sp = (/'TiO2   ', 'Al2O3  ', 'Fe     ', 'Mg2SiO4'/)
      rho_d = (/4.23_dp, 3.986_dp, 7.874_dp, 3.21_dp/)
      mol_w_sp = (/79.866_dp, 101.961_dp, 55.845_dp, 140.693_dp/)
      mol_w_v = (/79.866_dp, 26.98153860_dp, 55.845_dp, 24.305_dp/)

      m_seed = V_seed * rho_d(1)

      allocate(q_v(nlay,nsp), q_0(nlay), q_1(nlay,nsp), q_2(nlay), q0(2*nsp+2), q(nlay,2*nsp+2))
      allocate(r_c(nlay), m_c(nlay), vf(nlay,nsp+2), r_c_old(nlay), del(nlay))

      q_v(:,:) = 1e-30_dp
      q_0(:) = 1e-30_dp
      q_1(:,:) = 1e-30_dp
      q_2(:) = 1e-30_dp

      q0(:) = 1e-30_dp
      q0(1) = 9.33e-8_dp * mol_w_v(1)/mu(nlay)
      q0(2) = 2.69e-6_dp * mol_w_v(2)/mu(nlay)
      q0(3) = 2.88e-5_dp * mol_w_v(3)/mu(nlay)
      q0(4) = 3.55e-5_dp * mol_w_v(4)/mu(nlay)

      time = 0.0_dp
      n = 0
      r_c_old(:) = r_seed * 1e4_dp

      do n = 1, n_it

        !$omp parallel do default(shared), private(i, vf_tmp), schedule(dynamic)
        do i = 1, nlay
          call mini_cloud_vf(Tl(i), pl(i), grav, mu(i), VMR(i,:), rho_d(:), sp_bg, &
            & nsp, q_0(i), q_1(i,:), q_2(i), vf_tmp)
          vf(i,1) = vf_tmp(1)
          vf(i,2:nsp+1) = vf_tmp(2)
          vf(i,nsp+2) = vf_tmp(3)
        end do
        !$omp end parallel do

        q(:,1:nsp) = q_v(:,:)
        q(:,nsp+1) = q_0(:)
        q(:,nsp+2:2*nsp+1) = q_1(:,:)
        q(:,2*nsp+2) = q_2(:)

        call vert_adv_exp(nlay, nlev, t_step/2.0_dp, mu, grav, Tl, pl, pe, vf(:,:), nsp+2, q(:,nsp+1:))
        call vert_diff_imp(nlay, nlev, t_step/2.0_dp, mu, grav, Tl, pl, pe, Kzz, 2*nsp+2, q(:,:), q0(:))

        q_v(:,:) = q(:,1:nsp)
        q_0(:) = q(:,nsp+1)
        q_1(:,:) = q(:,nsp+2:2*nsp+1)
        q_2(:) = q(:,2*nsp+2)

        !$omp parallel do default(shared), private(i), schedule(dynamic)
        do i = 1, nlay
          call mini_cloud_3_lognormal_mix(i, Tl(i), pl(i), grav, mu(i), met, cp(i), VMR(i,:), t_step, sp, sp_bg, &
            & nsp, q_v(i,:), q_0(i), q_1(i,:), q_2(i), dTdt(i))
        end do
        !$omp end parallel do

        q(:,1:nsp) = q_v(:,:)
        q(:,nsp+1) = q_0(:)
        q(:,nsp+2:2*nsp+1) = q_1(:,:)
        q(:,2*nsp+2) = q_2(:)

        call vert_diff_imp(nlay, nlev, t_step/2.0_dp, mu, grav, Tl, pl, pe, Kzz, 2*nsp+2, q(:,:), q0(:))

        q_v(:,:) = q(:,1:nsp)
        q_0(:) = q(:,nsp+1)
        q_1(:,:) = q(:,nsp+2:2*nsp+1)
        q_2(:) = q(:,2*nsp+2)

        !$omp parallel do default(shared), private(i, vf_tmp), schedule(dynamic)
        do i = 1, nlay
          call mini_cloud_vf(Tl(i), pl(i), grav, mu(i), VMR(i,:), rho_d(:), sp_bg, &
            & nsp, q_0(i), q_1(i,:), q_2(i), vf_tmp)
          vf(i,1) = vf_tmp(1)
          vf(i,2:nsp+1) = vf_tmp(2)
          vf(i,nsp+2) = vf_tmp(3)
        end do
        !$omp end parallel do

        q(:,1:nsp) = q_v(:,:)
        q(:,nsp+1) = q_0(:)
        q(:,nsp+2:2*nsp+1) = q_1(:,:)
        q(:,2*nsp+2) = q_2(:)

        call vert_adv_exp(nlay, nlev, t_step/2.0_dp, mu, grav, Tl, pl, pe, vf(:,:), nsp+2, q(:,nsp+1:))

        q_v(:,:) = max(q(:,1:nsp), 1e-30_dp)
        q_0(:) = max(q(:,nsp+1), 1e-30_dp)
        q_1(:,:) = max(q(:,nsp+2:2*nsp+1), 1e-30_dp)
        q_2(:) = max(q(:,2*nsp+2), 1e-30_dp)

        !$omp parallel do default(shared), private(i), schedule(dynamic)
        do i = 1, nlay
          call opac_mie(nsp, sp, Tl(i), mu(i), pl(i), q_0(i), q_1(i,:), rho_d(:), n_wl, wl, k_ext(i,:), ssa(i,:), g(i,:))
        end do
        !$omp end parallel do

        do i = 1, nlay
          rho_c_tot = sum(q_1(i,:))*rho(i)
          Nc = q_0(i)*nd_atm(i)
          m_c(i) = max(rho_c_tot/max(Nc, 1.0e-300_dp), m_seed)
          if (rho_c_tot > 1e-30_dp) then
            rho_d_mean = rho_c_tot / sum((q_1(i,:) * rho(i)) / rho_d(:))
          else
            rho_d_mean = rho_d(1)
          end if
          r_c(i) = max(((3.0_dp*m_c(i))/(4.0_dp*pi*rho_d_mean))**(1.0_dp/3.0_dp),r_seed) * 1e4_dp
        end do

        end = .True.
        do i = 1, nlay
          del(i) = abs(r_c_old(i) - r_c(i))/r_c_old(i)
          if ((del(i) > 1e-4_dp) .or. (del(i)/t_step > 1e-4_dp)) end = .False.
        end do
        r_c_old(:) = r_c(:)

        time = time + t_step

        call progress_update(n, n_it, time, maxval(del(:)/t_step), width=100)

        if ((end .eqv. .True.) .and. (n > int(1e5))) then
          print*, 'exit: ', n, n_it, end
          exit
        end if

      end do

      do i = 1, nlay
        print*, i, pl(i)/1e5_dp, Tl(i), q_v(i,:), q_0(i), q_1(i,:), r_c(i), vf(i,1), q_0(i)*nd_atm(i)
      end do

      call output(n, time, nlay)

    case default
      print*, 'Invalid test case example: ', example
      stop
    end select

  call progress_end()

contains

  subroutine output(t, time, nlay)
    implicit none
    integer, intent(in) :: t, nlay
    double precision, intent(in) :: time
    integer, save :: u1, u2, u3, u4
    logical, save :: first_call = .True.

    if (first_call .eqv. .True.) then
      open(newunit=u1,file='results_3_lognormal_mix/tracers.txt',action='readwrite')
      write(u1,*) nsp
      write(u1,*) rho_d(:)
      write(u1,*) mol_w_sp(:)
      open(newunit=u2,file='results_3_lognormal_mix/opac_k.txt',action='readwrite')
      write(u2,*) wl(:)
      open(newunit=u3,file='results_3_lognormal_mix/opac_a.txt',action='readwrite')
      write(u3,*) wl(:)
      open(newunit=u4,file='results_3_lognormal_mix/opac_g.txt',action='readwrite')
      write(u4,*) wl(:)
      first_call = .False.
    end if

    do i = 1, nlay
      write(u1,*) t, time, Tl(i), pl(i), grav, mu(i), VMR(i,:), q_v(i,:), q_0(i), q_1(i,:), q_2(i), vf(i,:), dTdt(i)
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
    real(dp), intent(in) :: var
    integer :: jl, jm, ju

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

end program test_mini_cloud_3
