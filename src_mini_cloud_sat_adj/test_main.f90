program test_mini_cloud_sat_adj
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  use mini_cloud_sat_adj_mod, only : mini_cloud_sat_adj
  use mini_cloud_vf_sat_adj_mod, only : mini_cloud_vf_sat_adj
  use mini_cloud_opac_mie_sat_adj_mod, only : opac_mie_sat_adj
  use vert_diff_exp_mod, only : vert_diff_exp
  use vert_diff_imp_mod, only : vert_diff_imp
  use vert_adv_exp_mod, only : vert_adv_exp
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: kb = 1.380649e-16_dp ! erg K^-1 - Boltzmann's constant
  real(dp), parameter :: amu = 1.66053906660e-24_dp ! g - Atomic mass unit
  real(dp), parameter :: R = 8.31446261815324e7_dp

  integer :: example, n_it
  real(dp) :: t_step, time

  integer :: nlay, nlev, i, n, u
  character(len=20) :: sp
  character(len=20), allocatable, dimension(:) :: sp_bg
  real(dp), allocatable, dimension(:) :: Tl, pl, mu, pe, nd_atm, rho, Kzz
  real(dp), allocatable, dimension(:) :: q_1, q_v, q0, q_1_old, del
  real(dp), allocatable, dimension(:,:) :: VMR, q, vf
  real(dp) :: grav

  integer :: n_wl
  real(dp), allocatable, dimension(:) :: wl_e, wl
  real(dp), allocatable, dimension(:,:) :: k_ext, ssa, g

  integer :: nlines, io, idx, idx1
  real(dp) :: p_bot, p_top
  real(dp), allocatable, dimension(:) :: T_f, p_f

  real(dp) :: tau_cond, sigma, met, mol_w_sp, r_m, rho_d

  logical :: end

  !! time step
  t_step = 500.0_dp

  !! Number of iterations
  n_it = 100000

  !! Start time
  time = 6840.0_dp

  !! example select
  example = 2

  select case (example)

    case(1)

    case(2)

      !! In this example we perform a 1D test of mini_cloud with diffusion and settling included

      nlay = 192
      nlev = nlay + 1

      n_wl = 11
      allocate(wl_e(n_wl+1), wl(n_wl), k_ext(nlay,n_wl), ssa(nlay,n_wl), g(nlay,n_wl))

      !! Wavelengths to calculate opacity
      wl_e = (/0.260, 0.420, 0.610, 0.850, 1.320, 2.020,2.500,3.500,4.400,8.70,20.00,324.68 /)
      wl(:) = (wl_e(2:n_wl+1) +  wl_e(1:n_wl))/ 2.0_dp

      !! Allocate all variables, set constant values
      allocate(Tl(nlay), pl(nlay), pe(nlev), mu(nlay), Kzz(nlay), nd_atm(nlay), rho(nlay))

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
      open(newunit=u,file='Y_400K_paper/Gao_2018_400_325.txt',action='read')
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
      grav = (10.0_dp**(3.25_dp))/100.0_dp

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
      sp = 'KCl'
      rho_d = 1.99_dp
      mol_w_sp = 74.5513_dp

      met = 0.0_dp

      r_m = 10.0_dp * 1e-4_dp
      sigma = 2.0_dp
      tau_cond = 10.0_dp

      allocate(q_v(nlay), q_1(nlay), q0(2), q(nlay,2))
      allocate(vf(nlay,1), q_1_old(nlay), del(nlay))

      q_v(:) = 1e-30_dp
      q_1(:) = 1e-30_dp

      q0(1) = 1.17e-7_dp * mol_w_sp/mu(nlay)
      q0(2) = 1e-30_dp

      time = 0.0_dp
      n = 0

      do n = 1, n_it

        !$omp parallel do default(shared), private(i), schedule(dynamic)
        do i = 1, nlay
          !! Calculate settling velocity for this layer
          call mini_cloud_vf_sat_adj(Tl(i), pl(i), grav, mu(i), VMR(i,:), rho_d, sp_bg, r_m, sigma, vf(i,1))
        end do
        !$omp end parallel do  

        q(:,1) = q_v(:)
        q(:,2) = q_1(:)

        call vert_adv_exp(nlay, nlev, t_step/2.0_dp, mu, grav, Tl, pl, pe, vf(:,1), 1, q(:,2))

        call vert_diff_imp(nlay, nlev, t_step/2.0_dp, mu, grav, Tl, pl, pe, Kzz(:), 2, q(:,:), q0(:))

        q_v(:) = q(:,1)
        q_1(:) = q(:,2)

        !$omp parallel do default(shared), private(i), schedule(dynamic)
        do i = 1, nlay
          !! Call mini-cloud and perform integrations for a single layer
          call mini_cloud_sat_adj(i, t_step, mol_w_sp, sp, rho_d, Tl(i), pl(i), rho(i), met, tau_cond, q_v(i), q_1(i))
        end do
        !$omp end parallel do  

        q(:,1) = q_v(:)
        q(:,2) = q_1(:)

        call vert_diff_imp(nlay, nlev, t_step/2.0_dp, mu, grav, Tl, pl, pe, Kzz(:), 2, q(:,:), q0(:))

        q_v(:) = q(:,1)
        q_1(:) = q(:,2)

        !$omp parallel do default(shared), private(i), schedule(dynamic)
        do i = 1, nlay
          !! Calculate settling velocity for this layer
          call mini_cloud_vf_sat_adj(Tl(i), pl(i), grav, mu(i), VMR(i,:), rho_d, sp_bg, r_m, sigma, vf(i,1))
        end do
        !$omp end parallel do  

        q(:,1) = q_v(:)
        q(:,2) = q_1(:)
        
        call vert_adv_exp(nlay, nlev, t_step/2.0_dp, mu, grav, Tl, pl, pe, vf(:,1), 1, q(:,2))

        q_v(:) = q(:,1)
        q_1(:) = q(:,2)

        !$omp parallel do default(shared), private(i), schedule(dynamic)
        do i = 1, nlay
          !! Calculate the opacity at the wavelength grid
          call opac_mie_sat_adj(1, sp, Tl(i), mu(i), pl(i), q_1(i), r_m, rho_d, sigma, n_wl, wl, k_ext(i,:), ssa(i,:), g(i,:))
        end do
        !$omp end parallel do



        end = .True.

        do i = 1, nlay
          del(i) = abs(q_1_old(i) - q_1(i))/q_1_old(i)
          if ((del(i) > 1e-4_dp) .or.  (del(i)/t_step > 1e-4_dp)) then
            end = .False.
            exit
          end if
        end do
        q_1_old(:) = q_1(:)

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
        print*, i, pl(i)/1e5_dp, Tl(i), Kzz(i), q_v(i), q_1(i), r_m, vf(i,1)
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
    integer, save :: u1, u2
    logical, save :: first_call = .True.

    if (first_call .eqv. .True.) then
      open(newunit=u1,file='results_sat_adj/tracers.txt',action='readwrite')
      open(newunit=u2,file='results_sat_adj/opac.txt',action='readwrite')
      write(u2,*) wl(:)
      first_call = .False.
    end if

    do i = 1, nlay
      write(u1,*) t, time, Tl(i), pl(i), grav, mu(i), VMR(i,:), q_v(i), q_1(i), r_m, vf(i,:)
      write(u2,*) t, time, k_ext(i,:), ssa(i,:), g(i,:)
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

end program test_mini_cloud_sat_adj
