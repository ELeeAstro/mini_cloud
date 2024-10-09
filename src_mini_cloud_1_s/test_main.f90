program test_mini_cloud_1_simple
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  use mini_cloud_1_s_mod, only : mini_cloud_1_s, rho_c, tau_chem, rm, sigma, mol_w_sp, Kzz_deep, q_v_deep
  use mini_cloud_vf_mod, only : mini_cloud_vf
  use mini_cloud_opac_mie_mod, only : opac_mie
  implicit none

  integer, parameter :: dp = REAL64


  integer :: example, tt, n_it
  character(len=20) :: sp, sp_bg(3)
  logical :: deep_flag
  real(dp) :: T_in, P_in, VMR_in(3), mu_in, grav_in
  real(dp) :: q_c, q_v, v_f
  real(dp) :: t_step, time

  integer :: n_wl
  real(dp), allocatable, dimension(:) :: wl_e, wl, k_ext, ssa, g

  !! time step
  t_step = 100.0_dp

  !! Number of iterations
  n_it = 10000

  !! Start time
  time = 6840.0_dp

  !! example select
  example = 1

  !! Initial conditions
  q_v = 1.17e-7_dp  ! ~Mg abundance ratio at Solar (VMR)
  q_c = 1e-30_dp    ! ~Zero clouds at start 


  n_wl = 11
  allocate(wl_e(n_wl+1), wl(n_wl), k_ext(n_wl), ssa(n_wl), g(n_wl))

  !! Wavelengths to calculate opacity
  wl_e = (/0.260, 0.420, 0.610, 0.850, 1.320, 2.020,2.500,3.500,4.400,8.70,20.00,324.68 /)
  wl(:) = (wl_e(2:n_wl+1) +  wl_e(1:n_wl))/ 2.0_dp

 ! Start time iteration
  do tt = 1, n_it

    select case (example)

    case(1)

      !! In this example, we timestep a call to mini-cloud while slowly increasing the temperature

      !! Start sinusoid temperature variation [K]
      T_in = 1600.0_dp + 1000.0_dp * sin(2.0_dp * 3.14_dp * 0.01_dp *  time)

      !! Assume constant pressure [pa]
      P_in = 1e5_dp

       !! Assume constant H2, He and H background VMR @ approx solar
      sp_bg = (/'H2','He','H '/)
      VMR_in(1) = 0.85_dp
      VMR_in(2) = 0.15_dp
      VMR_in(3) = 1e-6_dp

      !! Assume constant background gas mean molecular weight [g mol-1] @ approx solar
      mu_in = 2.33_dp

      !! Assume constant gravity [m s-2]
      grav_in = 10.0_dp

      !! Assumed condensate species
      sp = 'KCl'

      !! Variables to be sent to mini-cloud module via global variables !!
      !! Inside the loop to demonstrate they can be changed with time if required !!
      !! In this specific code everything is in cgs
      !! --------------------------- !! 
      rho_c = 1.99_dp
      tau_chem = 10.0_dp

      rm = 1.0_dp * 1e-4_dp
      sigma = 2.0_dp

      mol_w_sp = 74.551_dp

      Kzz_deep = 1e8_dp
      q_v_deep = 1.17e-7_dp * mol_w_sp/mu_in !! MMR for deep mixing ratio
      
      !! --------------------------- !! 

      !! Convert vapour and condensate VMR to MMR at first iteration
      if (tt == 1) then
        q_v = q_v * mol_w_sp/mu_in
        q_c = q_c * mol_w_sp/mu_in
      end if

      !! Are we in the deep layer?
      deep_flag = .False.

      !! mini-cloud test output
      call output(tt, time)

      !! Call mini-cloud and perform integrations for a single layer
      call mini_cloud_1_s(deep_flag, T_in, P_in, grav_in, mu_in, t_step, sp, q_v, q_c)

      !! Calculate settling velocity 
      call mini_cloud_vf(T_in, P_in, grav_in, mu_in, VMR_in, rho_c, sp_bg, rm, sigma, v_f)

      !! Calculate the opacity at the weavelength grid
      call opac_mie(1, sp, T_in, mu_in, P_in, q_c, rm, rho_c, sigma, n_wl, wl, k_ext, ssa, g)

      !! increment time
      time = time + t_step

      !! Print to screen current progress
      print*, tt, time, P_in * 1e-5_dp, 'bar ', T_in, 'K ', mu_in, 'g mol-1 ',  trim(sp)

      print*, 'q', tt, q_v, q_c, v_f
      print*, 'r', tt, rm * 1e4_dp, sigma
      print*, 'o', tt, k_ext(1), ssa(1), g(1), k_ext(n_wl), ssa(n_wl), g(n_wl)

    case(2)
      !! In this example, we timestep a call to mini-cloud
      !! and use a sawtooth pressure wave in log space for 4 oscilations between
      !! PG_min and PG_max - this can be quite 'shocking' and numerically difficult

    case default 
      print*, 'Invalud test case example: ', example
      stop
    end select

  end do

contains

  subroutine output(t, time)
    implicit none

    integer, intent(in) :: t
    double precision, intent(in) :: time
    integer, save :: u1, u2
    logical, save :: first_call = .True.

    if (first_call .eqv. .True.) then
      open(newunit=u1,file='results_1_s/tracers.txt',action='readwrite')
      open(newunit=u2,file='results_1_s/opac.txt',action='readwrite')
      write(u2,*) wl(:)
      first_call = .False.
    end if

    write(u1,*) t, time, T_in, P_in, grav_in, mu_in, VMR_in(:), q_v, q_c, v_f
    write(u2,*) t, time, k_ext(:), ssa(:), g(:)

  end subroutine output

end program test_mini_cloud_1_simple