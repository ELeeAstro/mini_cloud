program test_mini_cloud_1_simple
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  use mini_cloud_1_s_mod, only : mini_cloud_1_s, rho_c, tau_chem, rm, sigma, mol_w_sp, Kzz_deep, q_v_deep
  implicit none

  integer, parameter :: dp = REAL64


  integer :: example, tt, n_it
  character(len=20) :: sp, sp_bg(3)
  logical :: deep_flag
  real(dp) :: T_in, P_in, VMR_in(3), mu_in, grav_in
  real(dp) :: q_c, q_v, v_f, k_ext, ssa, g
  real(dp) :: t_step, time


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
      call mini_cloud_1_s(deep_flag, T_in, P_in, grav_in, mu_in, VMR_in, t_step, sp, sp_bg, q_v, q_c, v_f)

      !! Call the ADT approximate opacity routine for a layer at a single wavelength
      !call adt_1_s
      k_ext = 0.0_dp
      ssa = 0.0_dp
      g = 0.0_dp

      !! increment time
      time = time + t_step

      !! Print to screen current progress
      print*, tt, time, P_in * 1e-5_dp, 'bar ', T_in, 'K ', mu_in, 'g mol-1 ',  trim(sp)

      print*, 'q', tt, q_v, q_c, v_f
      print*, 'r', tt, rm * 1e4_dp, sigma
      print*, 'o', tt, k_ext, ssa, g

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
      first_call = .False.
    end if

    write(u1,*) t, time, T_in, P_in, grav_in, mu_in, VMR_in(:), q_v, q_c, v_f
    write(u2,*) t, time, k_ext, ssa, g

  end subroutine output

end program test_mini_cloud_1_simple