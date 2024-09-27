program test_mini_cloud_2_mix
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  use mini_cloud_2_mix_mod, only : mini_cloud_2_mix, rho_d, mol_w_sp
  use mini_cloud_vf_mod, only : mini_cloud_vf
  use mini_cloud_opac_mie_mod, only : opac_mie
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: kb = 1.380649e-23_dp
  real(dp), parameter :: amu = 1.66053906892e-27_dp

  integer, parameter ::  n_bg = 3, n_dust = 4
  integer :: example, tt, n_it,

  real(dp) :: T_in, P_in, mu_in, grav_in, nd_atm, rho
  real(dp) :: q_0, v_f, r_c, m_c
  real(dp) :: t_step, time

  !! Background
  character(len=20), dimension(n_bg) :: sp_bg, VMR_bg

  !! Dust details
  character(len=20), dimension(n_dust) :: sp, q_1s, q_v
  real(dp), dimension(n_dust) :: rho_s, mw_s

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
  q_v(:) = (/1.041e-7_dp, 2.771e-6_dp, 2.231e-5_dp, 3.588e-5_dp/)  !Ti, Mg, Fe and Al abundances at start
  q_0 = 1.0e-30_dp    ! ~Zero clouds at start 
  q_1s(:) = 1.0e-30_dp


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
      sp_bg(:) = (/'H2','He','H '/)
      VMR_bg(:) = (/0.85_dp,0.15_dp,1e-6/)
  
      !! Assume constant background gas mean molecular weight [g mol-1] @ approx solar
      mu_in = 2.33_dp

      !! Assume constant gravity [m s-2]
      grav_in = 10.0_dp

      !! Assumed condensate species
      sp(:) = (/'TiO2', 'Al2O3', 'Fe', 'MgSiO3'/)

      !! Number density [m-3] of layer
      nd_atm = P_in/(kb*T_in)  

      !! Mass density of layer
      rho = (P_in*mu_in*amu)/(kb * T_in) ! Mass density [kg m-3]

      !! Call mini-cloud and perform integrations for a single layer
      call mini_cloud_2_mix(n_dust, T_in, P_in, grav_in, mu_in, VMR_bg, t_step, sp, sp_bg, q_v, q_0, q_1s)

      !! Calculate settling velocity for this layer
      call mini_cloud_vf(n_dust, T_in, P_in, grav_in, mu_in, VMR_bg, rho_d, sp_bg, q_0, q_1s, v_f)

      !! Calculate the opacity at the weavelength grid
      call opac_mie(n_dust, sp, T_in, mu_in, P_in, q_0, q_1s, rho_d, n_wl, wl, k_ext, ssa, g)

      !! increment time
      time = time + t_step

      !! Print to screen current progress
      print*, tt, time, P_in * 1e-5_dp, 'bar ', T_in, 'K ', mu_in, 'g mol-1 ',  trim(sp)

      !! Mean mass of particle
      m_c = (sum(q_1s*rho))/(q_0*nd_atm)

      !! Mass weighted mean radius of particle
      r_c = ((3.0_dp*m_c)/(4.0_dp*pi*rho_d*1000.0_dp))**(1.0_dp/3.0_dp) * 1e6_dp


      print*, 'q', tt, q_v, q_0, q_1, v_f
      print*, 'r', tt, m_c, r_c
      print*, 'o', tt, k_ext(1), ssa(1), g(1), k_ext(n_wl), ssa(n_wl), g(n_wl)

      !! mini-cloud test output
      call output(tt, time)

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
      open(newunit=u1,file='results_2_mix/tracers.txt',action='readwrite')
      open(newunit=u2,file='results_2_mix/opac.txt',action='readwrite')
      write(u2,*) wl(:)
      first_call = .False.
    end if

    write(u1,*) t, time, T_in, P_in, grav_in, mu_in, VMR_in(:), q_v, q_0, q_1, v_f
    write(u2,*) t, time, k_ext(:), ssa(:), g(:)

  end subroutine output

end program test_mini_cloud_2_mix
