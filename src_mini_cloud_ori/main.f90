program main
  use mini_cloud_precision
  use mini_cloud_class
  use mini_cloud_i_dlsode
  use mini_cloud_vf
  use mini_cloud_opac_mie_mod, only : opac_mie
  implicit none

  integer, parameter :: n_dust = 4
  real(dp) :: T_in, P_in, mu_in, t_step, grav

  character(len=11), dimension(n_dust) :: sp

  real(dp), dimension(4) :: k
  real(dp), dimension(n_dust) :: k3
  real(dp), dimension(n_dust) :: VMR, VMR0
  real(dp) :: vf

  integer, parameter :: n_wl = 11
  real(dp), dimension(n_wl+1) :: wl_e
  real(dp), dimension(n_wl) :: wl, k_ext, alb, gg

  integer :: tt, n_it, example
  real(dp) :: time

  !! Chose name of individual species
  sp(1) = 'TiO2   '
  sp(2) = 'Mg2SiO4'
  sp(3) = 'Fe     '
  sp(4) = 'Al2O3  '

  !! Initial moments and species VMR
  k(:) = 1.0e-30_dp
  k3(:) = 1.0e-30_dp
  VMR(1) = 1.041e-7_dp
  VMR(2) = 3.588e-5_dp
  VMR(3) = 2.231e-5_dp
  VMR(4) = 2.771e-6_dp

  VMR0(:) = VMR(:)

  !! Initialisation for v=arrays
  vf = 1e-30_dp
  k_ext(:) = 0.0_dp
  alb(:) = 0.0_dp
  gg(:) = 0.0_dp

  !! Wavelengths to calculate opacity
  wl_e = (/0.260, 0.420, 0.610, 0.850, 1.320, 2.020,2.500,3.500,4.400,8.70,20.00,324.68 /)
  wl(:) = (wl_e(2:n_wl+1) +  wl_e(1:n_wl))/ 2.0_dp

  !! input atmosphere properties
  T_in = 1500.0_dp
  P_in = 1.0e4_dp
  mu_in = 2.33_dp
  grav = 1000.0_dp

  !! time step
  t_step = 100.0_dp

  !! Number of iterations
  n_it = 10000

  !! Start time
  time = 6840.0_dp

  !! Select example
  example = 1

  ! Start time iteration
  do tt = 1, n_it

    select case (example)
    case(1)

      !! In this example, we timestep a call to mini-cloud while slowly increasing the temperature

      !! Start sinusoid temperature variation [K]
      T_in = 1600.0_dp + 1000.0_dp * sin(2.0_dp * 3.14_dp * 0.01_dp *  time)

      ! DIHRT test output
      call output(tt, time)

      ! Call DIHRT and perform integrations
      call mini_cloud_dlsode(n_dust, T_in, P_in, t_step, sp(:), k(:), k3(:), VMR(:), VMR0(:))

      ! Call the vertical falling rate routine
      call calc_vf(n_dust, T_in, P_in, mu_in, grav, k(:), k3(:), vf)

      !! Call the opaicity routine
      call opac_mie(n_dust, sp(:), T_in, mu_in, P_in, k(:), k3(:), n_wl, wl, k_ext, alb, gg)


      print*, 'k', tt, k(:), k(2)/k(1) * 1e4
      print*, 'k3', tt, k3(:), VMR(:)
      print*, 'o', tt, k(2)/k(1)*1e4_dp,k_ext(7), alb(7), gg(7)

      ! increment time
      time = time + t_step

      ! Print to screen current progress
      print*, tt, time, P_in * 1e-5_dp, 'bar', T_in, 'K'
    case(2)
      !! In this example, we timestep a call to mini-cloud
      !! and use a sawtooth pressure wave in log space for 4 oscilations between
      !! PG_min and PG_max - this can be quite 'shocking' and numerically difficult

    case default 
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
      open(newunit=u1,file='results_ori/tracers.txt',action='readwrite')
      open(newunit=u2,file='results_ori/opac.txt',action='readwrite')
      write(u2,*) wl(:)
      first_call = .False.
    end if

    write(u1,*) t, time, T_in, P_in, k(:), k3(:), VMR(:), vf 
    write(u2,*) t, time, k_ext(:), alb(:), gg(:)

  end subroutine output

end program main
