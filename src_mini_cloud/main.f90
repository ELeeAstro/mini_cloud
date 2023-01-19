program main
  use mini_cloud_precision
  use mini_cloud_class
  use mini_cloud_i_dvode
  use mini_cloud_vf
  implicit none

  integer, parameter :: n_dust = 4
  real(dp) :: T_in, P_in, mu_in, t_step, grav

  character(len=11), dimension(n_dust) :: sp

  real(dp), dimension(4) :: k
  real(dp), dimension(n_dust) :: k3
  real(dp), dimension(n_dust) :: VMR
  real(dp) :: vf

  integer :: tt, n_it, example
  real(dp) :: time

  sp(1) = 'TiO2  '
  sp(2) = 'Al2O3 '
  sp(3) = 'Fe    '
  sp(4) = 'Mg2SiO4'

  k(:) = 1.0e-30_dp
  k3(:) = 1.0e-30_dp
  VMR(1) = 1.041e-7_dp
  VMR(2) = 2.771e-6_dp
  VMR(3) = 2.231e-5_dp
  VMR(4) = 3.588e-5_dp

  ! Number of iterations and start time
  n_it = 3000!400
  time = 6840.0_dp !0.0e0_dp

  T_in = 1500.0_dp
  P_in = 1.0e4_dp
  mu_in = 2.33_dp
  grav = 2000.0_dp

  t_step = 10.0_dp

  example = 1

  ! Start time iteration
  do tt = 1, n_it

    select case (example)
    case(1)

      !! In this example, we timestep a call to mini-cloud while slowly increasing the temperature

      ! Start sinusoid temperature variation
      T_in = 1500.0_dp + 1000.0_dp * sin(2.0_dp * 3.14_dp * 0.5_dp *  time)

      ! DIHRT test output
      call output(tt, time)

      ! Call DIHRT and perform integrations
      call mini_cloud_dvode(n_dust, T_in, P_in, t_step, sp(:), k(:), k3(:), VMR(:))
      !call mini_cloud_dvode(n_dust, T_in, P_in, t_step, sp(:), k(:), k3(:), VMR(:))
      !call mini_cloud_dvode(n_dust, T_in, P_in, t_step, sp(:), k(:), k3(:), VMR(:))



      print*, k(:), k(2)/k(1) * 1e4
      print*, k3(:), VMR(:)

      ! Call the vertical falling rate routine
      call calc_vf(n_dust, T_in, P_in, mu_in, grav, k(:), k3(:), vf)

      ! increment time
      time = time + t_step

      ! Start sinusoid temperature variation
      T_in = 1500.0_dp + 1000.0_dp * sin(2.0_dp * 3.14_dp * 0.5_dp *  time)

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
    integer, save :: u1
    logical, save :: first_call = .True.

    if (first_call .eqv. .True.) then
      open(newunit=u1,file='tracers.txt',action='readwrite')
      first_call = .False.
    end if

    write(u1,*) t, time, T_in, P_in, k(:), k3(:), VMR(:), vf 

  end subroutine output

end program main
