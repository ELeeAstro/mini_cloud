program main
  use mini_cloud_precision
  use mini_cloud_class
  use mini_cloud_i_dvode
  implicit none

  integer, parameter :: n_dust = 3
  real(dp) :: T_in, P_in, t_step

  character(len=11), dimension(n_dust) :: sp

  real(dp), dimension(4) :: k
  real(dp), dimension(n_dust) :: k3
  real(dp), dimension(n_dust) :: VMR

  integer :: tt, n_it 
  real(dp) :: time

  sp(1) = 'TiO2'
  sp(2) = 'Fe'
  sp(3) = 'Mg2SiO4'

  k(:) = 1.0e-30_dp
  k3(:) = 1.0e-30_dp
  VMR(1) = 1e-7_dp
  VMR(2) = 2.2e-5_dp
  VMR(3) = 3.6e-5_dp

  ! Number of iterations and start time
  n_it = 1500
  time = 6840.0_dp !0.0e0_dp

  T_in = 1500.0_dp
  P_in = 1.0e4_dp

  t_step = 10.0_dp

  ! Start time iteration
  do tt = 1, n_it

    ! Start sinusoid temperature variation
    T_in = 1500.0_dp + 1000.0_dp * sin(2.0_dp * 3.14_dp * 0.5_dp *  time)

    ! DIHRT test output
    call output(tt, time)

    ! Call DIHRT and perform integrations
    call mini_cloud_dvode(n_dust, T_in, P_in, t_step, sp(:), k(:), k3(:), VMR(:))


    print*, k(:), k(2)/k(1) * 1e4
    print*, k3(:), VMR(:)

    ! increment time
    time = time + t_step

    ! Start sinusoid temperature variation
    T_in = 1500.0_dp + 1000.0_dp * sin(2.0_dp * 3.14_dp * 0.5_dp *  time)

    ! Print to screen current progress
    print*, tt, time, P_in * 1e-5_dp, 'bar', T_in, 'K'

  end do

contains

  subroutine output(t, time)
    implicit none

    integer, intent(in) :: t
    double precision, intent(in) :: time
    integer, save :: u1, u2, u3, u4, u5
    logical, save :: first_call = .True.

    if (first_call .eqv. .True.) then
      open(newunit=u1,file='tracers.txt',action='readwrite')

      first_call = .False.
    end if

    write(u1,*) t, time, T_in, P_in, k(:), k3(:), VMR(:) 

  end subroutine output

end program main
