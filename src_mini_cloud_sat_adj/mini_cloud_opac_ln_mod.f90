module mini_cloud_opac_ln_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none


  integer, parameter :: dp = REAL64

  real(dp), parameter :: pi = 4.0_dp*atan(1.0_dp) ! value of pi
  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: amu = 1.66053906660e-24_dp ! g mol-1 (note, has to be cgs g mol-1 !!!)


  logical :: first_call = .True.

  real(dp), allocatable, dimension(:) :: k_ext_f , alb_f, gg_f

  public :: ln_dist_opac
  private :: read_ln_table

  contains

  subroutine ln_dist_opac(n_dust, ln_file, T_in, mu_in, P_in, q_c, & 
    & rm, rho_d, sigma, n_wl, k_ext, alb, gg)
    implicit none

    integer, intent(in) :: n_dust, n_wl
    character(len=100), intent(in) :: ln_file
    real(dp), intent(in) :: T_in, mu_in, P_in
    real(dp), intent(in) :: q_c, rm, sigma, rho_d

    real(dp), dimension(n_wl), intent(out) :: k_ext, alb, gg

    real(dp) :: rho, n_d

    if (first_call .eqv. .True.) then
      call read_ln_table(ln_file,n_wl)
      first_call = .False.
    end if

    if (q_c < 1e-10_dp) then
      k_ext(:) = 0.0_dp
      alb(:) = 0.0
      gg(:) = 0.0_dp
      return
    end if

    rho = (P_in * 10.0_dp * mu_in * amu)/(kb * T_in)

    n_d = (3.0_dp * q_c * rho)/(4.0_dp*pi*rho_d*rm**3) * exp(-9.0_dp/2.0_dp * log(sigma)**2)! Find particle number density

    k_ext(:) = (k_ext_f(:) * n_d)/rho
    alb(:) = alb_f(:)
    gg(:) = gg_f(:)


  end subroutine ln_dist_opac

  subroutine read_ln_table(ln_file, n_wl)
    implicit none

    character(len=100), intent(in) :: ln_file
    integer, intent(in) :: n_wl

    integer :: u, i
    real(dp) :: wl_dum


    allocate(k_ext_f(n_wl), alb_f(n_wl), gg_f(n_wl))


    open(newunit=u,file=trim(ln_file),action='read')

    read(u,*)

    do i = 1, n_wl
      read(u,*) wl_dum, k_ext_f(i), alb_f(i), gg_f(i)
    end do

    close(u)

  end subroutine read_ln_table


end module mini_cloud_opac_ln_mod