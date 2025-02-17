program test_mini_cloud_2
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  use mini_cloud_2_mod, only : mini_cloud_2, rho_d, mol_w_sp
  use mini_cloud_vf_mod, only : mini_cloud_vf
  use mini_cloud_opac_mie_mod, only : opac_mie
  use vert_diff_exp_mod, only : vert_diff_exp
  use vert_diff_imp_mod, only : vert_diff_imp
  use vert_adv_exp_McCormack_mod, only : vert_adv_exp_McCormack
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: kb = 1.380649e-16_dp ! erg K^-1 - Boltzmann's constant
  real(dp), parameter :: amu = 1.66053906660e-24_dp ! g - Atomic mass unit
  real(dp), parameter :: r_seed = 1e-7_dp

  integer :: example, n_it
  real(dp) :: t_step, time

  integer :: nlay, nlev, i, n, u
  character(len=20) :: sp
  character(len=20), allocatable, dimension(:) :: sp_bg
  real(dp), allocatable, dimension(:) :: Tl, pl, mu, Kzz, pe, nd_atm, rho
  real(dp), allocatable, dimension(:) :: q_0, q_1, q_v, vf, r_c, m_c, q0, r_c_old, del
  real(dp), allocatable, dimension(:,:) :: VMR, q
  real(dp) :: grav

  integer :: n_wl
  real(dp), allocatable, dimension(:) :: wl_e, wl
  real(dp), allocatable, dimension(:,:) :: k_ext, ssa, g

  integer :: nlines, io, idx, idx1
  real(dp) :: p_bot, p_top
  real(dp), allocatable, dimension(:) :: T_f, p_f, Kzz_f
  real(dp) :: V_seed, m_seed

  logical :: end

  !! time step
  t_step = 1000.0_dp

  !! Number of iterations
  n_it = 1000000

  !! Start time
  time = 6840.0_dp

  !! example select
  example = 2

  select case (example)

    case(1)

      !! In this example, we timestep a call to mini-cloud while slowly increasing the temperature
      nlay = 1

      n_wl = 11
      allocate(wl_e(n_wl+1), wl(n_wl), k_ext(nlay,n_wl), ssa(nlay,n_wl), g(nlay,n_wl))

      !! Wavelengths to calculate opacity
      wl_e = (/0.260, 0.420, 0.610, 0.850, 1.320, 2.020,2.500,3.500,4.400,8.70,20.00,324.68 /)
      wl(:) = (wl_e(2:n_wl+1) +  wl_e(1:n_wl))/ 2.0_dp

      !! Allocate all variables, set constant values
      allocate(Tl(nlay), pl(nlay), mu(nlay), nd_atm(nlay), rho(nlay))

      !! Pressure in pa
      pl(1) = 1e5_dp

      !! Assume constant H2, He and H background VMR @ approx solar
      allocate(VMR(nlay,3),sp_bg(3))
      sp_bg = (/'H2','He','H '/)
      VMR(1,1) = 0.85_dp
      VMR(1,2) = 0.15_dp
      VMR(1,3) = 1e-6_dp

      !! Assume constant background gas mean molecular weight [g mol-1] @ approx solar
      mu(1) = 2.33_dp

      !! Assume constant gravity [m s-2]
      grav = 10.0_dp

      !! Assumed condensate species
      sp = 'KCl'
      rho_d = 1.99_dp
      mol_w_sp = 74.551_dp

      allocate(q_v(nlay),q_0(nlay),q_1(nlay))
      allocate(r_c(nlay), m_c(nlay), vf(nlay))

      !! Initial conditions
      q_v(1) = 1.17e-7_dp  ! ~K abundance ratio at Solar (VMR)
      q_0(1) = 1.0e-30_dp    ! ~Zero clouds at start 
      q_1(1) = 1.0e-30_dp


      !! Change vapur VMR to mass density ratio
      q_v(1) = q_v(1) * mol_w_sp/mu(1)

      V_seed = 4.0_dp/3.0_dp * pi * r_seed**3
      m_seed = V_seed * rho_d

      do n = 1, n_it

        !! Start sinusoid temperature variation [K]
        Tl(1) = 1600.0_dp + 1000.0_dp * sin(2.0_dp * 3.14_dp * 0.01_dp *  time)

        !! Number density [cm-3] of layer
        nd_atm(1) = (pl(1)*10.0_dp)/(kb*Tl(1))  

        !! Mass density of layer
        rho(1) = (pl(1)*10.0_dp*mu(1)*amu)/(kb * Tl(1)) ! Mass density [g cm-3]

        !! Call mini-cloud and perform integrations for a single layer
        call mini_cloud_2(Tl(1), pl(1), grav, mu(1), VMR(1,:), t_step, sp, sp_bg, q_v(1), q_0(1), q_1(1))

        !! Calculate settling velocity for this layer
        call mini_cloud_vf(Tl(1), pl(1), grav, mu(1), VMR(1,:), rho_d, sp_bg, q_0(1), q_1(1), vf(1))

        !! Calculate the opacity at the weavelength grid
        call opac_mie(1, sp, Tl(1), mu(1), pl(1), q_0(1), q_1(1), rho_d, n_wl, wl(:), k_ext(1,:), ssa(1,:), g(1,:))

        !! increment time
        time = time + t_step

        !! Print to screen current progress
        print*, n, time, pl(1) * 1e-5_dp, 'bar ', Tl(1), 'K ', mu(1), 'g mol-1 ',  trim(sp)

        !! Mean mass of particle [g]
        m_c(1) = (q_1(1)*rho(1))/(q_0(1)*nd_atm(1))

        !! Mass weighted mean radius of particle [um]
        r_c(1) = ((3.0_dp*m_c(1))/(4.0_dp*pi*rho_d))**(1.0_dp/3.0_dp) * 1e4_dp

        print*, 'q', n, q_v(1), q_0(1), q_1(1), vf(1)
        print*, 'r', n, m_c(1), r_c(1), rho_d
        print*, 'o', n, k_ext(1,1), ssa(1,1), g(1,1), k_ext(1,n_wl), ssa(1,n_wl), g(1,n_wl)

        !! mini-cloud test output
        call output(n, time, nlay)

      end do

    case(2)

      !! In this example we perform a 1D test of mini_cloud with diffusion and settling included

      nlay = 100
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
      p_bot = 100.0_dp * 1e5_dp

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

      !! Assume constant H2, He and H background VMR @ approx solar
      allocate(VMR(nlay,2),sp_bg(2))
      sp_bg = (/'H2','He'/)
      VMR(:,1) = 0.85_dp
      VMR(:,2) = 0.15_dp

      !! Assumed condensate species
      sp = 'KCl'
      rho_d = 1.99_dp
      mol_w_sp = 74.5513_dp

      allocate(q_v(nlay), q_0(nlay), q_1(nlay), q0(3), q(nlay,3))
      allocate(r_c(nlay), m_c(nlay), vf(nlay), r_c_old(nlay), del(nlay))

      q_v(:) = 1e-30_dp
      q_0(:) = 1e-30_dp
      q_1(:) = 1e-30_dp

      q0(1) = 1.17e-7_dp * mol_w_sp/mu(nlay)
      q0(2) = 1e-30_dp
      q0(3) = 1e-30_dp

      q_v(nlay) = q0(1)

      time = 0.0_dp
      n = 0

      r_c_old(:) = r_seed * 1e-4_dp

   
      do n = 1, n_it

        do i = 1, nlay

          !! Call mini-cloud and perform integrations for a single layer
          call mini_cloud_2(Tl(i), pl(i), grav, mu(i), VMR(i,:), t_step, sp, sp_bg, q_v(i), q_0(i), q_1(i))

          !! Calculate settling velocity for this layer
          call mini_cloud_vf(Tl(i), pl(i), grav, mu(i), VMR(i,:), rho_d, sp_bg, q_0(i), q_1(i), vf(i))

          !! Calculate the opacity at the wavelength grid
         !call opac_mie(1, sp, Tl(i), mu(i), pl(i), q_0(i), q_1(i), rho_d, n_wl, wl, k_ext(i,:), ssa(i,:), g(i,:))
        end do

        q(:,1) = q_v(:)
        q(:,2) = q_0(:)
        q(:,3) = q_1(:)

        call vert_adv_exp_McCormack(nlay, nlev, t_step, mu, grav, Tl, pl, pe, vf, 2, q(:,2:3), q0(2:3))
        !call vert_adv_exp_MUSCL(nlay, nlev, t_step, mu, grav, Tl, pl, pe, vf, 2, q(:,2:3), q0(2:3))


        call vert_diff_exp(nlay, nlev, t_step, mu, grav, Tl, pl, pe, Kzz, 3, q(:,:), q0(:))
        !call vert_diff_imp(nlay, nlev, t_step, mu, grav, Tl, pl, pe, Kzz, 3, q(:,:), q0(:))



        q_v(:) = q(:,1)
        q_0(:) = q(:,2)
        q_1(:) = q(:,3)

        !! Number density [cm-3] of layer
        nd_atm(:) = (pl(:)*10.0_dp)/(kb*Tl(:))  

        !! Mass density of layer
        rho(:) = (pl(:)*10.0_dp*mu(:)*amu)/(kb * Tl(:)) ! Mass density [g cm-3]

        !! Mean mass of particle [g]
        m_c(:) = max((q_1(:)*rho)/(q_0(:)*nd_atm),m_seed)

        !! Mass weighted mean radius of particle [um]
        r_c(:) = max(((3.0_dp*m_c(:))/(4.0_dp*pi*rho_d))**(1.0_dp/3.0_dp),r_seed) * 1e4_dp

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

        if ((end .eqv. .True.) .and. (n > int(1e6))) then
          print*, 'exit: ', n, n_it, end
          print*, del(:)
          print*, del(:)/t_step
          exit
        end if


      end do

      print*, del(:)
      print*, del(:)/t_step

      do i = 1, nlay
        print*, i, pl(i)/1e5_dp, Tl(i), Kzz(i), q_v(i), q_0(i), q_1(i), r_c(i), vf(i), q_0(i)*nd_atm(i)
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
      open(newunit=u1,file='results_2mom_mono/tracers.txt',action='readwrite')
      open(newunit=u2,file='results_2mom_mono/opac.txt',action='readwrite')
      write(u2,*) wl(:)
      first_call = .False.
    end if

    do i = 1, nlay
      write(u1,*) t, time, Tl(i), pl(i), grav, mu(i), VMR(i,:), q_v(i), q_0(i), q_1(i), vf(i)
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

end program test_mini_cloud_2
