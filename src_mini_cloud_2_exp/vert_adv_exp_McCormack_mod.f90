module vert_adv_exp_McCormack_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: CFL = 0.95_dp
  real(dp), parameter :: R = 8.31446261815324e7_dp
  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: amu = 1.66053906892e-24_dp



  real(dp), parameter :: q_min = 1e-99_dp

  character(len=*), parameter :: limiter = 'minmod'

  public :: vert_adv_exp_McCormack
  private :: minmod, superbee, vanleer, mc, koren

  contains
 
  subroutine vert_adv_exp_McCormack(nlay, nlev, t_end, mu, grav_in, Tl, pl_in, pe_in, vf, nq, q_in, q0)
    implicit none


    integer, intent(in) :: nlay, nlev, nq
    real(dp), intent(in) :: t_end, grav_in
    real(dp), dimension(nlay), intent(in) :: Tl, pl_in, vf, mu
    real(dp), dimension(nlev), intent(in) :: pe_in
    real(dp), dimension(nq), intent(in) :: q0

    real(dp), dimension(nlay,nq), intent(inout) :: q_in

    real(dp), dimension(nlay,nq) :: q

    integer :: k
    real(dp) :: grav
    real(dp), dimension(nlev) :: alte, vf_e, pe
    real(dp), dimension(nlay) :: altm,  pl, rho
    real(dp), dimension(nlay-1) :: delz_mid, c

    real(dp), dimension(nlay,nq) :: qc
    real(dp), dimension(nlay) :: sig

    integer :: n_it, n
    real(dp) :: dt, t_now


    pl(:) = pl_in(:) * 10.0_dp   ! Convert pascal to dyne cm-2
    pe(:) = pe_in(:) * 10.0_dp   ! Convert pascal to dyne cm-2

    grav = grav_in * 100.0_dp


    ! Find velocity at levels
    vf_e(1) = vf(1)
    do k = 1, nlay-1
      vf_e(k+1) = (vf(k+1) + vf(k))/2.0_dp
    end do
    vf_e(nlev) = vf(nlay)

    !! First calculate the vertical height (cm) assuming hydrostatic equilibrium and differences
    alte(nlev) = 0.0_dp
    do k = nlev-1, 1, -1
      alte(k) = alte(k+1) + (R*Tl(k))/(mu(k)*grav) * log(pe(k+1)/pe(k))
    end do

    rho(:) = pl(:) / ((R / mu(:)) * Tl(:))

    do n = 1, nq
      q(:,n) = q_in(:,n) * rho(:)
    end do

    !! Find differences between layers directly
    do k = 1, nlay
      altm(k) = 0.5_dp * (alte(k) + alte(k+1))
    end do

    do k = 1, nlay-1
      delz_mid(k) = altm(k) - altm(k+1)
    end do

    !! Find minimum timestep that allows the CFL condition
    dt = t_end
    do k = 1, nlay-1
      dt = min(dt,CFL*(delz_mid(k)/abs(vf_e(k+1))))
    enddo

    t_now = 0.0_dp
    n_it = 0

    do while ((t_now < t_end) .and. (n_it < 1000000))

      ! If next time step overshoots - last time step is equal tend
      if (t_now + dt > t_end) then
        dt = t_end - t_now
      end if

      !! Find the courant number
      c(:) = (abs(vf_e(2:nlay)) * dt) / delz_mid(:)

      do n = 1, nq



        !! Find the minmod limiter
        select case (trim(adjustl(limiter)))
        case ('minmod')
          call minmod(nlay, q(:,n), delz_mid(:), sig)
        case ('superbee')
          call superbee(nlay, q(:,n), delz_mid(:), sig)
        case ('vanleer')
          call vanleer(nlay, q(:,n), delz_mid(:), sig)
        case ('mc')
          call mc(nlay, q(:,n), delz_mid(:), sig)
        case ('koren')
          call koren(nlay, q(:,n), delz_mid(:), sig)
        case default
          print *, 'Error: Unknown limiter: ', trim(limiter)
          print*, 'STOP'
          stop
        end select

        !! Perform McCormack step - do not evolve index nlay as we have fixed boundary condition
        !! Predictor step (forward in space)
        qc(:,n) = q(:,n)
        qc(:nlay-1,n) = q(:nlay-1,n) - sig(:nlay-1)*c(:nlay-1)*(q(2:nlay,n) - q(:nlay-1,n))
        
        !! Corrector step (backward in space)
        q(2:nlay-1,n) = 0.5_dp * (q(2:nlay-1,n) + qc(2:nlay-1,n) - c(2:nlay-1)*(qc(2:nlay-1,n) - qc(1:nlay-2,n)))
        
        q(:,n) = max(q(:,n),q_min)

      end do

      t_now = t_now + dt
      n_it = n_it + 1

      ! Apply boundary conditions
      q(nlay,:) = q0(:)
      q(1,:) = q(2,:)

    end do

    do n = 1, nq
      q_in(:,n) = q(:,n)/rho(:)
      q_in(:,n) = max(q_in(:,n),q_min)
    end do
    
  end subroutine vert_adv_exp_McCormack

  subroutine minmod(nlay, q, dz, sig)
    implicit none
  
    integer, intent(in) :: nlay
    real(dp), dimension(nlay), intent(in) :: q
    real(dp), dimension(nlay-1), intent(in) :: dz
    real(dp), dimension(nlay), intent(out) :: sig
  
    integer :: i
    real(dp), parameter :: eps = 1.0e-10_dp
    real(dp) :: r, dq_minus, dq_plus
  
    sig(1) = 0.0_dp
    do i = 2, nlay-1
      dq_minus = (q(i) - q(i-1)) / dz(i-1)
      dq_plus  = (q(i+1) - q(i)) / dz(i)
      r = dq_minus / (dq_plus + sign(eps, dq_plus))
      sig(i) = max(0.0_dp, min(1.0_dp, r))
    end do
    sig(nlay) = 0.0_dp
  
  end subroutine minmod

  subroutine superbee(nlay, q, dz, sig)
    implicit none
  
    integer, intent(in) :: nlay
    real(dp), dimension(nlay), intent(in) :: q
    real(dp), dimension(nlay-1), intent(in) :: dz 
    real(dp), dimension(nlay), intent(out) :: sig
  
    integer :: i
    real(dp), parameter :: eps = 1.0e-10_dp
    real(dp) :: r, dq_minus, dq_plus
    
    sig(1) = 0.0_dp
    do i = 2, nlay-1
      dq_minus = (q(i) - q(i-1)) / dz(i-1)
      dq_plus  = (q(i+1) - q(i)) / dz(i)
      r = dq_minus / (dq_plus + sign(eps, dq_plus))
      sig(i) = max(0.0_dp, min(2.0_dp*r, 1.0_dp), min(r, 2.0_dp))
    end do
    sig(nlay) = 0.0_dp
  
  end subroutine superbee

  subroutine vanleer(nlay, q, dz, sig)
    implicit none
  
    integer, intent(in) :: nlay
    real(dp), dimension(nlay), intent(in) :: q
    real(dp), dimension(nlay-1), intent(in) :: dz
    real(dp), dimension(nlay), intent(out) :: sig
  
    integer :: i
    real(dp), parameter :: eps = 1.0e-10_dp
    real(dp) :: dq_minus, dq_plus, r
  
    sig(1) = 0.0_dp
    do i = 2, nlay-1
      dq_minus = (q(i) - q(i-1)) / dz(i-1)
      dq_plus  = (q(i+1) - q(i)) / dz(i)
      r = dq_minus / (dq_plus + sign(eps, dq_plus))
      sig(i) = (r + abs(r)) / (1.0_dp + abs(r))
    end do
    sig(nlay) = 0.0_dp
  
  end subroutine vanleer

  subroutine mc(nlay, q, dz, sig)
    implicit none
  
    integer, intent(in) :: nlay
    real(dp), dimension(nlay), intent(in) :: q
    real(dp), dimension(nlay-1), intent(in) :: dz
    real(dp), dimension(nlay), intent(out) :: sig
  
    integer :: i
    real(dp), parameter :: eps = 1.0e-10_dp
    real(dp) :: dq_minus, dq_plus, r
  
    sig(1) = 0.0_dp
    do i = 2, nlay-1
      dq_minus = (q(i) - q(i - 1)) / dz(i - 1)
      dq_plus  = (q(i + 1) - q(i)) / dz(i)
      r = dq_minus / (dq_plus + sign(eps, dq_plus))
      sig(i) = max(0.0_dp, min((1.0_dp + r) / 2.0_dp, 2.0_dp, r))
    end do
    sig(nlay) = 0.0_dp
  
  end subroutine mc
  
  subroutine koren(nlay, q, dz, sig)
    implicit none
  
    integer, intent(in) :: nlay
    real(dp), dimension(nlay), intent(in) :: q
    real(dp), dimension(nlay-1), intent(in) :: dz
    real(dp), dimension(nlay), intent(out) :: sig
  
    integer :: i
    real(dp), parameter :: eps = 1.0e-10_dp
    real(dp) :: dq_minus, dq_plus, r
  
    sig(1) = 0.0_dp
    do i = 2, nlay-1
      dq_minus = (q(i) - q(i-1)) / dz(i - 1)
      dq_plus  = (q(i+1) - q(i)) / dz(i)
      r = dq_minus / (dq_plus + sign(eps, dq_plus))
      sig(i) = max(0.0_dp, min(2.0_dp*r, (1.0_dp + 2.0_dp*r) / 3.0_dp, 2.0_dp))
    end do
    sig(nlay) = 0.0_dp
  
  end subroutine koren
  
end module vert_adv_exp_McCormack_mod
