module vert_adv_imp_mod
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: dp = REAL64

  real(dp), parameter :: QMIN  = 1.0e-30_dp
  real(dp), parameter :: EPS   = 1.0e-300_dp
  real(dp), parameter :: R_gas = 8.31446261815324e7_dp
  real(dp), parameter :: GAMMA = 2.0_dp - sqrt(2.0_dp)
  real(dp), parameter :: BDF_STAGE = 1.0_dp/(GAMMA*(2.0_dp - GAMMA))
  real(dp), parameter :: BDF_OLD = (1.0_dp - GAMMA)**2/(GAMMA*(2.0_dp - GAMMA))

  private
  public :: vert_adv_imp

contains

  subroutine vert_adv_imp(nlay, nlev, t_end, mu, grav_in, Tl, pl_in, pe_in, vf, nq, q)
    implicit none

    integer, intent(in)     :: nlay, nlev, nq
    real(dp), intent(in)    :: t_end, grav_in
    real(dp), intent(in)    :: Tl(nlay), pl_in(nlay), mu(nlay), pe_in(nlev)
    real(dp), intent(in)    :: vf(nlay, nq)
    real(dp), intent(inout) :: q(nlay, nq)

    real(dp) :: pl(nlay), pe(nlev)
    real(dp) :: alte(nlev), dz(nlay)
    real(dp) :: rho(nlay), rho_e(nlev)
    real(dp) :: v_e(nlev), A(nlev), lambda(nlay)
    real(dp) :: q_old(nlay), q_stage(nlay), rhs_vec(nlay)
    real(dp) :: grav, q_top
    integer :: i, k, n

    if (nlev /= nlay + 1) error stop 'vert_adv_imp: nlev must equal nlay + 1'

    if (t_end <= 0.0_dp) return

    q_top = 0.0_dp
    grav = 100.0_dp * grav_in

    pl = 10.0_dp * pl_in
    pe = 10.0_dp * pe_in

    alte(nlev) = 0.0_dp
    do k = nlev-1, 1, -1
      alte(k) = alte(k+1) + (R_gas*Tl(k))/(mu(k)*grav) * log(pe(k+1)/pe(k))
    end do

    do i = 1, nlay
      dz(i) = alte(i) - alte(i+1)
    end do

    rho(:) = pl(:) / ((R_gas / mu(:)) * Tl(:))
    rho_e(1) = rho(1)
    do k = 2, nlay
      rho_e(k) = 0.5_dp*(rho(k-1) + rho(k))
    end do
    rho_e(nlev) = rho(nlay)

    do i = 1, nlay
      lambda(i) = 0.5_dp*GAMMA*t_end / (rho(i)*dz(i) + EPS)
    end do

    do n = 1, nq
      v_e(1) = max(0.0_dp, vf(1,n))
      do i = 1, nlay-1
        v_e(i+1) = max(0.0_dp, 0.5_dp*(vf(i,n) + vf(i+1,n)))
      end do
      v_e(nlev) = max(0.0_dp, vf(nlay,n))

      do k = 1, nlev
        A(k) = rho_e(k) * v_e(k)
      end do

      q_old(:) = q(:,n)

      ! TR-BDF2 stage 1:
      ! (I - gamma*dt/2 L) q_stage = q_old + gamma*dt/2 L(q_old)
      call build_tr_rhs(nlay, A, lambda, q_top, q_old, rhs_vec)
      call solve_implicit_upwind(nlay, A, lambda, q_top, rhs_vec, q_stage)

      ! TR-BDF2 stage 2:
      ! (I - gamma*dt/2 L) q_new = c_stage*q_stage - c_old*q_old
      rhs_vec(:) = BDF_STAGE*q_stage(:) - BDF_OLD*q_old(:)
      call solve_implicit_upwind(nlay, A, lambda, q_top, rhs_vec, q(:,n))
    end do

  end subroutine vert_adv_imp

  subroutine build_tr_rhs(nlay, A, lambda, q_top, q_in, rhs)
    implicit none

    integer, intent(in) :: nlay
    real(dp), intent(in) :: A(nlay+1), lambda(nlay), q_top, q_in(nlay)
    real(dp), intent(out) :: rhs(nlay)

    integer :: i

    rhs(1) = q_in(1) + lambda(1)*(A(1)*q_top - A(2)*q_in(1))
    do i = 2, nlay
      rhs(i) = q_in(i) + lambda(i)*(A(i)*q_in(i-1) - A(i+1)*q_in(i))
    end do

  end subroutine build_tr_rhs

  subroutine solve_implicit_upwind(nlay, A, lambda, q_top, rhs, q_out)
    implicit none

    integer, intent(in) :: nlay
    real(dp), intent(in) :: A(nlay+1), lambda(nlay), q_top, rhs(nlay)
    real(dp), intent(out) :: q_out(nlay)

    integer :: i

    q_out(1) = (rhs(1) + lambda(1)*A(1)*q_top) / (1.0_dp + lambda(1)*A(2))
    if (q_out(1) < QMIN) q_out(1) = QMIN

    do i = 2, nlay
      q_out(i) = (rhs(i) + lambda(i)*A(i)*q_out(i-1)) / (1.0_dp + lambda(i)*A(i+1))
      if (q_out(i) < QMIN) q_out(i) = QMIN
    end do

  end subroutine solve_implicit_upwind

end module vert_adv_imp_mod
