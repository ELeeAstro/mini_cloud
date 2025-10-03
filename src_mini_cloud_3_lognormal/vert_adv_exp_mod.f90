module vert_adv_exp_mod
  use, intrinsic :: iso_fortran_env
  implicit none
  integer, parameter :: dp = REAL64

  ! Numerics (match Python-translated version)
  real(dp), parameter :: TINY  = 1.0e-30_dp      ! denom guard
  real(dp), parameter :: QMIN  = 1.0e-99_dp      ! positivity floor
  real(dp), parameter :: CFL   = 0.85_dp         ! SSPRK(3,3)

  private
  public :: vert_adv_exp

contains

  subroutine vert_adv_exp(nlay, nlev, t_end, mu, grav_in, Tl, pl_in, pe_in, vf, nq, q)
    ! Multi-tracer advection with shared velocity vf(nlay) for all tracers.
    ! Downward/positive flow only (A >= 0). Shared dt for all tracers.
    !
    ! Inputs:
    !   nlay, nlev=nlay+1
    !   t_end [s]
    !   mu(nlay) [g/mol], grav_in [m/s^2], Tl(nlay) [K], pl_in(nlay) [Pa], pe_in(nlev) [Pa]
    !   vf(nlay) [cm/s]  -- identical for all tracers
    !   nq: number of tracers
    ! InOut:
    !   q(nlay,nq): mixing ratios (cell-avg)
    implicit none
    integer, intent(in)            :: nlay, nlev, nq
    real(dp), intent(in)           :: t_end, grav_in
    real(dp), intent(in)           :: Tl(nlay), pl_in(nlay), mu(nlay), pe_in(nlev)
    real(dp), intent(in)           :: vf(nlay)          ! shared for all tracers
    real(dp), intent(inout)        :: q(nlay, nq)

    ! geometry/state
    real(dp), allocatable :: alte(:), altm(:), dz(:), dzm(:)
    real(dp), allocatable :: pl(:), pe(:), rho(:), rho_e(:)

    ! shared face velocities/fluxes
    real(dp), allocatable :: v_e(:), A(:)

    ! time integration buffers
    real(dp), allocatable :: q1(:,:), q2(:,:), rhs(:)
    real(dp), allocatable :: sigma(:), qR_face(:), F(:)

    ! scalars
    integer :: i, k, n
    real(dp) :: R_gas, grav, t, dt, dt_max, remaining, num, den
    real(dp) :: q_top, T_TOL

    q_top = 0.0_dp
    R_gas = 8.31446261815324e7_dp      ! erg/mol/K
    grav  = 100.0_dp * grav_in         ! cm/s^2
    T_TOL = 100.0_dp*epsilon(1.0_dp) * max(1.0_dp, abs(t_end))

    ! ---- allocate ----
    allocate(alte(nlev), altm(nlay), dz(nlay), dzm(nlay-1))
    allocate(pl(nlay), pe(nlev), rho(nlay), rho_e(nlev))
    allocate(v_e(nlev), A(nlev))
    allocate(q1(nlay,nq), q2(nlay,nq), rhs(nlay))
    allocate(sigma(nlay), qR_face(nlay), F(nlev))

    ! ---- hydrostatic altitudes (Pa -> dyn/cm^2) ----
    pl = 10.0_dp * pl_in
    pe = 10.0_dp * pe_in

    alte(nlev) = 0.0_dp
    do k = nlev-1, 1, -1
      alte(k) = alte(k+1) + (R_gas*Tl(k))/(mu(k)*grav) * log(pe(k+1)/pe(k))
    end do
    do i = 1, nlay
      dz(i)   = alte(i) - alte(i+1)
      altm(i) = 0.5_dp*(alte(i) + alte(i+1))
    end do
    do i = 1, nlay-1
      dzm(i) = altm(i) - altm(i+1)
    end do

    ! ---- densities at centres and faces (Python behaviour) ----
    rho(:) = pl(:) / ((R_gas / mu(:)) * Tl(:))
    rho_e(1) = rho(1)
    do k = 2, nlay
      rho_e(k) = 0.5_dp*(rho(k-1) + rho(k))
    end do
    rho_e(nlev) = rho(nlay)

    ! ---- face velocities from vf (shared across tracers) ----
    v_e(1) = vf(1)
    do i = 1, nlay-1
      v_e(i+1) = 0.5_dp*(vf(i) + vf(i+1))
    end do
    v_e(nlev) = vf(nlay)

    ! ---- face mass flux ----
    do k = 1, nlev
      A(k) = rho_e(k) * v_e(k)
    end do

    ! ---- time integration loop (shared dt for all tracers) ----
    t = 0.0_dp
    do while (t < t_end - T_TOL)

      ! Python-style CFL: interior faces only (k=2..nlay)
      dt_max = t_end - t
      do k = 2, nlay
        num = rho(k-1) * dz(k-1)
        den = A(k) + TINY
        if (den > 0.0_dp) dt_max = min(dt_max, CFL * num / den)
      end do

      remaining = t_end - t
      dt = min(dt_max, remaining)
      if (dt <= 0.0_dp) exit

      ! ===== Stage 1: q1 = q + dt*L(q) for ALL tracers =====
      do n = 1, nq
        call adv_rhs_muscl_down_nugrid(nlay, nlev, dz, dzm, rho, A, q(:,n), q_top, &
                                       sigma, qR_face, F, rhs)
        do i = 1, nlay
          q1(i,n) = q(i,n) + dt*rhs(i)
          if (q1(i,n) < QMIN) q1(i,n) = QMIN
        end do
      end do

      ! ===== Stage 2: q2 = 3/4 q + 1/4 (q1 + dt*L(q1)) =====
      do n = 1, nq
        call adv_rhs_muscl_down_nugrid(nlay, nlev, dz, dzm, rho, A, q1(:,n), q_top, &
                                       sigma, qR_face, F, rhs)
        do i = 1, nlay
          q2(i,n) = 0.75_dp*q(i,n) + 0.25_dp*( q1(i,n) + dt*rhs(i) )
          if (q2(i,n) < QMIN) q2(i,n) = QMIN
        end do
      end do

      ! ===== Stage 3: q = 1/3 q + 2/3 (q2 + dt*L(q2)) =====
      do n = 1, nq
        call adv_rhs_muscl_down_nugrid(nlay, nlev, dz, dzm, rho, A, q2(:,n), q_top, &
                                       sigma, qR_face, F, rhs)
        do i = 1, nlay
          q(i,n) = (1.0_dp/3.0_dp)*q(i,n) + (2.0_dp/3.0_dp)*( q2(i,n) + dt*rhs(i) )
          if (q(i,n) < QMIN) q(i,n) = QMIN
        end do
      end do

      t = t + dt
    end do

    ! ---- cleanup ----
    deallocate(alte, altm, dz, dzm, pl, pe, rho, rho_e, v_e, A)
    deallocate(q1, q2, rhs, sigma, qR_face, F)
  end subroutine vert_adv_exp


  ! ----- RHS: non-uniform MUSCL (downward flow, A >= 0), matches Python -----
  subroutine adv_rhs_muscl_down_nugrid(nlay, nlev, dz, dzm, rho, A, q, q_top, &
                                       sigma, qR_face, F, rhs)
    implicit none
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: dz(nlay), dzm(nlay-1), rho(nlay), A(nlev)
    real(dp), intent(in) :: q(nlay), q_top
    real(dp), intent(inout) :: sigma(nlay), qR_face(nlay), F(nlev)
    real(dp), intent(out) :: rhs(nlay)

    integer :: i, k
    real(dp) :: gL, gR, r, dzhalf_top, gL_top, r1
    real(dp) :: mmax, mmin

    sigma(:)   = 0.0_dp
    qR_face(:) = 0.0_dp
    F(:)       = 0.0_dp

    ! interior slopes
    do i = 2, nlay-1
      gL = (q(i)   - q(i-1)) / (dzm(i-1) + TINY)
      gR = (q(i+1) - q(i)  ) / (dzm(i)   + TINY)
      r  = gL / (gR + TINY)
      sigma(i) = koren_phi(r) * gR
    end do

    ! top boundary-aware slope
    dzhalf_top = 0.5_dp*dz(1)
    gL_top     = (q(1) - q_top) / (dzhalf_top + TINY)
    if (nlay > 1) then
      gR  = (q(2) - q(1)) / (dzm(1) + TINY)
      r1  = gL_top / (gR + TINY)
      sigma(1) = koren_phi(r1) * gR
    else
      sigma(1) = gL_top
    end if

    ! right-face reconstructions
    do i = 1, nlay-1
      qR_face(i) = q(i) + 0.5_dp*dzm(i)*sigma(i)
    end do

    ! clamp first interior face between local extrema (Barth–Jespersen style)
    if (nlay > 1) then
      mmax = q_top; if (q(1) > mmax) mmax = q(1); if (q(2) > mmax) mmax = q(2)
      mmin = q_top; if (q(1) < mmin) mmin = q(1); if (q(2) < mmin) mmin = q(2)
      if (qR_face(1) > mmax) qR_face(1) = mmax
      if (qR_face(1) < mmin) qR_face(1) = mmin
    end if

    ! upwind fluxes for A >= 0
    F(1) = A(1) * q_top
    do k = 2, nlay
      F(k) = A(k) * qR_face(k-1)
    end do
    F(nlev) = 0.0_dp

    ! conservative RHS
    do i = 1, nlay
      rhs(i) = - (F(i+1) - F(i)) / (rho(i)*dz(i) + TINY)
    end do
  end subroutine adv_rhs_muscl_down_nugrid


  pure real(dp) function koren_phi(r) result(phi)
    implicit none
    real(dp), intent(in) :: r
    real(dp) :: a, b, m1, m2
    a  = 2.0_dp * r
    b  = (1.0_dp + 2.0_dp * r) / 3.0_dp
    m1 = merge(a, b, a < b)             ! min(a,b)
    m2 = merge(2.0_dp, m1, 2.0_dp < m1) ! min(2,m1)
    if (m2 > 0.0_dp) then
      phi = m2
    else
      phi = 0.0_dp
    end if
  end function koren_phi

end module vert_adv_exp_mod

