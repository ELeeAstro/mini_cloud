module vert_adv_exp_mod
  use, intrinsic :: iso_fortran_env
  implicit none
  integer, parameter :: dp = REAL64

  private :: koren_phi, adv_rhs_muscl_down_nugrid
  public  :: vert_adv_exp

contains

  subroutine vert_adv_exp(nlay, nlev, t_end, mu, grav_in, Tl, pl_in, pe_in, vf, nq, q)
    ! TVD MUSCL advection with SSPRK(3,3) time integration (3rd order in time).
    ! Assumes downward/positive velocities (upwind = cell above).
    implicit none
    integer, intent(in)            :: nlay, nlev, nq
    real(dp), intent(in)           :: t_end, grav_in
    real(dp), intent(in)           :: Tl(nlay), pl_in(nlay), mu(nlay), pe_in(nlev)
    real(dp), intent(in)           :: vf(nlay,nq)      ! settling speed (cm/s), per tracer
    real(dp), intent(inout)        :: q(nlay,nq)       ! cell-avg mixing ratio

    integer  :: i, k, n
    real(dp) :: R_gas, grav, eps, CFL
    real(dp) :: t, dt, dt_max, num, den, Te, mue
    real(dp), dimension(nq) :: q_top
    real(dp), allocatable :: alte(:), altm(:), dz(:), dzm(:)
    real(dp), allocatable :: pl(:), pe(:), rho(:), rho_e(:), v_e(:), A(:)
    real(dp), allocatable :: rhs(:), q1(:), q2(:)

    CFL   = 0.85_dp
    q_top = 0.0_dp              ! top inflow per tracer (set if needed)
    R_gas = 8.31446261815324e7_dp
    grav  = 100.0_dp * grav_in
    eps   = 1.0e-30_dp

    allocate(alte(nlev), altm(nlay), dz(nlay), dzm(nlay-1))
    allocate(pl(nlay), pe(nlev), rho(nlay), rho_e(nlev))
    allocate(v_e(nlev), A(nlev))
    allocate(rhs(nlay), q1(nlay), q2(nlay))

    ! ---- Geometry (hydrostatic heights in cm) ----
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
      dzm(i) = altm(i) - altm(i+1)          ! center-to-center spacing
    end do

    ! ---- Densities at centers and faces ----
    rho(:) = pl(:) / ((R_gas / mu(:)) * Tl(:))
    rho_e(1) = pe(1) / ((R_gas / mu(1)) * Tl(1)) 
    do k = 2, nlay
      ! face between cells k-1 and k corresponds to pe(k)
      Te = 0.5_dp*(Tl(k-1) + Tl(k))
      mue = 0.5_dp*(mu(k-1) + mu(k))
      rho_e(k) = pe(k) / ((R_gas / mue) * Te)
    end do
    rho_e(nlev) = pe(nlev) / ((R_gas / mu(nlay)) * Tl(nlay))  ! bottom face: one-sided

    ! ---- Loop tracers (vf may differ per tracer) ----
    do n = 1, nq

      ! Face velocities (downward/positive)
      v_e(1) = vf(1,n)
      do i = 1, nlay-1
        v_e(i+1) = 0.5_dp*(vf(i,n) + vf(i+1,n))
      end do
      v_e(nlev) = vf(nlay,n)

      ! Mass flux at faces A = rho_e * v_e
      do k = 1, nlev
        A(k) = rho_e(k) * v_e(k)
      end do

      ! Time step from upwind face Courant: dt <= CFL * min( rho_up*dz_up / A_face )
      dt_max = t_end
      do k = 2, nlay
        num    = rho(k-1) * dz(k-1)
        den    = A(k) + eps
        dt_max = min(dt_max, CFL * num / den)
      end do

      t = 0.0_dp
      do while (t < t_end)
        dt = merge(dt_max, t_end - t, t + dt_max <= t_end)

        ! ===== SSPRK(3,3) stages =====

        ! Stage 1: q1 = q^n + dt * L(q^n)
        call adv_rhs_muscl_down_nugrid(nlay, nlev, dz, dzm, rho, A, q(:,n), q_top(n), eps, rhs)
        do i = 1, nlay
          q1(i) = q(i,n) + dt * rhs(i)
          if (q1(i) < 1.0e-99_dp) q1(i) = 1.0e-99_dp
        end do

        ! Stage 2: q2 = 3/4 q^n + 1/4 ( q1 + dt * L(q1) )
        call adv_rhs_muscl_down_nugrid(nlay, nlev, dz, dzm, rho, A, q1, q_top(n), eps, rhs)
        do i = 1, nlay
          q2(i) = 0.75_dp*q(i,n) + 0.25_dp*( q1(i) + dt*rhs(i) )
          if (q2(i) < 1.0e-99_dp) q2(i) = 1.0e-99_dp
        end do

        ! Stage 3: q^{n+1} = 1/3 q^n + 2/3 ( q2 + dt * L(q2) )
        call adv_rhs_muscl_down_nugrid(nlay, nlev, dz, dzm, rho, A, q2, q_top(n), eps, rhs)
        do i = 1, nlay
          q(i,n) = (1.0_dp/3.0_dp)*q(i,n) + (2.0_dp/3.0_dp)*( q2(i) + dt*rhs(i) )
          if (q(i,n) < 1.0e-99_dp) q(i,n) = 1.0e-99_dp
        end do

        t = t + dt
      end do
    end do

    deallocate(alte, altm, dz, dzm, pl, pe, rho, rho_e, v_e, A, rhs, q1, q2)
  end subroutine vert_adv_exp

  ! ==== Non-uniform grid MUSCL RHS (downward flow A>0) ====
  subroutine adv_rhs_muscl_down_nugrid(nlay, nlev, dz, dzm, rho, A, q, q_top, eps, rhs)
    implicit none
    integer, intent(in) :: nlay, nlev
    real(dp), intent(in) :: dz(nlay), dzm(nlay-1), rho(nlay), A(nlev)
    real(dp), intent(in) :: q(nlay), q_top, eps
    real(dp), intent(out):: rhs(nlay)

    integer :: i, k
    real(dp) :: gL, gR, r, phi
    real(dp), allocatable :: sigma(:), qR_face(:), F(:)
    real(dp) :: dzhalf_top, gL_top, r1, phi1

    allocate(sigma(nlay), qR_face(nlay), F(nlev))

    ! ---- Limited slopes sigma_i ~ dq/dz (use downwind gradient for positive flow) ----
    sigma(:) = 0.0_dp
    do i = 2, nlay-1
      gL = (q(i)   - q(i-1)) / (dzm(i-1) + eps)     ! left one-sided gradient
      gR = (q(i+1) - q(i)  ) / (dzm(i)   + eps)     ! right one-sided gradient
      r  = gL / (gR + eps)
      phi = koren_phi(r)
      sigma(i) = phi * gR
    end do

    ! ---- Boundary-aware reconstruction at top (inflow) ----
    dzhalf_top = 0.5_dp*dz(1)                       ! center(1) to top face
    gL_top     = (q(1) - q_top) / (dzhalf_top + eps)
    gR         = (q(2) - q(1)) / (dzm(1) + eps)
    r1         = gL_top / (gR + eps)
    phi1       = koren_phi(r1)
    sigma(1)   = phi1 * gR

    ! ---- Right face values qR_face(i) = q_i + 0.5*dzm(i)*sigma_i ----
    do i = 1, nlay-1
      qR_face(i) = q(i) + 0.5_dp*dzm(i)*sigma(i)
    end do
    qR_face(nlay) = q(nlay)   ! unused (no interior face to the right)

    ! Clamp the first right-face value between local extrema (Barth–Jespersen style)
    qR_face(1) = max( min(qR_face(1), max(q_top, q(1), q(2))), &
                      min(q_top, q(1), q(2)) )

    ! ---- Upwind fluxes (A>0 ⇒ use cell (i-1) right face at face k=i) ----
    F(1) = A(1) * q_top
    do k = 2, nlay
      F(k) = A(k) * qR_face(k-1)
    end do
    F(nlev) = 0.0_dp   ! outer face (not used in divergence)

    ! ---- Conservative RHS: -div(F)/(rho*dz) ----
    do i = 1, nlay
      rhs(i) = - (F(i+1) - F(i)) / (rho(i)*dz(i) + eps)
    end do

    deallocate(sigma, qR_face, F)
  end subroutine adv_rhs_muscl_down_nugrid

  pure real(dp) function koren_phi(r) result(phi)
    real(dp), intent(in) :: r
    phi = max(0.0_dp, min(2.0_dp*r, (1.0_dp + 2.0_dp*r)/3.0_dp, 2.0_dp))
  end function koren_phi

end module vert_adv_exp_mod
